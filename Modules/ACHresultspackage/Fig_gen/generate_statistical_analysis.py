import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


IMPORTANT_DESCRIPTORS = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_TPSA", "desc_NumHAcceptors", "desc_NumHDonors",
    "desc_NumRotatableBonds", "desc_NumAromaticRings",
    "desc_NumAliphaticRings", "desc_RingCount",
    "desc_HeavyAtomCount", "desc_NumHeteroatoms", "desc_Bertz",
]
SUMMARY_STATS = ["count", "mean", "std", "min", "q05", "q25", "median",
                 "q75", "q95", "max", "iqr", "skew", "kurtosis"]
CUTOFFS = [-5, -6, -7, -8, -9, -10, -11]


def parse_config(path):
    cfg = {}
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                cfg[key.strip().lower().replace("-", "_")] = value.strip()
    for key in ("prdir", "ppsdir", "outdir"):
        if not cfg.get(key):
            raise ValueError(f"Missing required config key: {key}")
    return cfg


def resolve(value, config_path, source=None):
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = Path(config_path).resolve().parent / path
    if path.is_dir() and source:
        path = path / f"{source}_scores.csv"
    return path


def describe_series(values):
    s = pd.to_numeric(values, errors="coerce").dropna()
    q = s.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    return {
        "count": len(s), "mean": s.mean(), "std": s.std(),
        "min": s.min(), "q05": q.loc[0.05], "q25": q.loc[0.25],
        "median": q.loc[0.5], "q75": q.loc[0.75], "q95": q.loc[0.95],
        "max": s.max(), "iqr": q.loc[0.75] - q.loc[0.25],
        "skew": s.skew(), "kurtosis": s.kurt(),
    }


def load_analysis_population(raw_path, source, cache_path):
    raw = pd.read_csv(raw_path)
    membership = pd.read_csv(
        cache_path, usecols=["source", "source_index", "canon_smiles"]
    )
    membership = membership[membership["source"] == source].copy()
    membership["source_index"] = membership["source_index"].astype(int)
    selected = raw.iloc[membership["source_index"].to_numpy()].copy()
    selected["source_index"] = membership["source_index"].to_numpy()
    selected["canon_smiles"] = membership["canon_smiles"].to_numpy()
    selected["source"] = source
    selected["docking_score"] = pd.to_numeric(
        selected[f"{source}_r_i_docking_score"], errors="coerce"
    )
    return raw, selected


def make_tables(raw_by_source, selected_by_source):
    overview = []
    distribution = []
    cutoff_rows = []
    descriptor_rows = []
    correlation_rows = []
    step_rows = []

    for source in ("PR", "PPS"):
        raw = raw_by_source[source]
        df = selected_by_source[source]
        score = df["docking_score"]
        overview.append({
            "source": source, "raw_rows": len(raw),
            "valid_rows": int(raw["valid"].astype(str).str.lower().eq("true").sum()),
            "unique_rows": int(raw["unique"].astype(str).str.lower().eq("true").sum()),
            "nonzero_docking_rows": int(pd.to_numeric(
                raw[f"{source}_r_i_docking_score"], errors="coerce"
            ).ne(0).sum()),
            "analysis_rows": len(df),
            "distinct_steps": df["step"].nunique(),
            "min_step": df["step"].min(), "max_step": df["step"].max(),
        })

        dist = {"source": source, **describe_series(score)}
        distribution.append(dist)
        for cutoff in CUTOFFS:
            count = int(score.le(cutoff).sum())
            cutoff_rows.append({
                "source": source, "cutoff": cutoff, "count": count,
                "fraction": count / len(score) if len(score) else np.nan,
            })

        for descriptor in IMPORTANT_DESCRIPTORS:
            stats = describe_series(df[descriptor])
            descriptor_rows.append({
                "source": source, "descriptor": descriptor, **stats
            })
            pair = df[["docking_score", descriptor]].apply(
                pd.to_numeric, errors="coerce"
            ).dropna()
            correlation_rows.append({
                "source": source,
                "descriptor": descriptor,
                "n": len(pair),
                "pearson_r": pair["docking_score"].corr(
                    pair[descriptor], method="pearson"
                ),
                "spearman_rho": pair["docking_score"].rank().corr(
                    pair[descriptor].rank(), method="pearson"
                ),
            })

        grouped = df.groupby("step")["docking_score"]
        for step, values in grouped:
            step_rows.append({
                "source": source, "step": step, "count": values.count(),
                "mean": values.mean(), "std": values.std(),
                "median": values.median(), "min": values.min(),
                "q25": values.quantile(0.25), "q75": values.quantile(0.75),
            })

    descriptor_df = pd.DataFrame(descriptor_rows)
    comparison = descriptor_df.pivot(
        index="descriptor", columns="source", values=["mean", "median", "std"]
    )
    comparison.columns = [f"{stat}_{source}" for stat, source in comparison.columns]
    comparison = comparison.reset_index()
    comparison["mean_PR_minus_PPS"] = comparison["mean_PR"] - comparison["mean_PPS"]
    comparison["median_PR_minus_PPS"] = (
        comparison["median_PR"] - comparison["median_PPS"]
    )

    return {
        "overview": pd.DataFrame(overview),
        "docking_distribution": pd.DataFrame(distribution),
        "docking_cutoffs": pd.DataFrame(cutoff_rows),
        "descriptor_summary": descriptor_df,
        "score_descriptor_correlation": pd.DataFrame(correlation_rows),
        "descriptor_comparison": comparison,
        "docking_by_step": pd.DataFrame(step_rows),
    }


def main(config_path):
    cfg = parse_config(config_path)
    outdir = resolve(cfg["outdir"], config_path)
    cache_path = outdir / "embedding_cache" / "joint_umap_fp.csv.gz"
    if not cache_path.exists():
        raise FileNotFoundError(
            f"Required population cache not found: {cache_path}"
        )
    raw_by_source = {}
    selected_by_source = {}
    for source, key in (("PR", "prdir"), ("PPS", "ppsdir")):
        raw_path = resolve(cfg[key], config_path, source)
        raw, selected = load_analysis_population(raw_path, source, cache_path)
        raw_by_source[source] = raw
        selected_by_source[source] = selected

    tables = make_tables(raw_by_source, selected_by_source)
    analysis_dir = outdir / "statistical_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(analysis_dir / f"{name}.csv", index=False)
        print(f"[saved] {name}: {len(table)} rows")
    json_tables = {
        name: {
            "columns": table.columns.tolist(),
            "rows": json.loads(table.to_json(orient="records")),
        }
        for name, table in tables.items()
    }
    (analysis_dir / "tables.json").write_text(
        json.dumps(json_tables, indent=2), encoding="utf-8"
    )
    manifest = {
        "population_definition": (
            "Same population as joint_umap_fp cache: valid, nonzero docking "
            "score, canonicalized, best score per canonical SMILES within source."
        ),
        "tables": {name: len(table) for name, table in tables.items()},
    }
    (analysis_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate comprehensive PR/PPS statistical tables."
    )
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)
