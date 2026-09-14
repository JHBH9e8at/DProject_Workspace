import argparse
import json
import os
from pathlib import Path

import pandas as pd


DESCRIPTORS = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_TPSA", "desc_NumHAcceptors", "desc_NumHDonors",
    "desc_NumRotatableBonds", "desc_NumAromaticRings",
    "desc_NumAliphaticRings", "desc_RingCount",
    "desc_HeavyAtomCount", "desc_NumHeteroatoms", "desc_Bertz",
]
METRICS = ["docking_score", "ligand_efficiency"]


def parse_config(path):
    cfg = {}
    with open(path, encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(f"{path}:{lineno}: expected key=value")
            key, value = line.split("=", 1)
            cfg[key.strip().lower().replace("-", "_")] = value.strip()
    for key in ("prdir", "ppsdir", "outdir"):
        if not cfg.get(key):
            raise ValueError(f"Missing required config key: {key}")
    cfg["last_n"] = int(cfg.get("last_n", "50"))
    if cfg["last_n"] < 1:
        raise ValueError("last_n must be at least 1")
    return cfg


def resolve(value, config_path, source=None):
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = Path(config_path).resolve().parent / path
    if path.is_dir() and source:
        path = path / f"{source}_scores.csv"
    return path


def describe(values):
    series = pd.to_numeric(values, errors="coerce").dropna()
    quantiles = series.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    return {
        "count": len(series),
        "mean": series.mean(),
        "std": series.std(),
        "min": series.min(),
        "q05": quantiles.loc[0.05],
        "q25": quantiles.loc[0.25],
        "median": quantiles.loc[0.5],
        "q75": quantiles.loc[0.75],
        "q95": quantiles.loc[0.95],
        "max": series.max(),
        "iqr": quantiles.loc[0.75] - quantiles.loc[0.25],
        "skew": series.skew(),
        "kurtosis": series.kurt(),
    }


def load_population(raw_path, source, cache_path):
    raw = pd.read_csv(raw_path)
    membership = pd.read_csv(
        cache_path, usecols=["source", "source_index", "canon_smiles"]
    )
    membership = membership[membership["source"] == source].copy()
    membership["source_index"] = membership["source_index"].astype(int)
    df = raw.iloc[membership["source_index"].to_numpy()].copy()
    df["source_index"] = membership["source_index"].to_numpy()
    df["canon_smiles"] = membership["canon_smiles"].to_numpy()
    df["source"] = source
    df["docking_score"] = pd.to_numeric(
        df[f"{source}_r_i_docking_score"], errors="coerce"
    )
    le_col = f"{source}_r_i_glide_ligand_efficiency"
    if le_col not in df:
        raise ValueError(f"{source}: missing ligand-efficiency column {le_col}")
    df["ligand_efficiency"] = pd.to_numeric(df[le_col], errors="coerce")
    df["step"] = pd.to_numeric(df["step"], errors="coerce")
    return df


def correlation(x, y, method):
    pair = pd.DataFrame({"x": x, "y": y}).apply(
        pd.to_numeric, errors="coerce"
    ).dropna()
    if len(pair) < 3:
        return len(pair), None
    if method == "spearman":
        value = pair["x"].rank().corr(pair["y"].rank(), method="pearson")
    else:
        value = pair["x"].corr(pair["y"], method="pearson")
    return len(pair), value


def make_tables(populations, last_n):
    overview_rows = []
    distribution_rows = []
    target_corr_rows = []
    matrix_rows = []
    step_rows = []

    for source, full in populations.items():
        max_step = int(full["step"].max())
        first_step = max(int(full["step"].min()), max_step - last_n + 1)
        recent = full[full["step"].between(first_step, max_step)].copy()
        overview_rows.append({
            "source": source,
            "full_population_rows": len(full),
            "last_n": last_n,
            "first_last_step": first_step,
            "max_step": max_step,
            "last_window_rows": len(recent),
            "distinct_last_steps": recent["step"].nunique(),
        })

        for scope, frame in (("full", full), ("last_window", recent)):
            for metric in METRICS:
                distribution_rows.append({
                    "source": source, "scope": scope, "metric": metric,
                    **describe(frame[metric]),
                })

        variables = METRICS + DESCRIPTORS
        for target in METRICS:
            for variable in variables:
                if variable == target:
                    continue
                n_p, pearson = correlation(
                    recent[target], recent[variable], "pearson"
                )
                n_s, spearman = correlation(
                    recent[target], recent[variable], "spearman"
                )
                target_corr_rows.append({
                    "source": source,
                    "target": target,
                    "variable": variable,
                    "n": min(n_p, n_s),
                    "pearson_r": pearson,
                    "spearman_rho": spearman,
                })

        for row_variable in variables:
            row = {"source": source, "variable": row_variable}
            for column_variable in variables:
                _, value = correlation(
                    recent[row_variable], recent[column_variable], "spearman"
                )
                row[column_variable] = value
            matrix_rows.append(row)

        grouped = recent.groupby("step")
        for step, frame in grouped:
            step_rows.append({
                "source": source,
                "step": step,
                "count": len(frame),
                "docking_mean": frame["docking_score"].mean(),
                "docking_median": frame["docking_score"].median(),
                "ligand_efficiency_mean": frame["ligand_efficiency"].mean(),
                "ligand_efficiency_median": frame["ligand_efficiency"].median(),
            })

    distribution = pd.DataFrame(distribution_rows)
    comparison = distribution.pivot(
        index=["scope", "metric"], columns="source",
        values=["mean", "median", "std", "q25", "q75"],
    )
    comparison.columns = [
        f"{stat}_{source}" for stat, source in comparison.columns
    ]
    comparison = comparison.reset_index()
    comparison["mean_PR_minus_PPS"] = (
        comparison["mean_PR"] - comparison["mean_PPS"]
    )
    comparison["median_PR_minus_PPS"] = (
        comparison["median_PR"] - comparison["median_PPS"]
    )

    return {
        "last50_overview": pd.DataFrame(overview_rows),
        "metric_distribution": distribution,
        "metric_comparison": comparison,
        "last50_target_correlations": pd.DataFrame(target_corr_rows),
        "last50_spearman_matrix": pd.DataFrame(matrix_rows),
        "last50_by_step": pd.DataFrame(step_rows),
    }


def main(config_path):
    cfg = parse_config(config_path)
    outdir = resolve(cfg["outdir"], config_path)
    cache_path = outdir / "embedding_cache" / "joint_umap_fp.csv.gz"
    if not cache_path.exists():
        raise FileNotFoundError(f"Population cache not found: {cache_path}")

    populations = {}
    for source, key in (("PR", "prdir"), ("PPS", "ppsdir")):
        populations[source] = load_population(
            resolve(cfg[key], config_path, source), source, cache_path
        )

    tables = make_tables(populations, cfg["last_n"])
    analysis_dir = outdir / "ligand_efficiency_last50_analysis"
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
    (analysis_dir / "manifest.json").write_text(
        json.dumps({
            "population": (
                "Same canonical-deduplicated, nonzero-docking population as "
                "joint_umap_fp cache."
            ),
            "last_n": cfg["last_n"],
            "ligand_efficiency_columns": {
                "PR": "PR_r_i_glide_ligand_efficiency",
                "PPS": "PPS_r_i_glide_ligand_efficiency",
            },
            "tables": {name: len(table) for name, table in tables.items()},
        }, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compare PR/PPS ligand efficiency and analyze last-N-step "
            "correlations among docking score, ligand efficiency, and descriptors."
        )
    )
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)
