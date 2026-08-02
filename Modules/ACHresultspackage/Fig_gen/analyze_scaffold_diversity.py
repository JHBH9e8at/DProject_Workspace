import argparse
import json
import math
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


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
    cfg["fp_radius"] = int(cfg.get("fp_radius", "2"))
    cfg["fp_size"] = 16384
    cfg["diversity_sample_size"] = int(
        cfg.get("diversity_sample_size", "3000")
    )
    cfg["novelty_reference_size"] = int(
        cfg.get("novelty_reference_size", "5000")
    )
    cfg["random_state"] = int(cfg.get("random_state", "42"))
    return cfg


def resolve(value, config_path, source=None):
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = Path(config_path).resolve().parent / path
    if path.is_dir() and source:
        path = path / f"{source}_scores.csv"
    return path


def load_population(raw_path, source, cache_path):
    raw = pd.read_csv(raw_path, usecols=lambda c: c in {"step", "smiles"})
    membership = pd.read_csv(
        cache_path, usecols=["source", "source_index", "canon_smiles"]
    )
    membership = membership[membership["source"] == source].copy()
    membership["source_index"] = membership["source_index"].astype(int)
    df = raw.iloc[membership["source_index"].to_numpy()].copy()
    df["source_index"] = membership["source_index"].to_numpy()
    df["canon_smiles"] = membership["canon_smiles"].to_numpy()
    df["source"] = source
    df["step"] = pd.to_numeric(df["step"], errors="coerce")
    return df.dropna(subset=["step", "canon_smiles"]).reset_index(drop=True)


def murcko_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return "[NO_SCAFFOLD]"
    return Chem.MolToSmiles(scaffold, canonical=True)


def pairwise_similarity(fps):
    if len(fps) < 2:
        return {"pairs": 0, "mean_similarity": None, "median_similarity": None,
                "q95_similarity": None, "diversity": None}
    similarities = []
    for index in range(1, len(fps)):
        similarities.extend(
            DataStructs.BulkTanimotoSimilarity(fps[index], fps[:index])
        )
    values = np.asarray(similarities, dtype=float)
    return {
        "pairs": len(values),
        "mean_similarity": float(values.mean()),
        "median_similarity": float(np.median(values)),
        "q95_similarity": float(np.quantile(values, 0.95)),
        "diversity": float(1.0 - values.mean()),
    }


def scaffold_metrics(scaffolds):
    counts = Counter(scaffolds)
    total = sum(counts.values())
    unique = len(counts)
    probabilities = np.asarray(list(counts.values()), dtype=float) / total
    entropy = float(-(probabilities * np.log(probabilities)).sum())
    normalized = entropy / math.log(unique) if unique > 1 else 0.0
    ordered = sorted(counts.values(), reverse=True)
    return {
        "molecules": total,
        "unique_scaffolds": unique,
        "unique_scaffold_fraction": unique / total if total else None,
        "singleton_scaffolds": sum(value == 1 for value in counts.values()),
        "singleton_scaffold_fraction": (
            sum(value == 1 for value in counts.values()) / unique
            if unique else None
        ),
        "top1_scaffold_share": ordered[0] / total if total else None,
        "top5_scaffold_share": sum(ordered[:5]) / total if total else None,
        "shannon_entropy": entropy,
        "normalized_entropy": normalized,
    }


def sample_frame(frame, size, rng):
    if len(frame) <= size:
        return frame
    indices = rng.choice(frame.index.to_numpy(), size=size, replace=False)
    return frame.loc[indices]


def analyze_source(df, source, cfg):
    rng = np.random.default_rng(cfg["random_state"])
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["fp_radius"], fpSize=cfg["fp_size"]
    )
    df = df.copy()
    df["mol"] = df["canon_smiles"].map(Chem.MolFromSmiles)
    df = df[df["mol"].notna()].copy()
    df["fp"] = df["mol"].map(generator.GetFingerprint)
    df["murcko_scaffold"] = df["canon_smiles"].map(murcko_smiles)
    df = df[df["murcko_scaffold"].notna()].copy()

    max_step = int(df["step"].max())
    first_last_step = max(int(df["step"].min()), max_step - cfg["last_n"] + 1)
    recent = df[df["step"].between(first_last_step, max_step)].copy()

    fp_step_rows = []
    scaffold_step_rows = []
    seen_scaffolds = set()
    previous = pd.DataFrame()
    for step in sorted(df["step"].unique()):
        current = df[df["step"] == step]
        fp_metrics = pairwise_similarity(current["fp"].tolist())

        novelty_mean = None
        novelty_median = None
        fraction_below_04 = None
        fraction_below_06 = None
        if not previous.empty:
            reference = sample_frame(
                previous, cfg["novelty_reference_size"], rng
            )
            reference_fps = reference["fp"].tolist()
            nearest = np.asarray([
                max(DataStructs.BulkTanimotoSimilarity(fp, reference_fps))
                for fp in current["fp"]
            ])
            novelty_mean = float(nearest.mean())
            novelty_median = float(np.median(nearest))
            fraction_below_04 = float((nearest < 0.4).mean())
            fraction_below_06 = float((nearest < 0.6).mean())

        fp_step_rows.append({
            "source": source, "step": int(step), "molecules": len(current),
            **fp_metrics,
            "mean_nearest_similarity_to_prior": novelty_mean,
            "median_nearest_similarity_to_prior": novelty_median,
            "fraction_nearest_below_0_4": fraction_below_04,
            "fraction_nearest_below_0_6": fraction_below_06,
        })

        scaffolds = current["murcko_scaffold"].tolist()
        metrics = scaffold_metrics(scaffolds)
        unique_current = set(scaffolds)
        new_scaffolds = unique_current - seen_scaffolds
        scaffold_step_rows.append({
            "source": source, "step": int(step), **metrics,
            "new_scaffolds": len(new_scaffolds),
            "new_scaffold_fraction": (
                len(new_scaffolds) / len(unique_current)
                if unique_current else None
            ),
        })
        seen_scaffolds.update(unique_current)
        previous = pd.concat([previous, current], ignore_index=True)

    aggregate_rows = []
    for scope, frame in (("full", df), ("last_window", recent)):
        sampled = sample_frame(frame, cfg["diversity_sample_size"], rng)
        aggregate_rows.append({
            "source": source,
            "scope": scope,
            "first_step": int(frame["step"].min()),
            "last_step": int(frame["step"].max()),
            "population_molecules": len(frame),
            "sampled_molecules": len(sampled),
            **pairwise_similarity(sampled["fp"].tolist()),
            **scaffold_metrics(frame["murcko_scaffold"].tolist()),
        })

    top_scaffold_rows = []
    counts = Counter(recent["murcko_scaffold"])
    for rank, (scaffold, count) in enumerate(counts.most_common(25), 1):
        top_scaffold_rows.append({
            "source": source, "rank": rank, "murcko_scaffold": scaffold,
            "count": count, "share": count / len(recent),
        })

    return {
        "fp_by_step": pd.DataFrame(fp_step_rows),
        "scaffold_by_step": pd.DataFrame(scaffold_step_rows),
        "aggregate": pd.DataFrame(aggregate_rows),
        "top_scaffolds_last_window": pd.DataFrame(top_scaffold_rows),
    }


def main(config_path):
    cfg = parse_config(config_path)
    outdir = resolve(cfg["outdir"], config_path)
    cache_path = outdir / "embedding_cache" / "joint_umap_fp.csv.gz"
    if not cache_path.exists():
        raise FileNotFoundError(f"Population cache not found: {cache_path}")

    combined = {}
    for source, key in (("PR", "prdir"), ("PPS", "ppsdir")):
        population = load_population(
            resolve(cfg[key], config_path, source), source, cache_path
        )
        print(f"[{source}] analyzing {len(population)} molecules")
        results = analyze_source(population, source, cfg)
        for name, table in results.items():
            combined.setdefault(name, []).append(table)

    tables = {
        name: pd.concat(parts, ignore_index=True)
        for name, parts in combined.items()
    }
    analysis_dir = outdir / "scaffold_diversity_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(analysis_dir / f"{name}.csv", index=False)
        print(f"[saved] {name}: {len(table)} rows")
    (analysis_dir / "manifest.json").write_text(
        json.dumps({
            "fingerprint": {
                "type": "Morgan count-independent bit fingerprint",
                "radius": cfg["fp_radius"],
                "fp_size": cfg["fp_size"],
                "similarity": "Tanimoto",
            },
            "scaffold": "Bemis-Murcko scaffold",
            "last_n": cfg["last_n"],
            "diversity_sample_size": cfg["diversity_sample_size"],
            "novelty_reference_size": cfg["novelty_reference_size"],
            "random_state": cfg["random_state"],
            "tables": {name: len(table) for name, table in tables.items()},
        }, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Analyze PR/PPS Morgan-fingerprint Tanimoto and Bemis-Murcko "
            "scaffold diversity over training and the final N steps."
        )
    )
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)
