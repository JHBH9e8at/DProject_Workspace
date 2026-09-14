"""Morgan-fingerprint and Bemis-Murcko scaffold diversity analysis."""

from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
SCRIPT_ID = "PY-109"


def load_population(raw_path, source: str, membership_path) -> pd.DataFrame:
    raw = pd.read_csv(raw_path, usecols=lambda column: column in {"step", "smiles"})
    membership = pd.read_csv(membership_path, usecols=["source", "source_index", "canon_smiles"])
    membership = membership[membership["source"] == source].copy()
    membership["source_index"] = membership["source_index"].astype(int)
    frame = raw.iloc[membership["source_index"].to_numpy()].copy()
    frame["source_index"] = membership["source_index"].to_numpy()
    frame["canon_smiles"] = membership["canon_smiles"].to_numpy()
    frame["source"] = source
    frame["step"] = pd.to_numeric(frame["step"], errors="coerce")
    return frame.dropna(subset=["step", "canon_smiles"]).reset_index(drop=True)


def murcko_smiles(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return "[NO_SCAFFOLD]" if scaffold.GetNumAtoms() == 0 else Chem.MolToSmiles(scaffold, canonical=True)


def pairwise_similarity(fingerprints) -> dict[str, float | int | None]:
    if len(fingerprints) < 2:
        return {"pairs": 0, "mean_similarity": None, "median_similarity": None,
                "q95_similarity": None, "diversity": None,
                "mean_nearest_similarity_within_scope": None,
                "median_nearest_similarity_within_scope": None,
                "fraction_nearest_ge_0_8_within_scope": None}
    similarities: list[float] = []
    nearest = np.zeros(len(fingerprints), dtype=float)
    for index in range(1, len(fingerprints)):
        current = np.asarray(DataStructs.BulkTanimotoSimilarity(
            fingerprints[index], fingerprints[:index]), dtype=float)
        similarities.extend(current.tolist()); nearest[index] = current.max()
        nearest[:index] = np.maximum(nearest[:index], current)
    values = np.asarray(similarities)
    # Fingerprint diversity: D_FP = 1 - mean[T(A_i, A_j)] for all pairs i < j.
    return {"pairs": len(values), "mean_similarity": float(values.mean()),
            "median_similarity": float(np.median(values)),
            "q95_similarity": float(np.quantile(values, .95)),
            "diversity": float(1 - values.mean()),
            "mean_nearest_similarity_within_scope": float(nearest.mean()),
            "median_nearest_similarity_within_scope": float(np.median(nearest)),
            "fraction_nearest_ge_0_8_within_scope": float((nearest >= .8).mean())}


def scaffold_metrics(scaffolds) -> dict[str, float | int | None]:
    counts = Counter(scaffolds); total = sum(counts.values()); unique = len(counts)
    # p_s = count(s)/N; H = -Σ p_s ln(p_s); H_norm = H/ln(unique scaffolds).
    probabilities = np.asarray(list(counts.values()), dtype=float) / total
    entropy = float(-(probabilities * np.log(probabilities)).sum())
    ordered = sorted(counts.values(), reverse=True)
    singletons = sum(value == 1 for value in counts.values())
    return {"molecules": total, "unique_scaffolds": unique,
            "unique_scaffold_fraction": unique / total if total else None,
            "singleton_scaffolds": singletons,
            "singleton_scaffold_fraction": singletons / unique if unique else None,
            "top1_scaffold_share": ordered[0] / total if total else None,
            "top5_scaffold_share": sum(ordered[:5]) / total if total else None,
            "shannon_entropy": entropy,
            "normalized_entropy": entropy / math.log(unique) if unique > 1 else 0.0}


def _sample(frame: pd.DataFrame, size: int, rng) -> pd.DataFrame:
    if len(frame) <= size:
        return frame
    return frame.loc[rng.choice(frame.index.to_numpy(), size=size, replace=False)]


def analyze_source(frame: pd.DataFrame, source: str, *, last_n: int = 50,
                   fp_radius: int = 2, fp_size: int = 16384,
                   diversity_sample_size: int = 3000,
                   novelty_reference_size: int = 5000,
                   random_state: int = 42) -> dict[str, pd.DataFrame]:
    rng = np.random.default_rng(random_state)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=fp_radius, fpSize=fp_size)
    work = frame.copy(); work["mol"] = work["canon_smiles"].map(Chem.MolFromSmiles)
    work = work[work["mol"].notna()].copy(); work["fp"] = work["mol"].map(generator.GetFingerprint)
    work["murcko_scaffold"] = work["canon_smiles"].map(murcko_smiles)
    work = work[work["murcko_scaffold"].notna()].copy()
    max_step = int(work["step"].max()); first = max(int(work["step"].min()), max_step - last_n + 1)
    recent = work[work["step"].between(first, max_step)].copy()
    fp_rows = []; scaffold_rows = []; seen = set(); previous = pd.DataFrame()
    for step in sorted(work["step"].unique()):
        current = work[work["step"] == step]; similarity = pairwise_similarity(current["fp"].tolist())
        novelty = {"mean_nearest_similarity_to_prior": None, "median_nearest_similarity_to_prior": None,
                   "fraction_nearest_below_0_4": None, "fraction_nearest_below_0_6": None}
        if not previous.empty:
            refs = _sample(previous, novelty_reference_size, rng)["fp"].tolist()
            # Novelty proxy for molecule i: max_j T(fp_i, prior_fp_j).
            # Lower nearest-prior similarity means greater fingerprint novelty.
            nearest = np.asarray([max(DataStructs.BulkTanimotoSimilarity(fp, refs)) for fp in current["fp"]])
            novelty = {"mean_nearest_similarity_to_prior": float(nearest.mean()),
                       "median_nearest_similarity_to_prior": float(np.median(nearest)),
                       "fraction_nearest_below_0_4": float((nearest < .4).mean()),
                       "fraction_nearest_below_0_6": float((nearest < .6).mean())}
        fp_rows.append({"source": source, "step": int(step), "molecules": len(current), **similarity, **novelty})
        scaffolds = current["murcko_scaffold"].tolist(); unique_current = set(scaffolds)
        scaffold_rows.append({"source": source, "step": int(step), **scaffold_metrics(scaffolds),
                              "new_scaffolds": len(unique_current - seen),
                              "new_scaffold_fraction": len(unique_current - seen) / len(unique_current) if unique_current else None})
        seen.update(unique_current); previous = pd.concat([previous, current], ignore_index=True)
    aggregate = []
    for scope, subset in (("full", work), ("last_window", recent)):
        sampled = _sample(subset, diversity_sample_size, rng)
        aggregate.append({"source": source, "scope": scope, "first_step": int(subset["step"].min()),
                          "last_step": int(subset["step"].max()), "population_molecules": len(subset),
                          "sampled_molecules": len(sampled), **pairwise_similarity(sampled["fp"].tolist()),
                          **scaffold_metrics(subset["murcko_scaffold"].tolist())})
    counts = Counter(recent["murcko_scaffold"]); top = [
        {"source": source, "rank": rank, "murcko_scaffold": scaffold,
         "count": count, "share": count / len(recent)}
        for rank, (scaffold, count) in enumerate(counts.most_common(25), 1)]
    return {"fp_by_step": pd.DataFrame(fp_rows), "scaffold_by_step": pd.DataFrame(scaffold_rows),
            "aggregate": pd.DataFrame(aggregate), "top_scaffolds_last_window": pd.DataFrame(top)}
