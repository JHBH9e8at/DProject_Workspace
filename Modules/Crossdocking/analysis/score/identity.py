"""Canonical molecular identity and within-population score aggregation."""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem


def canonicalize_isomeric_smiles(smiles):
    """Return canonical isomeric SMILES without neutralizing or tautomerizing."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None, None
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return canonical, Chem.GetFormalCharge(mol)


def _valid_boolean(series):
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1"})


def prepare_population_instances(
    input_csv,
    population,
    score_col,
    smiles_col="smiles",
    valid_col="valid",
):
    """Load one raw population and return valid canonicalized and excluded rows."""
    input_csv = Path(input_csv).resolve()
    if not input_csv.is_file():
        raise FileNotFoundError(f"Population scores CSV not found: {input_csv}")

    frame = pd.read_csv(input_csv)
    required = [smiles_col, score_col]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{population} scores CSV missing columns: {missing}")

    working = frame.copy()
    if "Unnamed: 0" in working.columns:
        working["source_row"] = working["Unnamed: 0"]
    else:
        working["source_row"] = working.index

    working["_score"] = pd.to_numeric(working[score_col], errors="coerce")
    reasons = pd.Series("", index=working.index, dtype="object")

    missing_smiles = working[smiles_col].isna() | working[smiles_col].astype(str).str.strip().eq("")
    reasons.loc[missing_smiles] = "missing_smiles"

    if valid_col in working.columns:
        invalid_flag = ~_valid_boolean(working[valid_col])
        reasons.loc[(reasons == "") & invalid_flag] = "invalid_flag"

    invalid_score = ~np.isfinite(working["_score"])
    reasons.loc[(reasons == "") & invalid_score] = "invalid_score"

    nonnegative_score = working["_score"] >= 0
    reasons.loc[(reasons == "") & nonnegative_score] = "nonnegative_score"

    canonical_values = {}
    charge_values = {}
    eligible = working.index[reasons == ""]
    for index in eligible:
        canonical, formal_charge = canonicalize_isomeric_smiles(working.at[index, smiles_col])
        if canonical is None:
            reasons.at[index] = "smiles_parse_failed"
        else:
            canonical_values[index] = canonical
            charge_values[index] = formal_charge

    valid = working.loc[reasons == ""].copy()
    valid["population"] = population
    valid["canonical_isomeric_smiles"] = valid.index.map(canonical_values)
    valid["canonical_formal_charge"] = valid.index.map(charge_values)
    valid["docking_score"] = valid["_score"]
    valid["input_smiles"] = valid[smiles_col]

    audit_columns = [
        "population",
        "source_row",
        "input_smiles",
        "canonical_isomeric_smiles",
        "canonical_formal_charge",
        "docking_score",
    ]
    for optional in ("step", "batch_idx", "task", "model", "PPS_best_variant", "PR_best_variant"):
        if optional in valid.columns and optional not in audit_columns:
            audit_columns.append(optional)
    valid_audit = valid[audit_columns].copy()

    excluded = working.loc[reasons != ""].copy()
    excluded.insert(0, "population", population)
    excluded.insert(1, "exclusion_reason", reasons.loc[excluded.index])
    excluded_columns = ["population", "exclusion_reason", "source_row", smiles_col, score_col]
    excluded = excluded[[column for column in excluded_columns if column in excluded.columns]]

    return valid, valid_audit, excluded, len(frame)


def aggregate_population(valid_instances, population, score_col):
    """Collapse repeated generations to one row per canonical molecule."""
    rows = []
    group_col = "canonical_isomeric_smiles"

    for canonical_smiles, group in valid_instances.groupby(group_col, sort=False):
        best_index = group["_score"].idxmin()
        best_row = group.loc[best_index]
        scores = group["_score"]

        row = {
            "canonical_isomeric_smiles": canonical_smiles,
            "canonical_formal_charge": int(best_row["canonical_formal_charge"]),
            f"{population}_occurrence_count": len(group),
            f"{population}_raw_smiles_count": group["input_smiles"].nunique(),
            f"{population}_best_score": float(scores.min()),
            f"{population}_median_score": float(scores.median()),
            f"{population}_mean_score": float(scores.mean()),
            f"{population}_score_std": float(scores.std(ddof=1)) if len(scores) > 1 else 0.0,
            f"{population}_worst_score": float(scores.max()),
            f"{population}_score_range": float(scores.max() - scores.min()),
            f"{population}_best_source_row": best_row["source_row"],
            f"{population}_best_raw_smiles": best_row["input_smiles"],
        }

        if "step" in group.columns:
            row[f"{population}_first_step"] = group["step"].min()
            row[f"{population}_last_step"] = group["step"].max()
            row[f"{population}_best_step"] = best_row["step"]
        if "batch_idx" in group.columns:
            row[f"{population}_best_batch_idx"] = best_row["batch_idx"]

        variant_col = f"{population}_best_variant"
        if variant_col in group.columns:
            row[f"{population}_best_variant_id"] = best_row[variant_col]

        rows.append(row)

    aggregated = pd.DataFrame(rows)
    if aggregated.empty:
        raise ValueError(f"No technically valid {population} molecules remain")
    return aggregated.sort_values(f"{population}_best_score").reset_index(drop=True)
