"""Stable residue-hotspot tables extracted from mixed legacy Figure scripts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


SCRIPT_ID = "PY-122"
ABBREVIATIONS = {"HBDonor": "HBD", "HBAcceptor": "HBA", "PiStacking": "pi",
                 "PiCation": "pi-cat", "CationPi": "cat-pi", "Cationic": "cat",
                 "Anionic": "ani", "XBDonor": "XB"}


def summarize_interactions(interactions_csv: str | Path, denominators: dict[str, int],
                           *, chunksize: int = 250_000):
    required = ["population_state", "molecule_id", "residue_name",
                "original_residue_number", "interaction_type"]
    residue_parts = []; typed_parts = []
    for chunk in pd.read_csv(interactions_csv, usecols=required, chunksize=chunksize):
        chunk = chunk.dropna(subset=required[:-1]).copy()
        chunk["residue_number"] = chunk["original_residue_number"].astype(int)
        chunk["residue"] = chunk["residue_name"].astype(str) + chunk["residue_number"].astype(str)
        residue_parts.append(chunk[["population_state", "molecule_id", "residue"]].drop_duplicates())
        typed_parts.append(chunk[["population_state", "molecule_id", "residue", "interaction_type"]].drop_duplicates())
    if not residue_parts: raise ValueError("Interaction table contains no rows")
    residue_unique = pd.concat(residue_parts, ignore_index=True).drop_duplicates()
    typed_unique = pd.concat(typed_parts, ignore_index=True).drop_duplicates()
    unknown = sorted(set(residue_unique.population_state) - set(denominators))
    if unknown: raise ValueError(f"Missing pose denominators for states: {unknown}")
    if any(value <= 0 for value in denominators.values()): raise ValueError("Pose denominators must be positive")
    residue = residue_unique.groupby(["residue", "population_state"])["molecule_id"].nunique().rename("poses").reset_index()
    residue["prevalence"] = residue.apply(lambda row: row.poses / denominators[row.population_state], axis=1)
    typed = typed_unique.groupby(["residue", "population_state", "interaction_type"])["molecule_id"].nunique().rename("poses").reset_index()
    typed["prevalence"] = typed.apply(lambda row: row.poses / denominators[row.population_state], axis=1)
    dominant = typed.loc[typed.groupby(["residue", "population_state"])["prevalence"].idxmax()].copy()
    dominant["dominant_type"] = dominant.interaction_type.map(ABBREVIATIONS).fillna(dominant.interaction_type)
    return residue, typed, dominant[["residue", "population_state", "dominant_type"]]


def read_pocket_residues(path: str | Path) -> set[str]:
    residues = set()
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name, number = line[17:20].strip(), line[22:26].strip()
            if name and number: residues.add(f"{name}{number}")
    return residues

