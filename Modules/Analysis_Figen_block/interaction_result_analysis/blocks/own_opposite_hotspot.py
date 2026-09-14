"""Own-state versus opposite-state residue prevalence analysis."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


SCRIPT_ID = "PY-125"
ARM_SPECS = {"L_PPS_to_R_PR": {"ligand": "PPS"}, "L_PR_to_R_PPS": {"ligand": "PR"}}


def _residue(frame):
    number = pd.to_numeric(frame["original_residue_number"], errors="coerce")
    return frame["residue_name"].astype(str) + number.astype("Int64").astype(str)


def load_cohorts(metadata_csv: str | Path) -> dict[str, set[str]]:
    meta = pd.read_csv(metadata_csv, usecols=["crossdock_arm", "molecule_id", "pose_read_status"])
    meta = meta[meta.pose_read_status.eq("valid")]
    return {arm: set(meta.loc[meta.crossdock_arm.eq(arm), "molecule_id"].astype(str)) for arm in ARM_SPECS}


def select_best_own_poses(metadata_csv: str | Path, cohorts: dict[str, set[str]]) -> pd.DataFrame:
    cols = ["population_state", "molecule_id", "pose_index", "docking_score", "pose_read_status"]
    meta = pd.read_csv(metadata_csv, usecols=cols); meta = meta[meta.pose_read_status.eq("valid")].copy()
    meta["molecule_id"] = meta.molecule_id.astype(str)
    meta["docking_score"] = pd.to_numeric(meta.docking_score, errors="coerce")
    selected = []
    for arm, spec in ARM_SPECS.items():
        part = meta[meta.population_state.eq(spec["ligand"]) & meta.molecule_id.isin(cohorts[arm])].dropna(subset=["docking_score"])
        selected.append(part.sort_values("docking_score").drop_duplicates("molecule_id").assign(crossdock_arm=arm))
    return pd.concat(selected, ignore_index=True)


def summarize_own(interactions_csv, selected, denominators):
    cols = ["population_state", "pose_index", "molecule_id", "residue_name", "original_residue_number", "interaction_type"]
    interactions = pd.read_csv(interactions_csv, usecols=cols); interactions["molecule_id"] = interactions.molecule_id.astype(str)
    summaries = []
    for arm, spec in ARM_SPECS.items():
        keys = selected.loc[selected.crossdock_arm.eq(arm), ["population_state", "pose_index", "molecule_id"]]
        part = interactions[interactions.population_state.eq(spec["ligand"])].merge(keys,
            on=["population_state", "pose_index", "molecule_id"], how="inner", validate="many_to_one")
        part["residue"] = _residue(part)
        counts = part.dropna(subset=["residue"])[["molecule_id", "residue"]].drop_duplicates().groupby("residue").molecule_id.nunique()
        summary = counts.rename("molecule_count").reset_index(); summary["prevalence"] = summary.molecule_count / denominators[arm]
        summary["crossdock_arm"] = arm; summary["condition"] = "own"; summaries.append(summary)
    return pd.concat(summaries, ignore_index=True)


def summarize_opposite(interactions_csv, denominators):
    cols = ["crossdock_arm", "molecule_id", "residue_name", "original_residue_number", "interaction_type"]
    interactions = pd.read_csv(interactions_csv, usecols=cols); interactions["molecule_id"] = interactions.molecule_id.astype(str)
    interactions["residue"] = _residue(interactions)
    unique = interactions.dropna(subset=["residue"])[["crossdock_arm", "molecule_id", "residue"]].drop_duplicates()
    summary = unique.groupby(["crossdock_arm", "residue"]).molecule_id.nunique().rename("molecule_count").reset_index()
    summary["prevalence"] = summary.apply(lambda row: row.molecule_count / denominators[row.crossdock_arm], axis=1)
    summary["condition"] = "opposite"; return summary

def summarize_own_types(interactions_csv, selected, denominators):
    cols = ["population_state", "pose_index", "molecule_id", "residue_name", "original_residue_number", "interaction_type"]
    interactions = pd.read_csv(interactions_csv, usecols=cols); interactions["molecule_id"] = interactions.molecule_id.astype(str)
    parts = []
    for arm, spec in ARM_SPECS.items():
        keys = selected.loc[selected.crossdock_arm.eq(arm), ["population_state", "pose_index", "molecule_id"]]
        part = interactions[interactions.population_state.eq(spec["ligand"])].merge(keys,
            on=["population_state", "pose_index", "molecule_id"], how="inner", validate="many_to_one")
        part["residue"] = _residue(part)
        part = part.dropna(subset=["residue", "interaction_type"]).drop_duplicates(
            ["population_state", "molecule_id", "residue", "interaction_type"])
        counts = part.groupby(["residue", "population_state", "interaction_type"]).molecule_id.nunique().rename("poses").reset_index()
        counts["prevalence"] = counts.poses / denominators[arm]; parts.append(counts)
    return pd.concat(parts, ignore_index=True)


def build_own_opposite_tables(*, own_metadata, own_interactions, opposite_metadata, opposite_interactions):
    cohorts = load_cohorts(opposite_metadata); selected = select_best_own_poses(own_metadata, cohorts)
    denominators = {arm: len(cohorts[arm].intersection(set(selected.loc[selected.crossdock_arm.eq(arm), "molecule_id"]))) for arm in ARM_SPECS}
    if any(value <= 0 for value in denominators.values()): raise ValueError(f"No matched molecules for an arm: {denominators}")
    summary = pd.concat([summarize_own(own_interactions, selected, denominators),
                         summarize_opposite(opposite_interactions, denominators)], ignore_index=True)
    typed = summarize_own_types(own_interactions, selected, denominators)
    return selected, summary, typed, denominators
