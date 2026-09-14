from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import generate_crossdock_own_opposite_hotspots as crossdock
import generate_crossdock_own_opposite_hotspots_fullstyle as fullstyle


OUTPUT_STEM = "top5_state_selective_own_opposite_residue_hotspots"


def clean_residue(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.[A-Za-z0-9]+$", "", regex=True)


def load_top5(input_dir: Path):
    summary = pd.read_csv(input_dir / "top5_candidate_summary.csv")
    transitions = pd.read_csv(input_dir / "top5_interaction_transitions.csv")
    summary["own_state"] = summary["own_state"].astype(str)
    transitions["residue"] = clean_residue(transitions["residue"])
    transitions["molecule_id"] = transitions["molecule_id"].astype(str)

    denominators = summary.groupby("own_state")["molecule_id"].nunique().astype(int).to_dict()
    prevalence_frames = []
    for condition, flag in (("own", "own_present"), ("opposite", "opposite_present")):
        present = transitions.loc[transitions[flag].astype(str).str.lower().eq("true")].copy()
        unique = present[["own_state", "molecule_id", "residue"]].drop_duplicates()
        counts = unique.groupby(["residue", "own_state"])["molecule_id"].nunique().rename("poses").reset_index()
        counts["prevalence"] = counts.apply(
            lambda row: row["poses"] / denominators[row["own_state"]], axis=1
        )
        counts = counts.rename(columns={"own_state": "population_state"})
        counts["condition"] = condition
        prevalence_frames.append(counts)

    own = prevalence_frames[0]
    opposite = prevalence_frames[1]
    own_typed = transitions.loc[
        transitions["own_present"].astype(str).str.lower().eq("true")
    ].drop_duplicates(["own_state", "molecule_id", "residue", "interaction_type"])
    typed = (
        own_typed.groupby(["residue", "own_state", "interaction_type"])["molecule_id"]
        .nunique().rename("poses").reset_index()
        .rename(columns={"own_state": "population_state"})
    )
    typed["prevalence"] = typed.apply(
        lambda row: row["poses"] / denominators[row["population_state"]], axis=1
    )
    return own, opposite, typed, denominators


def main(config_path: Path) -> None:
    cfg = crossdock.parse_config(config_path)
    input_dir = crossdock.resolve(
        cfg.get("top5_input_dir", str(crossdock.resolve(cfg["opposite_result_dir"], config_path) / "analysis" / "top5_state_selective")),
        config_path,
    )
    output_dir = crossdock.resolve(cfg["output_dir"], config_path)
    own, opposite, typed, denominators = load_top5(input_dir)
    pocket_root = crossdock.resolve(cfg["fpocket_root"], config_path)
    pockets = {
        f"{state} {label}": crossdock.read_pocket_residues(pocket_root / relative)
        for state, specs in crossdock.POCKET_SPECS.items()
        for label, relative, _ in specs
    }
    reference = fullstyle.base.read_reference_interactions()
    fullstyle.OUTPUT_STEM = OUTPUT_STEM
    fullstyle.draw_figure(
        own, opposite, typed, pockets, reference, denominators, output_dir,
        "Top-5 state-selective candidates: own vs opposite receptor state",
    )
    print(f"Top-5 denominators: PPS={denominators.get('PPS', 0)}, PR={denominators.get('PR', 0)}")
    print(f"Output: {output_dir / (OUTPUT_STEM + '.png')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    main(args.config)
