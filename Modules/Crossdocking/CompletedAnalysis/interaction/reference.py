"""Generate stable PR/PPS reference-interaction inputs for hotspot figures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .common.single_complex import run_single_complex


ABBREVIATIONS = {
    "HBDonor": "HBD", "HBAcceptor": "HBA", "PiStacking": "pi",
    "PiCation": "pi-cat", "CationPi": "cat-pi", "Cationic": "cat",
    "Anionic": "ani", "XBDonor": "XB",
}


def combine_reference_interactions(pps_csv: str | Path, pr_csv: str | Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    for state, source in (("PPS", pps_csv), ("PR", pr_csv)):
        frame = pd.read_csv(
            source,
            usecols=["residue_name", "original_residue_number", "interaction_type"],
        ).dropna()
        numbers = pd.to_numeric(frame["original_residue_number"], errors="raise").astype(int)
        frame["residue"] = frame["residue_name"].astype(str) + numbers.astype(str)
        frame["label"] = frame["interaction_type"].map(ABBREVIATIONS).fillna(
            frame["interaction_type"].astype(str)
        )
        for residue, group in frame.groupby("residue", sort=False):
            labels = list(dict.fromkeys(group["label"].astype(str)))
            records.append({"residue": residue, "state": state, "reference": "/".join(labels)})
    if not records:
        return pd.DataFrame(columns=["residue", "PPS", "PR"])
    table = pd.DataFrame(records).pivot(index="residue", columns="state", values="reference")
    table = table.reindex(columns=["PPS", "PR"]).fillna("-").reset_index()
    table.columns.name = None
    return table


def run_reference_interactions(
    *, pps_receptor, pr_receptor, pps_ligand, pr_ligand, results_root,
    run_id, include_secondary_interactions=False, prolif_workers=1, resume=False,
) -> Path:
    output_root = Path(results_root) / "analysis" / "reference_interactions" / str(run_id)
    if output_root.exists() and not resume:
        raise FileExistsError(
            f"Reference-interaction run already exists: {output_root}. "
            "Use resume=on to reuse it or choose a new run_id."
        )
    state_inputs = {
        "PPS": (pps_receptor, pps_ligand),
        "PR": (pr_receptor, pr_ligand),
    }
    interaction_files: dict[str, Path] = {}
    for state, (protein, ligand) in state_inputs.items():
        state_dir = output_root / state
        interaction_file = state_dir / "interactions_long.csv"
        if not (resume and interaction_file.is_file()):
            run_single_complex(
                protein_file=protein,
                pose_file=ligand,
                receptor_state=state,
                output_dir=state_dir,
                analysis="prolif",
                include_secondary_interactions=include_secondary_interactions,
                prolif_workers=prolif_workers,
                resume=resume,
            )
        interaction_files[state] = interaction_file
    combined = combine_reference_interactions(
        interaction_files["PPS"], interaction_files["PR"]
    )
    output_file = output_root / "reference_interactions.csv"
    output_root.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, index=False)
    return output_file
