"""PoseCheck adapter and transparent pose-quality thresholding."""

from pathlib import Path

import numpy as np
import pandas as pd


def _require_posecheck():
    try:
        from posecheck import PoseCheck
    except ImportError as error:
        raise RuntimeError(
            "PoseCheck is required for pose-quality analysis but is not installed. "
            "Install it in a dedicated compatible environment before running this step."
        ) from error
    return PoseCheck


def run_posecheck(protein_pdb, ligand_sdf, pose_metadata=None):
    """Calculate raw clash and strain metrics without constructing a composite PQI."""
    protein_pdb = Path(protein_pdb).resolve()
    ligand_sdf = Path(ligand_sdf).resolve()
    if protein_pdb.suffix.lower() != ".pdb":
        raise ValueError("PoseCheck currently requires the protein input as .pdb")

    PoseCheck = _require_posecheck()
    checker = PoseCheck()
    checker.load_protein_from_pdb(str(protein_pdb))
    checker.load_ligands_from_sdf(str(ligand_sdf))

    clashes = list(checker.calculate_clashes())
    strain = list(checker.calculate_strain_energy())
    if len(clashes) != len(strain):
        raise RuntimeError(
            "PoseCheck returned different numbers of clash and strain results"
        )

    quality = pd.DataFrame(
        {
            "pose_index": range(len(clashes)),
            "clash_count": pd.to_numeric(clashes, errors="coerce"),
            "strain_energy": pd.to_numeric(strain, errors="coerce"),
        }
    )
    quality["posecheck_status"] = np.where(
        np.isfinite(quality["clash_count"])
        & np.isfinite(quality["strain_energy"]),
        "completed",
        "invalid_metric",
    )

    if pose_metadata is not None:
        quality = pose_metadata.merge(
            quality, on="pose_index", how="left", validate="one_to_one"
        )
        valid_count = int((pose_metadata["pose_read_status"] == "valid").sum())
        if valid_count != len(clashes):
            raise RuntimeError(
                "PoseCheck result count does not match the number of valid SDF poses: "
                f"{len(clashes)} versus {valid_count}"
            )

    if "heavy_atom_count" in quality.columns:
        heavy_atoms = pd.to_numeric(quality["heavy_atom_count"], errors="coerce")
        quality["clashes_per_heavy_atom"] = quality["clash_count"] / heavy_atoms
        quality.loc[heavy_atoms <= 0, "clashes_per_heavy_atom"] = np.nan
    if "rotatable_bond_count" in quality.columns:
        rotors = pd.to_numeric(quality["rotatable_bond_count"], errors="coerce")
        quality["strain_per_rotatable_bond"] = quality["strain_energy"] / rotors
        quality.loc[rotors <= 0, "strain_per_rotatable_bond"] = np.nan
    return quality


def apply_pose_quality_thresholds(
    quality,
    max_clashes=None,
    max_clashes_per_heavy_atom=None,
    max_strain_energy=None,
):
    """Apply explicit QC thresholds while preserving every raw metric."""
    assessed = quality.copy()
    assessed["pose_quality_pass"] = assessed["posecheck_status"].eq("completed")
    failure_reasons = pd.Series("", index=assessed.index, dtype="object")

    criteria = [
        ("clash_count", max_clashes, "clash_count"),
        (
            "clashes_per_heavy_atom",
            max_clashes_per_heavy_atom,
            "clashes_per_heavy_atom",
        ),
        ("strain_energy", max_strain_energy, "strain_energy"),
    ]
    for column, threshold, label in criteria:
        if threshold is None:
            continue
        if threshold < 0:
            raise ValueError(f"{label} threshold must be nonnegative")
        if column not in assessed.columns:
            raise ValueError(f"Quality table does not contain {column}")
        failed = pd.to_numeric(assessed[column], errors="coerce") > threshold
        assessed.loc[failed, "pose_quality_pass"] = False
        failure_reasons.loc[failed] = failure_reasons.loc[failed].map(
            lambda current: f"{current};{label}" if current else label
        )

    invalid = ~assessed["posecheck_status"].eq("completed")
    failure_reasons.loc[invalid] = "invalid_posecheck_metric"
    assessed["pose_quality_failure_reason"] = failure_reasons
    return assessed
