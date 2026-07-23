"""Input validation and common pose metadata for interaction analysis."""

from dataclasses import dataclass
import gzip
from pathlib import Path
import shutil

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Lipinski


SUPPORTED_PROTEIN_SUFFIXES = {".pdb", ".mol2"}
SUPPORTED_POSE_SUFFIXES = {".sdf", ".sdfgz"}


@dataclass(frozen=True)
class SingleComplexInput:
    """Validated files and state labels for one protein and one pose collection."""

    protein_file: Path
    pose_file: Path
    receptor_state: str

    @classmethod
    def from_paths(cls, protein_file, pose_file, receptor_state):
        protein = Path(protein_file).resolve()
        poses = Path(pose_file).resolve()
        state = str(receptor_state).strip().upper()

        if not protein.is_file():
            raise FileNotFoundError(f"Protein structure not found: {protein}")
        if protein.suffix.lower() not in SUPPORTED_PROTEIN_SUFFIXES:
            raise ValueError(
                f"Unsupported protein format: {protein}. Expected .pdb or .mol2."
            )
        if not poses.is_file():
            raise FileNotFoundError(f"Docking pose file not found: {poses}")
        if poses.suffix.lower() not in SUPPORTED_POSE_SUFFIXES:
            raise ValueError(
                f"Unsupported pose format: {poses}. Expected .sdf or .sdfgz."
            )
        if state not in {"PPS", "PR"}:
            raise ValueError("receptor_state must be PPS or PR")

        return cls(protein, poses, state)


def materialize_sdf(pose_file, output_dir):
    """Return an uncompressed SDF, expanding SDFGZ into an auditable directory."""
    pose_file = Path(pose_file).resolve()
    if pose_file.suffix.lower() == ".sdf":
        return pose_file
    if pose_file.suffix.lower() != ".sdfgz":
        raise ValueError(f"Expected .sdf or .sdfgz: {pose_file}")

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / pose_file.with_suffix(".sdf").name

    if output_file.is_file() and output_file.stat().st_size > 0:
        return output_file
    if output_file.exists():
        output_file.unlink()

    try:
        with gzip.open(pose_file, "rb") as source:
            with output_file.open("wb") as destination:
                shutil.copyfileobj(source, destination)
    except Exception:
        if output_file.exists():
            output_file.unlink()
        raise
    return output_file


def load_pose_metadata(sdf_file):
    """Read stable identifiers and basic ligand properties from every valid pose."""
    sdf_file = Path(sdf_file).resolve()
    rows = []
    supplier = Chem.ForwardSDMolSupplier(str(sdf_file), removeHs=False)

    for record_index, mol in enumerate(supplier, start=1):
        if mol is None:
            rows.append(
                {
                    "pose_index": record_index - 1,
                    "sdf_record": record_index,
                    "pose_read_status": "invalid",
                }
            )
            continue

        properties = mol.GetPropsAsDict()
        rows.append(
            {
                "pose_index": record_index - 1,
                "sdf_record": record_index,
                "pose_read_status": "valid",
                "pose_title": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
                "docking_score": pd.to_numeric(
                    properties.get("r_i_docking_score"), errors="coerce"
                ),
                "glide_variant": properties.get("s_lp_Variant", ""),
                "heavy_atom_count": mol.GetNumHeavyAtoms(),
                "rotatable_bond_count": Lipinski.NumRotatableBonds(mol),
                "formal_charge": Chem.GetFormalCharge(mol),
            }
        )

    metadata = pd.DataFrame(rows)
    if metadata.empty:
        raise ValueError(f"No SDF records found: {sdf_file}")
    return metadata
