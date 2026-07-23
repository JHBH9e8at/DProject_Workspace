"""ProLIF adapter producing an auditable long-format interaction table."""

import json
from pathlib import Path

import pandas as pd


PRIMARY_INTERACTIONS = (
    "HBDonor",
    "HBAcceptor",
    "Anionic",
    "Cationic",
    "PiStacking",
    "CationPi",
    "PiCation",
    "XBDonor",
    "XBAcceptor",
)

SECONDARY_INTERACTIONS = ("Hydrophobic", "VdWContact")


def _require_prolif():
    try:
        import MDAnalysis as mda
        import prolif as plf
    except ImportError as error:
        raise RuntimeError(
            "ProLIF and MDAnalysis are required for interaction fingerprinting. "
            "Install both in a compatible analysis environment."
        ) from error
    return mda, plf


def _json_value(value):
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, tuple):
        value = list(value)
    return json.dumps(
        value,
        sort_keys=True,
        default=lambda item: item.item() if hasattr(item, "item") else str(item),
    )


def _normalize_occurrences(details):
    if isinstance(details, tuple):
        return details
    if isinstance(details, list):
        return tuple(details)
    return (details,)


def _protein_has_explicit_hydrogens(universe):
    try:
        elements = [str(value).upper() for value in universe.atoms.elements]
        if elements:
            return "H" in elements
    except (AttributeError, KeyError, TypeError):
        pass
    names = [str(value).upper() for value in universe.atoms.names]
    return any(name.startswith("H") or (name[:1].isdigit() and "H" in name) for name in names)


def run_prolif(
    protein_file,
    ligand_sdf,
    interactions=PRIMARY_INTERACTIONS,
    include_secondary=False,
    count=True,
):
    """Generate ProLIF fingerprints and return occurrence-level interaction rows."""
    protein_file = Path(protein_file).resolve()
    ligand_sdf = Path(ligand_sdf).resolve()
    mda, plf = _require_prolif()

    selected = list(interactions)
    if include_secondary:
        selected.extend(
            interaction
            for interaction in SECONDARY_INTERACTIONS
            if interaction not in selected
        )

    universe = mda.Universe(str(protein_file))
    if not _protein_has_explicit_hydrogens(universe):
        raise ValueError(
            "The receptor contains no explicit hydrogens. ProLIF interaction "
            "typing requires the prepared hydrogen-containing receptor structure."
        )
    protein = plf.Molecule.from_mda(universe)
    poses = list(plf.sdf_supplier(str(ligand_sdf)))
    if not poses:
        raise ValueError(f"ProLIF found no valid ligand poses: {ligand_sdf}")

    fingerprint = plf.Fingerprint(interactions=selected, count=count)
    fingerprint.run_from_iterable(poses, protein)

    rows = []
    for pose_index, residue_pairs in fingerprint.ifp.items():
        for (ligand_residue, protein_residue), interaction_map in residue_pairs.items():
            for interaction_type, details in interaction_map.items():
                occurrences = _normalize_occurrences(details)
                for occurrence_index, detail in enumerate(occurrences, start=1):
                    detail = detail or {}
                    indices = detail.get("indices", {})
                    geometry = {
                        key: value
                        for key, value in detail.items()
                        if key not in {"indices", "parent_indices"}
                    }
                    rows.append(
                        {
                            "pose_index": int(pose_index),
                            "ligand_residue": str(ligand_residue),
                            "protein_residue": str(protein_residue),
                            "interaction_type": str(interaction_type),
                            "interaction_present": True,
                            "interaction_occurrence": occurrence_index,
                            "ligand_atom_indices": _json_value(
                                indices.get("ligand", [])
                            ),
                            "protein_atom_indices": _json_value(
                                indices.get("protein", [])
                            ),
                            "interaction_geometry_json": _json_value(geometry),
                        }
                    )

    columns = [
        "pose_index",
        "ligand_residue",
        "protein_residue",
        "interaction_type",
        "interaction_present",
        "interaction_occurrence",
        "ligand_atom_indices",
        "protein_atom_indices",
        "interaction_geometry_json",
    ]
    return pd.DataFrame(rows, columns=columns), fingerprint
