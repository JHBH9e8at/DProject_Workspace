"""Run pose QC and/or residue-interaction analysis for one docking pose file."""

import argparse
from pathlib import Path

import pandas as pd

try:
    from .input_models import SingleComplexInput, load_pose_metadata, materialize_sdf
    from .fingerprint import run_prolif
    from .pose_quality import apply_pose_quality_thresholds, run_posecheck
    from .receptor_mapping import annotate_interactions, load_residue_mapping
except ImportError:
    from input_models import SingleComplexInput, load_pose_metadata, materialize_sdf
    from fingerprint import run_prolif
    from pose_quality import apply_pose_quality_thresholds, run_posecheck
    from receptor_mapping import annotate_interactions, load_residue_mapping


def run_single_complex(
    protein_file,
    pose_file,
    receptor_state,
    output_dir,
    analysis="both",
    residue_map=None,
    include_secondary_interactions=False,
    max_clashes=None,
    max_clashes_per_heavy_atom=None,
    max_strain_energy=None,
):
    """Analyze one receptor with all poses in one SDF/SDFGZ file."""
    inputs = SingleComplexInput.from_paths(
        protein_file, pose_file, receptor_state
    )
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    sdf_file = materialize_sdf(inputs.pose_file, output_dir / "intermediate")
    metadata = load_pose_metadata(sdf_file)
    invalid_records = metadata["pose_read_status"].ne("valid")
    if invalid_records.any():
        records = metadata.loc[invalid_records, "sdf_record"].tolist()
        raise ValueError(
            f"Invalid SDF record(s) prevent stable pose indexing: {records}"
        )
    metadata.insert(0, "receptor_state", inputs.receptor_state)
    metadata.to_csv(output_dir / "pose_metadata.csv", index=False)

    statuses = []
    if analysis in {"posecheck", "both"}:
        quality = run_posecheck(inputs.protein_file, sdf_file, metadata)
        quality = apply_pose_quality_thresholds(
            quality,
            max_clashes=max_clashes,
            max_clashes_per_heavy_atom=max_clashes_per_heavy_atom,
            max_strain_energy=max_strain_energy,
        )
        quality.to_csv(output_dir / "pose_quality.csv", index=False)
        statuses.append(("posecheck", "completed", len(quality)))

    if analysis in {"prolif", "both"}:
        interactions, _ = run_prolif(
            inputs.protein_file,
            sdf_file,
            include_secondary=include_secondary_interactions,
        )
        mapping = load_residue_mapping(residue_map) if residue_map else None
        interactions = annotate_interactions(
            interactions, inputs.receptor_state, mapping
        )
        interactions = metadata.merge(
            interactions,
            on=["pose_index", "receptor_state"],
            how="inner",
            validate="one_to_many",
        )
        interactions.to_csv(output_dir / "interactions_long.csv", index=False)
        statuses.append(("prolif", "completed", len(interactions)))

    summary = pd.DataFrame(
        statuses, columns=["analysis_component", "status", "output_rows"]
    )
    summary.insert(0, "receptor_state", inputs.receptor_state)
    summary["protein_file"] = str(inputs.protein_file)
    summary["pose_file"] = str(inputs.pose_file)
    summary["materialized_sdf"] = str(sdf_file)
    summary.to_csv(output_dir / "run_summary.csv", index=False)

    print("Single-complex interaction analysis completed")
    print(f"Analysis: {analysis}")
    print(f"Pose records: {len(metadata)}")
    print(f"Output directory: {output_dir}")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze pose quality and residue interactions for one docked SDF"
    )
    parser.add_argument("--protein", required=True, help="Prepared receptor .pdb or .mol2")
    parser.add_argument("--poses", required=True, help="Docked ligand .sdf or .sdfgz")
    parser.add_argument("--receptor-state", required=True, choices=["PPS", "PR"])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--analysis",
        default="both",
        choices=["posecheck", "prolif", "both"],
    )
    parser.add_argument("--residue-map", help="Optional PPS/PR common residue CSV")
    parser.add_argument(
        "--include-secondary-interactions",
        action="store_true",
        help="Also include Hydrophobic and VdWContact fingerprints",
    )
    parser.add_argument("--max-clashes", type=float)
    parser.add_argument("--max-clashes-per-heavy-atom", type=float)
    parser.add_argument("--max-strain-energy", type=float)
    args = parser.parse_args()

    run_single_complex(
        protein_file=args.protein,
        pose_file=args.poses,
        receptor_state=args.receptor_state,
        output_dir=args.output_dir,
        analysis=args.analysis,
        residue_map=args.residue_map,
        include_secondary_interactions=args.include_secondary_interactions,
        max_clashes=args.max_clashes,
        max_clashes_per_heavy_atom=args.max_clashes_per_heavy_atom,
        max_strain_energy=args.max_strain_energy,
    )
