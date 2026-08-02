"""Run pose QC and/or residue-interaction analysis for one docking pose file."""

import argparse
from pathlib import Path
import time

import pandas as pd

try:
    from .input_models import SingleComplexInput, load_pose_metadata, materialize_sdf
    from .fingerprint import run_prolif
    from .pose_quality import (
        apply_pose_quality_thresholds,
        run_posecheck,
        run_posecheck_parallel,
    )
    from .receptor_mapping import annotate_interactions, load_residue_mapping
except ImportError:
    from input_models import SingleComplexInput, load_pose_metadata, materialize_sdf
    from fingerprint import run_prolif
    from pose_quality import (
        apply_pose_quality_thresholds,
        run_posecheck,
        run_posecheck_parallel,
    )
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
    posecheck_workers=1,
    posecheck_chunk_size=500,
    prolif_workers=None,
    resume=False,
):
    """Analyze one receptor with all poses in one SDF/SDFGZ file."""
    inputs = SingleComplexInput.from_paths(
        protein_file, pose_file, receptor_state
    )
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    print(
        f"[{inputs.receptor_state}] Reading pose file: {inputs.pose_file}",
        flush=True,
    )
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
    print(
        f"[{inputs.receptor_state}] Valid pose records: {len(metadata)}",
        flush=True,
    )

    statuses = []
    if analysis in {"posecheck", "both"}:
        component_started = time.monotonic()
        print(f"[{inputs.receptor_state}] PoseCheck started", flush=True)
        if int(posecheck_workers) > 1:
            quality = run_posecheck_parallel(
                inputs.protein_file,
                sdf_file,
                metadata,
                output_dir / "posecheck_chunks",
                workers=posecheck_workers,
                chunk_size=posecheck_chunk_size,
                resume=resume,
            )
        else:
            quality = run_posecheck(inputs.protein_file, sdf_file, metadata)
        quality = apply_pose_quality_thresholds(
            quality,
            max_clashes=max_clashes,
            max_clashes_per_heavy_atom=max_clashes_per_heavy_atom,
            max_strain_energy=max_strain_energy,
        )
        quality.to_csv(output_dir / "pose_quality.csv", index=False)
        statuses.append(
            (
                "posecheck",
                "completed",
                len(quality),
                time.monotonic() - component_started,
            )
        )
        print(f"[{inputs.receptor_state}] PoseCheck completed", flush=True)

    if analysis in {"prolif", "both"}:
        component_started = time.monotonic()
        print(f"[{inputs.receptor_state}] ProLIF started", flush=True)
        interactions, _ = run_prolif(
            inputs.protein_file,
            sdf_file,
            include_secondary=include_secondary_interactions,
            n_jobs=prolif_workers,
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
        statuses.append(
            (
                "prolif",
                "completed",
                len(interactions),
                time.monotonic() - component_started,
            )
        )
        print(
            f"[{inputs.receptor_state}] ProLIF completed: "
            f"{len(interactions)} interaction occurrence(s)",
            flush=True,
        )

    summary = pd.DataFrame(
        statuses,
        columns=[
            "analysis_component",
            "status",
            "output_rows",
            "elapsed_seconds",
        ],
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
    parser.add_argument("--posecheck-workers", type=int, default=1)
    parser.add_argument("--posecheck-chunk-size", type=int, default=500)
    parser.add_argument("--prolif-workers", type=int)
    parser.add_argument("--resume", action="store_true")
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
        posecheck_workers=args.posecheck_workers,
        posecheck_chunk_size=args.posecheck_chunk_size,
        prolif_workers=args.prolif_workers,
        resume=args.resume,
    )



# example...
# nice -n 15 taskset -c 0 \
# python -m analysis.MP_interaction.single_complex \
#   --protein /home/andy/proj/DProject_Workspace/AHC_Related/ref_structures/Structures/PPS_Protein_no_conect.pdb \
#   --poses /home/andy/proj/DProject_Workspace/AHC_Related/ref_structures/Structures/PPS_ligand.sdfgz \
#   --receptor-state PPS \
#   --output-dir /home/andy/proj/Results/reference_interaction/PPS \
#   --analysis both \
#   --posecheck-workers 1 \
#   --prolif-workers 1



# nice -n 15 taskset -c 1 \
# python -m analysis.MP_interaction.single_complex \
#   --protein /home/andy/proj/DProject_Workspace/AHC_Related/ref_structures/Structures/PR_Protein_no_conect.pdb \
#   --poses /home/andy/proj/DProject_Workspace/AHC_Related/ref_structures/Structures/PR_ligand.sdfgz \
#   --receptor-state PR \
#   --output-dir /home/andy/proj/Results/reference_interaction/PR \
#   --analysis both \
#   --posecheck-workers 1 \
#   --prolif-workers 1