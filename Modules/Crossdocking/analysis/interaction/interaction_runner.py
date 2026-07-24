"""Run the complete pose-quality and interaction-analysis block."""

import argparse

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", required=True, help="Prepared receptor PDB or MOL2")
    parser.add_argument("--poses", required=True, help="Docked ligand SDF or SDFGZ")
    parser.add_argument("--receptor-state", required=True, choices=["PPS", "PR"])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--analysis",
        default="both",
        choices=["posecheck", "prolif", "both"],
    )
    parser.add_argument("--residue-map", help="Optional common-residue mapping CSV")
    parser.add_argument("--include-secondary-interactions", action="store_true")
    parser.add_argument("--max-clashes", type=float)
    parser.add_argument("--max-clashes-per-heavy-atom", type=float)
    parser.add_argument("--max-strain-energy", type=float)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    from .single_complex import run_single_complex

    return run_single_complex(
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


if __name__ == "__main__":
    main()
