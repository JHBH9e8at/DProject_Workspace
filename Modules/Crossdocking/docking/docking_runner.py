"""Run the complete Glide docking block for one prepared ligand file."""

import argparse
from pathlib import Path

def run_docking(
    ligand_file,
    grid_file,
    output_dir,
    precision="SP",
    resume=False,
):
    """Create Glide input, run docking, and write all-pose and best-score CSVs."""
    from .glide_execution import run_glide
    from .glide_input import make_glide_input
    from .scores import extract_scores

    ligand_file = Path(ligand_file).resolve()
    grid_file = Path(grid_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    glide_input = output_dir / "docking.in"
    score_file = output_dir / "docking_scores.csv"
    best_score_file = output_dir / "best_docking_score.csv"

    pose_files = sorted(output_dir.glob("*_lib.sdfgz"))
    if resume and pose_files:
        pose_file = pose_files[0]
        print(f"Reusing Glide pose file: {pose_file}")
    else:
        make_glide_input(
            prepared_ligand=ligand_file,
            grid_file=grid_file,
            output_in=glide_input,
            precision=precision,
        )
        pose_file = run_glide(glide_input=glide_input, output_dir=output_dir)

    if resume and score_file.is_file():
        import pandas as pd

        scores = pd.read_csv(score_file)
        print(f"Reusing docking scores: {score_file}")
    else:
        scores = extract_scores(pose_file=pose_file, output_csv=score_file)

    best = scores.nsmallest(1, "r_i_docking_score")
    best.to_csv(best_score_file, index=False)

    print("Docking block completed")
    print(f"Pose file: {pose_file}")
    print(f"Scores: {score_file}")
    print(f"Best score: {best_score_file}")
    return {
        "glide_input": glide_input,
        "pose_file": pose_file,
        "score_file": score_file,
        "best_score_file": best_score_file,
    }


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligand", required=True, help="Prepared ligand SDF")
    parser.add_argument("--grid", required=True, help="Glide receptor grid ZIP")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--precision", default="SP", choices=["HTVS", "SP", "XP"])
    parser.add_argument("--resume", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_docking(
        ligand_file=args.ligand,
        grid_file=args.grid,
        output_dir=args.output_dir,
        precision=args.precision,
        resume=args.resume,
    )


if __name__ == "__main__":
    main()
