"""Run LigPrep, Glide docking, and score extraction for one ligand file."""

import argparse
from pathlib import Path

from extract_scores import extract_scores
from make_glide_input import make_glide_input
from prepare_ligand import prepare_ligand
from run_glide import run_glide


def single_dock(ligand_file, grid_file, output_dir, precision="SP"):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    prepared_sdf = output_dir / "prepared_ligand.sdf"
    glide_input = output_dir / "single_docking.in"
    scores_csv = output_dir / "docking_scores.csv"

    print("[1/4] LigPrep")
    prepare_ligand(ligand_file, prepared_sdf)

    print("[2/4] Create Glide input")
    make_glide_input(prepared_sdf, grid_file, glide_input, precision)

    print("[3/4] Glide docking")
    pose_file = run_glide(glide_input, output_dir)

    print("[4/4] Extract scores")
    scores = extract_scores(pose_file, scores_csv)
    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dock one ligand file to one Glide grid")
    parser.add_argument("--ligand", required=True, help="SMI, SDF, or MAE ligand file")
    parser.add_argument("--grid", required=True, help="Glide receptor grid ZIP")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--precision", default="SP", choices=["HTVS", "SP", "XP"])
    args = parser.parse_args()

    single_dock(args.ligand, args.grid, args.output_dir, args.precision)

