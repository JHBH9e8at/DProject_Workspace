"""Dock one prepared SDF or SDFGZ ligand file against one Glide grid.
This is a utility rnr for small tests. It does not use Dask or LigPrep!!.
"""

import argparse
from pathlib import Path

from docking.glide_input import make_glide_input
from docking.glide_execution import run_glide
from docking.scores import extract_scores
from prep.unpack_sdfgz import unpack_sdfgz


def prepare_input_sdf(ligand_file, output_dir):
    """Return an SDF path, unpacking an SDFGZ input when necessary."""
    ligand_file = Path(ligand_file).resolve()

    if not ligand_file.is_file():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")

    suffix = ligand_file.suffix.lower()
    if suffix == ".sdf":
        print(f"Using SDF directly: {ligand_file}")
        return ligand_file

    if suffix == ".sdfgz":
        print(f"Unpacking SDFGZ: {ligand_file}")
        return unpack_sdfgz(
            input_file=ligand_file,
            output_dir=output_dir / "input",
            log_dir=output_dir / "logs",
        )

    raise ValueError(f"Unsupported lig form err: {ligand_file}. Expected .sdf or .sdfgz.")


def single_dock(ligand_file, grid_file, output_dir, precision="SP"):
    """Run a single Glide docking job and extract its scores."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # run_glide identifies the generated pose by this pattern. Refuse to mix a
    # new run with an old result so that it cannot select a stale pose file.
    existing_poses = sorted(output_dir.glob("*_lib.sdfgz"))
    if existing_poses:
        raise FileExistsError(
            f"Output directory already contains a Glide pose file: "
            f"{existing_poses[0]}. Use a new --output-dir."
        )

    print("[1/4] Preparing ligand input")
    input_sdf = prepare_input_sdf(ligand_file, output_dir)

    print("[2/4] Creating Glide input")
    glide_input = make_glide_input(
        prepared_ligand=input_sdf,
        grid_file=grid_file,
        output_in=output_dir / "single_docking.in",
        precision=precision,
    )

    print("[3/4] Running Glide")
    pose_file = run_glide(glide_input=glide_input, output_dir=output_dir)

    print("[4/4] Extracting docking scores")
    score_file = output_dir / "docking_scores.csv"
    scores = extract_scores(pose_file=pose_file, output_csv=score_file)

    print("Single docking completed")
    print(f"Pose file: {pose_file}")
    print(f"Score file: {score_file}")
    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dock one prepared SDF/SDFGZ ligand file to one Glide grid"
    )
    parser.add_argument("--ligand", required=True, help="Input .sdf or .sdfgz")
    parser.add_argument("--grid", required=True, help="Glide receptor grid ZIP")
    parser.add_argument("--output-dir", required=True, help="New result directory")
    parser.add_argument(
        "--precision",
        default="SP",
        choices=["HTVS", "SP", "XP"],
        help="Glide docking precision (default: SP)",
    )
    args = parser.parse_args()

    single_dock(args.ligand, args.grid, args.output_dir, args.precision)
