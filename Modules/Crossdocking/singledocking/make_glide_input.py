"""Create a minimal Glide docking input file."""

import argparse
from pathlib import Path


def make_glide_input(prepared_ligand, grid_file, output_in, precision="SP"):
    prepared_ligand = Path(prepared_ligand).resolve()
    grid_file = Path(grid_file).resolve()
    output_in = Path(output_in).resolve()

    if not prepared_ligand.is_file():
        raise FileNotFoundError(f"Prepared ligand not found: {prepared_ligand}")
    if not grid_file.is_file():
        raise FileNotFoundError(f"Receptor grid not found: {grid_file}")

    precision = precision.upper()
    if precision not in {"HTVS", "SP", "XP"}:
        raise ValueError(" prec typeer... must be one of HTVS, SP, XP")

    lines = [
        "FORCEFIELD   OPLS4",
        f"GRIDFILE   {grid_file}",
        f"LIGANDFILE   {prepared_ligand}",
        f"PRECISION   {precision}",
        "POSE_OUTTYPE   ligandlib_sd",
        f"OUTPUTDIR   {output_in.parent}",
    ]

    output_in.parent.mkdir(parents=True, exist_ok=True)
    output_in.write_text("\n".join(lines) + "\n")
    print(f"Glide input: {output_in}")
    return output_in


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a Glide .in file")
    parser.add_argument("--ligand", required=True, help="LigPrep output SDF")
    parser.add_argument("--grid", required=True, help="Glide receptor grid ZIP")
    parser.add_argument("--output", required=True, help="Output .in file")
    parser.add_argument("--precision", default="SP", choices=["HTVS", "SP", "XP"])
    args = parser.parse_args()

    make_glide_input(args.ligand, args.grid, args.output, args.precision)
