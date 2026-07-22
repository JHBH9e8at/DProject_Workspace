"""Prepare one ligand file with Schrodinger LigPrep."""

import argparse
import os
import subprocess
from pathlib import Path


INPUT_FLAGS = {".smi": "-ismi",
               ".smiles": "-ismi",
               ".sdf": "-isd",
               ".sd": "-isd",
               ".mae": "-imae",
               ".maegz": "-imae",
}


def prepare_ligand(ligand_file, output_sdf):
    ligand_file = Path(ligand_file).resolve()
    output_sdf = Path(output_sdf).resolve()

    if not ligand_file.is_file():
        raise FileNotFoundError(f"no lig file err: {ligand_file}")

    input_flag = INPUT_FLAGS.get(ligand_file.suffix.lower())
    if input_flag is None:
        raise ValueError(f"lig Ftype err: {ligand_file.suffix}")

    schrodinger = os.environ.get("SCHRODINGER")
    if not schrodinger:
        raise EnvironmentError("crit SCHRODINGER environ is none")

    ligprep = Path(schrodinger) / "ligprep"
    if not ligprep.is_file():
        raise FileNotFoundError(f"LigPrep module not detec: {ligprep}")

    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(ligprep),
        input_flag,
        str(ligand_file),
        "-osd",
        output_sdf.name,
        "-WAIT",
    ]

    subprocess.run(
        command,
        cwd=output_sdf.parent,
        check=True,
    )

    print("Running:", " ".join(command))
    subprocess.run(command, check=True)

    if not output_sdf.is_file() or output_sdf.stat().st_size == 0:
        raise RuntimeError(f"LigPrep fail: {output_sdf}")

    print(f"Prepared ligand: {output_sdf}")
    return output_sdf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare one ligand with LigPrep")
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    prepare_ligand(args.ligand, args.output)

