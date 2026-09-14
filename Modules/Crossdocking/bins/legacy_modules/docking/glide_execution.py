"""Execute one Glide input file and wait for completion."""

import argparse
import os
import subprocess
from pathlib import Path


def run_glide(glide_input, output_dir):
    glide_input = Path(glide_input).resolve()
    output_dir = Path(output_dir).resolve()

    if not glide_input.is_file():
        raise FileNotFoundError(f"Glide input not found: {glide_input}")

    schrodinger = os.environ.get("SCHRODINGER")
    if not schrodinger:
        raise EnvironmentError("crit SCHRODINGER environ is none")

    glide = Path(schrodinger) / "glide"
    if not glide.is_file():
        raise FileNotFoundError(f"Glide module not detec: {glide}")

    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "glide_console.log"
    command = [str(glide), "-WAIT", "-NOJOBID", "-NOLOCAL", str(glide_input)]

    print("Running:", " ".join(command))
    with log_file.open("w") as log:
        subprocess.run(
            command,
            cwd=output_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )

    pose_files = sorted(output_dir.glob("*_lib.sdfgz"))
    if not pose_files:
        raise RuntimeError(f"Glide produced no *_lib.sdfgz file; see {log_file}")

    print(f"Glide pose file: {pose_files[0]}")
    return pose_files[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a Glide docking input")
    parser.add_argument("--input", required=True, help="Glide .in file")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    run_glide(args.input, args.output_dir)
