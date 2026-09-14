"""Unpack a gzip-compressed SDF file without changing its contents."""

import argparse
import gzip
import shutil
from pathlib import Path


def unpack_sdfgz(input_file, output_file=None):
    input_file = Path(input_file)

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if input_file.suffix.lower() != ".sdfgz":
        raise ValueError(f"Expected an .sdfgz file: {input_file}")

    if output_file is None:
        output_file = input_file.with_suffix(".sdf")
    else:
        output_file = Path(output_file)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(input_file, "rb") as compressed:
        with output_file.open("wb") as uncompressed:
            shutil.copyfileobj(compressed, uncompressed)

    sdf_text = output_file.read_text(errors="replace")
    molecule_count = sum(line.strip() == "$$$$" for line in sdf_text.splitlines())

    if molecule_count == 0:
        output_file.unlink()
        raise ValueError("The decompressed file contains no complete SDF records")

    print(f"Input: {input_file.resolve()}")
    print(f"Output: {output_file.resolve()}")
    print(f"Molecule records: {molecule_count}")
    return output_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unpack .sdfgz to .sdf")
    parser.add_argument("input", help="Input .sdfgz file")
    parser.add_argument("-o", "--output", help="Output .sdf path")
    args = parser.parse_args()

    unpack_sdfgz(args.input, args.output)

