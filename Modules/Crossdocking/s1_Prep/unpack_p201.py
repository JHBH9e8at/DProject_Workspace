"""Unpack a gzip-compressed SDF file without changing its contents."""

import argparse
import gzip
import logging
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)


def _setup_logging(log_file):
    logger.setLevel(logging.INFO)
    logger.propagate = False

    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)   # Close previous handlers to prev leak

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)         # Append mode preserves the history but probs useless.


def unpack_sdfgz(input_file, output_dir, log_dir=None):
    input_file = Path(input_file)

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if input_file.suffix.lower() != ".sdfgz":
        raise ValueError(f"Expected an .sdfgz file: {input_file}")

    output_dir = Path(output_dir)
    log_dir = output_dir if log_dir is None else Path(log_dir)

    output_file = output_dir / input_file.with_suffix(".sdf").name
    log_file = log_dir / input_file.with_suffix(".log").name

    log_dir.mkdir(parents=True, exist_ok=True)
    _setup_logging(log_file)

    if output_file.exists():
        logger.error("Output file already exists: %s", output_file)
        raise FileExistsError(
            f"Output file already exists: {output_file}. "
            "Specify a different -o/--output-dir."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Decompressing %s -> %s", input_file, output_file)
    try:
        with gzip.open(input_file, "rb") as compressed:
            with output_file.open("wb") as uncompressed:
                shutil.copyfileobj(compressed, uncompressed)
    except Exception:
        logger.exception("Decompression failed, removing partial output: %s", output_file)
        if output_file.exists():
            output_file.unlink()
        raise

    molecule_count = 0
    with output_file.open("r", errors="replace") as f:
        for line in f:
            if line.strip() == "$$$$":
                molecule_count += 1

    if molecule_count == 0:
        logger.error("No SDF records found, removing output: %s", output_file)
        output_file.unlink()
        raise ValueError("The decompressed file contains no complete SDF records")

    logger.info("Input: %s", input_file.resolve())
    logger.info("Output: %s", output_file.resolve())
    logger.info("Molecule records: %d", molecule_count)
    return output_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unpack .sdfgz to .sdf")
    parser.add_argument("input", help="Input .sdfgz file")
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory to write the decompressed .sdf file to",
    )
    parser.add_argument(
        "--log-dir",
        help="Directory to write the log file to (defaults to --output-dir)",
    )
    args = parser.parse_args()

    unpack_sdfgz(args.input, args.output_dir, args.log_dir)

