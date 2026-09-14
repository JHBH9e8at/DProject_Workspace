from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def index_metadata(path: Path, state_column: str) -> dict[tuple[str, str], list[dict[str, str]]]:
    output: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in read_csv(path):
        if row.get("pose_read_status") != "valid":
            continue
        key = (row[state_column], row["molecule_id"])
        output.setdefault(key, []).append(row)
    return output


def select_metadata(
    indexed: dict[tuple[str, str], list[dict[str, str]]],
    state: str,
    molecule_id: str,
    expected_variant: str,
    expected_score: str,
) -> dict[str, str]:
    rows = indexed.get((state, molecule_id), [])
    if not rows:
        raise RuntimeError(f"No valid metadata for {state}/{molecule_id}")
    exact = [row for row in rows if row.get("glide_variant") == expected_variant]
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        rows = exact
    target = float(expected_score)
    rows = sorted(rows, key=lambda row: abs(float(row["docking_score"]) - target))
    if abs(float(rows[0]["docking_score"]) - target) > 1e-4:
        raise RuntimeError(
            f"Metadata score mismatch for {state}/{molecule_id}: "
            f"expected {target}, nearest {rows[0]['docking_score']}"
        )
    return rows[0]


def open_pose_source(path: Path):
    if path.name.lower().endswith((".sdfgz", ".sdf.gz")):
        return gzip.open(path, "rb")
    return path.open("rb")


def extract_sdf_record(source: Path, record_number: int, destination: Path) -> None:
    if record_number < 1:
        raise ValueError(f"SDF record must be 1-based: {record_number}")
    current = 0
    record: list[bytes] = []
    with open_pose_source(source) as handle:
        for line in handle:
            record.append(line)
            if line.strip() == b"$$$$":
                current += 1
                if current == record_number:
                    destination.write_bytes(b"".join(record))
                    return
                record = []
    raise RuntimeError(f"Record {record_number} not found in {source}; records read: {current}")


def validate_source(path: Path, allowed_roots: list[Path]) -> None:
    resolved = path.resolve()
    if not any(resolved.is_relative_to(root.resolve()) for root in allowed_roots):
        raise RuntimeError(f"Pose source lies outside allowed roots: {resolved}")
    if not resolved.is_file():
        raise FileNotFoundError(resolved)


def main(args: argparse.Namespace) -> None:
    candidates = read_csv(args.summary)
    own_index = index_metadata(args.own_metadata, "population_state")
    opposite_index = index_metadata(args.opposite_metadata, "ligand_state")
    allowed_roots = [args.pps_raw_root, args.pr_raw_root, args.crossdock_root]
    manifest: list[dict[str, str]] = []

    for candidate in candidates:
        state = candidate["own_state"]
        rank = int(candidate["candidate_rank"])
        molecule_id = candidate["molecule_id"]
        own_meta = select_metadata(
            own_index, state, molecule_id,
            candidate["own_glide_variant"], candidate["own_state_score"],
        )
        opposite_meta = select_metadata(
            opposite_index, state, molecule_id,
            candidate["opposite_glide_variant"], candidate["opposite_state_score"],
        )

        for condition, metadata, score in (
            ("Own", own_meta, candidate["own_state_score"]),
            ("Opposite", opposite_meta, candidate["opposite_state_score"]),
        ):
            source = Path(metadata["source_pose_file"])
            validate_source(source, allowed_roots)
            record_number = int(metadata["source_sdf_record"])
            output_dir = args.output_dir / state / condition
            output_dir.mkdir(parents=True, exist_ok=True)
            destination = output_dir / f"rank_{rank:02d}_{molecule_id}.sdf"
            extract_sdf_record(source, record_number, destination)
            manifest.append({
                "own_state": state,
                "candidate_rank": str(rank),
                "molecule_id": molecule_id,
                "condition": condition,
                "receptor_state": metadata["receptor_state"],
                "docking_score": score,
                "glide_variant": metadata["glide_variant"],
                "source_pose_file": str(source),
                "source_sdf_record": str(record_number),
                "extracted_pose_file": str(destination),
            })

    expected = {(state, condition): 5 for state in ("PR", "PPS") for condition in ("Own", "Opposite")}
    actual = {
        key: sum(row["own_state"] == key[0] and row["condition"] == key[1] for row in manifest)
        for key in expected
    }
    if actual != expected:
        raise RuntimeError(f"Unexpected extracted-pose counts: {actual}")

    manifest_path = args.output_dir / "extracted_pose_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    (args.output_dir / "extraction_summary.json").write_text(
        json.dumps({f"{state}/{condition}": count for (state, condition), count in actual.items()}, indent=2),
        encoding="utf-8",
    )
    print("Extracted exactly one pose per candidate and receptor condition")
    for key, count in actual.items():
        print(f"{key[0]}/{key[1]}: {count}")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract exact own/opposite poses for top state-selective candidates.")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--own-metadata", type=Path, required=True)
    parser.add_argument("--opposite-metadata", type=Path, required=True)
    parser.add_argument("--pps-raw-root", type=Path, required=True)
    parser.add_argument("--pr-raw-root", type=Path, required=True)
    parser.add_argument("--crossdock-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    main(parser.parse_args())
