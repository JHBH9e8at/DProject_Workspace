from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def index_opposite_metadata(path: Path, selected: set[tuple[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    output: dict[tuple[str, str], dict[str, str]] = {}
    for row in read_csv(path):
        key = (row["ligand_state"], row["molecule_id"])
        if key in selected:
            output[key] = row
    return output


def find_job(jobs_root: Path, molecule_id: str) -> Path:
    matches = sorted(path for path in jobs_root.glob(f"*_{molecule_id}") if path.is_dir())
    if len(matches) != 1:
        raise RuntimeError(f"Expected one job for {molecule_id}, found {len(matches)}: {matches}")
    return matches[0]


def copy_required(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    shutil.copy2(source, destination)


def main(args: argparse.Namespace) -> None:
    candidates = read_csv(args.summary)
    selected = {(row["own_state"], row["molecule_id"]) for row in candidates}
    opposite_meta = index_opposite_metadata(args.opposite_metadata, selected)
    missing = selected - set(opposite_meta)
    if missing:
        raise RuntimeError(f"Missing opposite-state metadata: {sorted(missing)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str]] = []
    for row in candidates:
        state = row["own_state"]
        rank = int(row["candidate_rank"])
        molecule_id = row["molecule_id"]
        arm = row["crossdock_arm"]
        metadata = opposite_meta[(state, molecule_id)]
        job = find_job(args.crossdock_root / arm / "jobs", molecule_id)
        candidate_dir = args.output_dir / state / f"rank_{rank:02d}_{molecule_id}"
        candidate_dir.mkdir(parents=True, exist_ok=True)

        smi_source = job / "input" / f"{molecule_id}.smi"
        prepared_source = job / "ligprep" / f"{molecule_id}_prepared.sdf"
        docking_source = job / "docking_lib.sdfgz"
        copied = {
            "input_smiles": (smi_source, candidate_dir / f"{molecule_id}.smi"),
            "prepared_ligand": (prepared_source, candidate_dir / f"{molecule_id}_prepared.sdf"),
            "opposite_state_docking_library": (
                docking_source,
                candidate_dir / f"{molecule_id}_opposite_state_docking_lib.sdfgz",
            ),
        }
        for role, (source, destination) in copied.items():
            copy_required(source, destination)
            manifest_rows.append({
                "own_state": state,
                "candidate_rank": str(rank),
                "molecule_id": molecule_id,
                "crossdock_arm": arm,
                "file_role": role,
                "copied_file": str(destination),
                "source_file": str(source),
                "source_sdf_record": metadata.get("source_sdf_record", "") if role == "opposite_state_docking_library" else "",
                "opposite_glide_variant": metadata.get("glide_variant", "") if role == "opposite_state_docking_library" else "",
            })

        (candidate_dir / "candidate_metadata.json").write_text(
            json.dumps({**row, "opposite_pose_metadata": metadata}, indent=2),
            encoding="utf-8",
        )

    manifest_path = args.output_dir / "ligand_file_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Collected {len(candidates)} candidates and {len(manifest_rows)} ligand files")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect local ligand files for top state-selective candidates.")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--opposite-metadata", type=Path, required=True)
    parser.add_argument("--crossdock-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    main(parser.parse_args())
