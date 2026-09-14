"""Run serial or multiprocessing analysis over best cross-docked poses."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd


ARMS = (
    {
        "arm": "L_PPS_to_R_PR",
        "ligand_state": "PPS",
        "receptor_state": "PR",
        "receptor_key": "pr_receptor",
    },
    {
        "arm": "L_PR_to_R_PPS",
        "ligand_state": "PR",
        "receptor_state": "PPS",
        "receptor_key": "pps_receptor",
    },
)
VALID_ANALYSIS = {"posecheck", "prolif", "both"}
VALID_ONOFF = {"on", "off"}


def parse_config(path: str | Path) -> dict:
    path = Path(path)
    cfg: dict[str, object] = {}
    with path.open(encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(f"{path}:{lineno}: expected key=value")
            key, value = line.split("=", 1)
            key = key.strip().lower()
            if key in cfg:
                raise ValueError(f"{path}:{lineno}: duplicate key: {key}")
            cfg[key] = value.strip()

    required = ("crossdocking_dir", "pps_receptor", "pr_receptor", "output_dir")
    missing = [key for key in required if not cfg.get(key)]
    if missing:
        raise ValueError(f"{path}: missing required key(s): {missing}")

    cfg["analysis"] = str(cfg.get("analysis", "both")).lower()
    if cfg["analysis"] not in VALID_ANALYSIS:
        raise ValueError(f"analysis must be one of {sorted(VALID_ANALYSIS)}")
    for key in ("include_secondary_interactions", "resume", "fail_fast"):
        value = str(cfg.get(key, "off")).lower()
        if value not in VALID_ONOFF:
            raise ValueError(f"{key} must be on or off")
        cfg[key] = value == "on"

    top_n = str(cfg.get("testmode_top_n", "")).strip()
    cfg["testmode_top_n"] = int(top_n) if top_n else None
    if cfg["testmode_top_n"] is not None and cfg["testmode_top_n"] <= 0:
        raise ValueError("testmode_top_n must be positive")
    cfg["posecheck_workers"] = int(str(cfg.get("posecheck_workers", "30")))
    chunk_size = str(cfg.get("posecheck_chunk_size", "auto")).strip().lower()
    cfg["posecheck_chunk_size"] = "auto" if chunk_size == "auto" else int(chunk_size)
    prolif = str(cfg.get("prolif_workers", "")).strip()
    cfg["prolif_workers"] = int(prolif) if prolif else 30
    if cfg["posecheck_workers"] <= 0:
        raise ValueError("posecheck_workers must be positive")
    if cfg["posecheck_chunk_size"] != "auto" and cfg["posecheck_chunk_size"] <= 0:
        raise ValueError("posecheck_chunk_size must be auto or a positive integer")
    if cfg["prolif_workers"] <= 0:
        raise ValueError("prolif_workers must be positive")
    cfg["residue_map"] = str(cfg.get("residue_map", "")).strip() or None
    return cfg


def build_crossdock_manifest(crossdocking_dir: str | Path, arm_spec: dict,
                             receptor_file: str | Path,
                             testmode_top_n: int | None = None) -> pd.DataFrame:
    """Resolve completed cross-docking jobs and their best 1-based pose records."""
    arm_dir = Path(crossdocking_dir) / arm_spec["arm"]
    manifest_path = arm_dir / "target_manifest.csv"
    table = pd.read_csv(manifest_path)
    required = {"molecule_id", "pipeline_status", "best_docking_score", "job_dir"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{manifest_path}: missing columns: {missing}")

    table = table.loc[table["pipeline_status"].eq("completed")].copy()
    table["best_docking_score"] = pd.to_numeric(
        table["best_docking_score"], errors="coerce"
    )
    table = table.dropna(subset=["best_docking_score"])
    table = table.sort_values("best_docking_score", ascending=True)
    if testmode_top_n is not None:
        table = table.head(testmode_top_n).copy()

    records = []
    for rank, (_, row) in enumerate(table.iterrows(), start=1):
        job_name = Path(str(row["job_dir"])).name
        job_dir = arm_dir / "jobs" / job_name
        pose_file = job_dir / "docking_lib.sdfgz"
        best_csv = job_dir / "best_docking_score.csv"
        best = pd.read_csv(best_csv)
        if len(best) != 1 or "pose" not in best.columns:
            raise ValueError(f"Expected one best pose row with a pose column: {best_csv}")
        best_pose = int(best.iloc[0]["pose"])
        records.append({
            "interaction_selection_rank": rank,
            "crossdock_arm": arm_spec["arm"],
            "ligand_state": arm_spec["ligand_state"],
            "receptor_state": arm_spec["receptor_state"],
            "molecule_id": str(row["molecule_id"]),
            "best_docking_score": float(row["best_docking_score"]),
            "best_pose_record": best_pose,
            "receptor_file": str(Path(receptor_file)),
            "source_pose_file": str(pose_file),
            "pose_file_exists": pose_file.is_file(),
            "best_score_file": str(best_csv),
        })
    manifest = pd.DataFrame(records)
    if manifest.empty:
        raise ValueError(f"No completed cross-docking jobs found: {manifest_path}")
    return manifest


def _read_sdf_records(path: Path) -> list[bytes]:
    with gzip.open(path, "rb") as handle:
        payload = handle.read()
    chunks = payload.split(b"$$$$")
    return [chunk + b"$$$$\n" for chunk in chunks if chunk.strip()]


def materialize_best_population(manifest: pd.DataFrame, output_sdf: str | Path,
                                source_map_csv: str | Path) -> pd.DataFrame:
    """Extract one best pose record per manifest row into a population SDF."""
    output_sdf = Path(output_sdf)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    source_rows = []
    with output_sdf.open("wb") as destination:
        for pose_index, (_, row) in enumerate(manifest.iterrows()):
            source_file = Path(row["source_pose_file"])
            if not source_file.is_file():
                raise FileNotFoundError(source_file)
            records = _read_sdf_records(source_file)
            record_number = int(row["best_pose_record"])
            if record_number < 1 or record_number > len(records):
                raise IndexError(
                    f"Best pose {record_number} outside 1..{len(records)}: {source_file}"
                )
            destination.write(records[record_number - 1])
            source_rows.append({
                "pose_index": pose_index,
                "crossdock_arm": row["crossdock_arm"],
                "ligand_state": row["ligand_state"],
                "receptor_state": row["receptor_state"],
                "molecule_id": row["molecule_id"],
                "best_docking_score": row["best_docking_score"],
                "source_pose_file": row["source_pose_file"],
                "source_sdf_record": record_number,
            })
    source_map = pd.DataFrame(source_rows)
    source_map.to_csv(source_map_csv, index=False)
    return source_map


def _collect_outputs(job_dir: Path, source_map: pd.DataFrame,
                     aggregates: dict[str, list[pd.DataFrame]]) -> None:
    for key, filename in (
        ("pose_metadata", "pose_metadata.csv"),
        ("pose_quality", "pose_quality.csv"),
        ("interactions_long", "interactions_long.csv"),
        ("run_summary", "run_summary.csv"),
    ):
        path = job_dir / filename
        if not path.is_file():
            continue
        table = pd.read_csv(path)
        if "pose_index" in table.columns:
            annotations = source_map[
                [column for column in source_map.columns
                 if column == "pose_index" or column not in table.columns]
            ]
            table = table.merge(annotations, on="pose_index", how="left",
                                validate="many_to_one")
        else:
            arm = source_map.iloc[0]
            for column in ("crossdock_arm", "ligand_state", "receptor_state"):
                if column not in table.columns:
                    table.insert(0, column, arm[column])
        aggregates[key].append(table)


def run_crossdock_interactions(cfg: dict, analyzer=None) -> None:
    if analyzer is None:
        try:
            from .common.single_complex import run_single_complex as analyzer
        except ImportError:
            from common.single_complex import run_single_complex as analyzer

    output_dir = Path(str(cfg["output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "manifests").mkdir(exist_ok=True)
    aggregates: dict[str, list[pd.DataFrame]] = {
        "pose_metadata": [], "pose_quality": [],
        "interactions_long": [], "run_summary": [],
    }
    all_manifests = []

    for arm_spec in ARMS:
        receptor = cfg[arm_spec["receptor_key"]]
        manifest = build_crossdock_manifest(
            cfg["crossdocking_dir"], arm_spec, receptor, cfg["testmode_top_n"]
        )
        all_manifests.append(manifest)
        arm_dir = output_dir / "populations" / arm_spec["arm"]
        arm_dir.mkdir(parents=True, exist_ok=True)
        source_map = materialize_best_population(
            manifest, arm_dir / "best_crossdocked_poses.sdf",
            arm_dir / "source_pose_map.csv",
        )
        required = [arm_dir / "pose_metadata.csv", arm_dir / "run_summary.csv"]
        if cfg["analysis"] in {"posecheck", "both"}:
            required.append(arm_dir / "pose_quality.csv")
        if cfg["analysis"] in {"prolif", "both"}:
            required.append(arm_dir / "interactions_long.csv")
        if not (cfg["resume"] and all(path.is_file() for path in required)):
            analyzer(
                protein_file=receptor,
                pose_file=arm_dir / "best_crossdocked_poses.sdf",
                receptor_state=arm_spec["receptor_state"],
                output_dir=arm_dir,
                analysis=cfg["analysis"],
                residue_map=cfg["residue_map"],
                include_secondary_interactions=cfg["include_secondary_interactions"],
                posecheck_workers=cfg["posecheck_workers"],
                posecheck_chunk_size=cfg["posecheck_chunk_size"],
                prolif_workers=cfg["prolif_workers"],
                resume=cfg["resume"],
            )
        _collect_outputs(arm_dir, source_map, aggregates)

    combined_manifest = pd.concat(all_manifests, ignore_index=True)
    combined_manifest.to_csv(
        output_dir / "manifests" / "crossdock_interaction_manifest.csv", index=False
    )
    aggregate_dir = output_dir / "aggregate"
    aggregate_dir.mkdir(exist_ok=True)
    for key, tables in aggregates.items():
        if tables:
            pd.concat(tables, ignore_index=True).to_csv(
                aggregate_dir / f"{key}.csv", index=False
            )
    summary = {
        "status": "completed",
        "pose_scope": "one best cross-docked pose per molecule",
        "molecule_count": len(combined_manifest),
        "analysis": cfg["analysis"],
        "output_dir": str(output_dir),
    }
    (output_dir / "crossdock_interaction_run_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    run_crossdock_interactions(parse_config(args.config))


if __name__ == "__main__":
    main()
