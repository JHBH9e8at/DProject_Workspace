"""Configuration-driven entry point for completed cross-docking analysis."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from Modules.Analysis_Figen_block.common.paths import validate_results_root

from .interaction.crossdocked import run_crossdock_interactions
from .score.score_workflow import run_score_analysis


def _on_off(value: str, key: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"on", "off"}:
        raise ValueError(f"{key} must be on or off")
    return normalized == "on"


def _path(value: str, base: Path) -> Path:
    path = Path(value).expanduser()
    return (base / path).resolve(strict=False) if not path.is_absolute() else path.resolve(strict=False)


def parse_config(config_path: str | Path) -> dict[str, object]:
    path = Path(config_path).expanduser().resolve(strict=True)
    raw: dict[str, str] = {}
    for lineno, source in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = source.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{lineno}: expected key=value")
        key, value = line.split("=", 1)
        key = key.strip().lower()
        if key in raw:
            raise ValueError(f"{path}:{lineno}: duplicate key {key!r}")
        raw[key] = value.strip()

    missing = [key for key in ("results_root", "run_id") if not raw.get(key)]
    if missing:
        raise ValueError(f"{path}: missing required keys: {missing}")
    base = path.parent
    cfg: dict[str, object] = dict(raw)
    cfg["results_root"] = validate_results_root(_path(raw["results_root"], base))
    for key in ("calculation_manifest", "pps_input", "pr_input", "pps_receptor", "pr_receptor"):
        cfg[key] = _path(raw[key], base) if raw.get(key) else None
    if cfg["calculation_manifest"] is None and raw.get("output_dir"):
        cfg["calculation_manifest"] = _path(raw["output_dir"], base) / "crossdocking_calculation_manifest.json"

    cfg["score_analysis"] = _on_off(raw.get("score_analysis", "on"), "score_analysis")
    cfg["interaction_analysis"] = _on_off(raw.get("interaction_analysis", "off"), "interaction_analysis")
    cfg["multiprocessing"] = _on_off(raw.get("multiprocessing", "off"), "multiprocessing")
    cfg["resume"] = _on_off(raw.get("resume", "off"), "resume")
    cfg["include_secondary_interactions"] = _on_off(
        raw.get("include_secondary_interactions", "off"), "include_secondary_interactions"
    )
    cfg["posecheck_workers"] = int(raw.get("posecheck_workers", "30"))
    cfg["prolif_workers"] = int(raw.get("prolif_workers", "30"))
    chunk = raw.get("posecheck_chunk_size", "auto").strip().lower()
    cfg["posecheck_chunk_size"] = "auto" if chunk == "auto" else int(chunk)
    cfg["testmode_top_n"] = int(raw["testmode_top_n"]) if raw.get("testmode_top_n") else None
    cfg["selectivity_threshold"] = float(raw.get("selectivity_threshold", "2.0"))
    cfg["strong_score_threshold"] = float(raw.get("strong_score_threshold", "-8.0"))
    cfg["analysis"] = raw.get("interaction_mode", "both").strip().lower()
    if cfg["analysis"] not in {"posecheck", "prolif", "both"}:
        raise ValueError("interaction_mode must be posecheck, prolif, or both")
    if cfg["posecheck_workers"] <= 0 or cfg["prolif_workers"] <= 0:
        raise ValueError("worker counts must be positive")
    if cfg["posecheck_chunk_size"] != "auto" and cfg["posecheck_chunk_size"] <= 0:
        raise ValueError("posecheck_chunk_size must be auto or a positive integer")
    if (cfg["score_analysis"] or cfg["interaction_analysis"]) and cfg["calculation_manifest"] is None:
        raise ValueError("analysis requires calculation_manifest or output_dir")
    if cfg["interaction_analysis"] and (cfg["pps_receptor"] is None or cfg["pr_receptor"] is None):
        raise ValueError("interaction_analysis=on requires pps_receptor and pr_receptor")
    cfg["config_path"] = path
    return cfg


def run(config_path: str | Path) -> dict[str, object]:
    cfg = parse_config(config_path)
    completed: list[str] = []
    if cfg["score_analysis"]:
        run_score_analysis(
            calculation_manifest=cfg["calculation_manifest"],
            run_id=str(cfg["run_id"]), results_root=cfg["results_root"],
            resume=bool(cfg["resume"]), pps_scores=cfg["pps_input"],
            pr_scores=cfg["pr_input"],
            selectivity_threshold=float(cfg["selectivity_threshold"]),
            strong_score_threshold=float(cfg["strong_score_threshold"]),
        )
        completed.append("score")

    if cfg["interaction_analysis"]:
        payload = json.loads(Path(cfg["calculation_manifest"]).read_text(encoding="utf-8"))
        workers = bool(cfg["multiprocessing"])
        interaction_cfg = {
            "crossdocking_dir": payload["double_arm_dir"],
            "pps_receptor": cfg["pps_receptor"], "pr_receptor": cfg["pr_receptor"],
            "output_dir": Path(cfg["results_root"]) / "analysis" / "crossdocking_interactions" / str(cfg["run_id"]),
            "testmode_top_n": cfg["testmode_top_n"], "analysis": cfg["analysis"],
            "include_secondary_interactions": cfg["include_secondary_interactions"],
            "resume": cfg["resume"], "fail_fast": True, "residue_map": None,
            "posecheck_workers": cfg["posecheck_workers"] if workers else 1,
            "posecheck_chunk_size": cfg["posecheck_chunk_size"],
            "prolif_workers": cfg["prolif_workers"] if workers else 1,
        }
        run_crossdock_interactions(interaction_cfg)
        completed.append("interaction_mp" if workers else "interaction_standard")

    summary = {
        "status": "completed", "run_id": cfg["run_id"], "blocks": completed,
        "multiprocessing": bool(cfg["multiprocessing"]),
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    summary_dir = Path(cfg["results_root"]) / "analysis" / "completed_crossdocking" / str(cfg["run_id"])
    summary_dir.mkdir(parents=True, exist_ok=True)
    (summary_dir / "crossdocking_analysis.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args(argv)
    run(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
