"""Legacy configuration for the retired interaction execution adapter."""

from __future__ import annotations

from pathlib import Path


SCRIPT_ID = "PY-119"
REQUIRED = {"pps_results", "pr_results", "pps_run_dir", "pr_run_dir",
            "pps_receptor", "pr_receptor"}


def parse_int_assignments(values: list[str], *, label: str) -> dict[str, int]:
    """Parse repeatable NAME=POSITIVE_INTEGER CLI values."""

    parsed = {}
    for raw in values:
        if "=" not in raw:
            raise ValueError(f"{label} must use NAME=VALUE: {raw!r}")
        key, value = (part.strip() for part in raw.split("=", 1))
        if not key or key in parsed:
            raise ValueError(f"Invalid or duplicate {label} key: {key!r}")
        parsed[key] = int(value)
        if parsed[key] <= 0:
            raise ValueError(f"{label} values must be positive: {raw!r}")
    if not parsed:
        raise ValueError(f"At least one {label} is required")
    return parsed


def parse_path_assignments(values: list[str], *, label: str) -> dict[str, Path]:
    """Parse repeatable NAME=PATH CLI values."""
    parsed = {}
    for raw in values:
        if "=" not in raw: raise ValueError(f"{label} must use NAME=PATH: {raw!r}")
        key, value = (part.strip() for part in raw.split("=", 1))
        if not key or not value or key in parsed:
            raise ValueError(f"Invalid or duplicate {label}: {raw!r}")
        parsed[key] = Path(value)
    if not parsed: raise ValueError(f"At least one {label} is required")
    return parsed


def parse_interaction_config(path: str | Path) -> dict:
    path = Path(path); cfg = {}
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"): continue
        if "=" not in line: raise ValueError(f"{path}:{lineno}: expected key=value")
        key, value = line.split("=", 1); key = key.strip().lower()
        if key in cfg: raise ValueError(f"{path}:{lineno}: duplicate key: {key}")
        cfg[key] = value.strip()
    missing = sorted(key for key in REQUIRED if not cfg.get(key))
    if missing: raise ValueError(f"{path}: missing required key(s): {missing}")
    defaults = {"pps_run_name": "PPS", "pr_run_name": "PR", "step_col": "step",
        "pps_variant_col": "PPS_best_variant", "pr_variant_col": "PR_best_variant",
        "smiles_col": "smiles", "analysis": "both"}
    for key, value in defaults.items(): cfg[key] = cfg.get(key, "").strip() or value
    cfg["analysis"] = cfg["analysis"].lower()
    if cfg["analysis"] not in {"posecheck", "prolif", "both"}: raise ValueError("invalid analysis")
    cfg["execution_mode"] = cfg.get("execution_mode", "mp").strip().lower() or "mp"
    if cfg["execution_mode"] not in {"standard", "mp"}:
        raise ValueError("execution_mode must be standard or mp")
    for key in ("include_secondary_interactions", "allow_unavailable", "resume", "fail_fast"):
        value = cfg.get(key, "off").lower()
        if value not in {"on", "off"}: raise ValueError(f"{key} must be on or off")
        cfg[key] = value == "on"
    cfg["posecheck_workers"] = int(cfg.get("posecheck_workers", "30") or 30)
    cfg["posecheck_chunk_size"] = int(cfg.get("posecheck_chunk_size", "250") or 250)
    workers = cfg.get("prolif_workers", "30").strip(); cfg["prolif_workers"] = int(workers) if workers else None
    if cfg["posecheck_workers"] <= 0 or cfg["posecheck_chunk_size"] <= 0 or (
        cfg["prolif_workers"] is not None and cfg["prolif_workers"] <= 0): raise ValueError("worker settings must be positive")
    top_n = cfg.get("testmode_top_n", "").strip(); cfg["testmode_top_n"] = int(top_n) if top_n else None
    cfg["residue_map"] = cfg.get("residue_map", "").strip() or None
    for key in ("max_clashes", "max_clashes_per_heavy_atom", "max_strain_energy"):
        value = cfg.get(key, "").strip(); cfg[key] = float(value) if value else None
    cfg.pop("output_dir", None)
    return cfg
