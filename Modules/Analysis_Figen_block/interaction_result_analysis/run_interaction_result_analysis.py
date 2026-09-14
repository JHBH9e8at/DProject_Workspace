"""Analyze existing interaction result files without running ProLIF or PoseCheck."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKING_ROOT = HERE.parents[2]
if str(WORKING_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKING_ROOT))

from Modules.Analysis_Figen_block.common.config import parse_on_off, read_key_value_config, require_values  # noqa: E402
from Modules.Analysis_Figen_block.common.paths import resolve_config_path, validate_results_root  # noqa: E402
from Modules.Analysis_Figen_block.interaction_result_analysis.blocks.config_helpers import parse_int_assignments  # noqa: E402
from Modules.Analysis_Figen_block.interaction_result_analysis.blocks.hotspot_variants_runner import run_full_last, run_top_candidates  # noqa: E402
from Modules.Analysis_Figen_block.interaction_result_analysis.blocks.own_opposite_hotspot_runner import run_own_opposite_hotspot  # noqa: E402
from Modules.Analysis_Figen_block.interaction_result_analysis.blocks.residue_hotspot_runner import run_residue_hotspot  # noqa: E402

TOGGLES = ("residue_hotspot", "full_last", "top_candidates", "own_opposite")
PATH_KEYS = (
    "interactions_csv", "metadata_csv", "summary_csv", "transitions_csv",
    "own_metadata", "own_interactions", "opposite_metadata", "opposite_interactions",
)
ALLOWED_KEYS = {"results_root", "run_id", "resume", "denominators", "total_poses", "last_start", *TOGGLES, *PATH_KEYS}


def _assignments(value: str, *, label: str) -> dict[str, int]:
    return parse_int_assignments(
        [item.strip() for item in value.split(",") if item.strip()], label=label
    )


def parse_config(config_path: str | Path) -> dict[str, object]:
    document = read_key_value_config(config_path, allowed_keys=ALLOWED_KEYS)
    raw = document.values
    require_values(raw, ("results_root", "run_id"), source=document.path)
    cfg: dict[str, object] = {
        "config_path": document.path,
        "results_root": validate_results_root(resolve_config_path(raw["results_root"], document.path)),
        "run_id": raw["run_id"],
        "resume": parse_on_off(raw.get("resume", "off"), key="resume"),
    }
    for key in TOGGLES:
        cfg[key] = parse_on_off(raw.get(key, "off"), key=key)
    if not any(cfg[key] for key in TOGGLES):
        raise ValueError("At least one interaction-result block must be ON")
    requirements = {
        "residue_hotspot": ("interactions_csv", "denominators"),
        "full_last": ("interactions_csv", "metadata_csv", "total_poses", "last_start"),
        "top_candidates": ("summary_csv", "transitions_csv"),
        "own_opposite": ("own_metadata", "own_interactions", "opposite_metadata", "opposite_interactions"),
    }
    for block, keys in requirements.items():
        if cfg[block]:
            require_values(raw, keys, source=document.path)
    for key in PATH_KEYS:
        cfg[key] = resolve_config_path(raw[key], document.path) if raw.get(key) else None
    for key in ("denominators", "total_poses", "last_start"):
        cfg[key] = _assignments(raw[key], label=key) if raw.get(key) else None
    return cfg


def run(config_path: str | Path, *, run_id: str | None = None,
        results_root: str | Path | None = None) -> dict[str, object]:
    cfg = parse_config(config_path)
    if run_id is not None:
        cfg["run_id"] = run_id
    if results_root is not None:
        cfg["results_root"] = validate_results_root(results_root)
    common = {"results_root": cfg["results_root"], "resume": cfg["resume"]}
    completed: list[str] = []
    outputs: dict[str, object] = {}
    if cfg["residue_hotspot"]:
        outputs["residue_hotspot"] = run_residue_hotspot(
            interactions_csv=cfg["interactions_csv"], denominators=cfg["denominators"],
            run_id=f"{cfg['run_id']}_residue_hotspot", **common,
        )
        completed.append("residue_hotspot")
    if cfg["full_last"]:
        outputs["full_last"] = run_full_last(
            interactions_csv=cfg["interactions_csv"], metadata_csv=cfg["metadata_csv"],
            total_poses=cfg["total_poses"], last_start=cfg["last_start"],
            run_id=f"{cfg['run_id']}_full_last", **common,
        )
        completed.append("full_last")
    if cfg["top_candidates"]:
        outputs["top_candidates"] = run_top_candidates(
            summary_csv=cfg["summary_csv"], transitions_csv=cfg["transitions_csv"],
            run_id=f"{cfg['run_id']}_top_candidates", **common,
        )
        completed.append("top_candidates")
    if cfg["own_opposite"]:
        outputs["own_opposite"] = run_own_opposite_hotspot(
            own_metadata=cfg["own_metadata"], own_interactions=cfg["own_interactions"],
            opposite_metadata=cfg["opposite_metadata"], opposite_interactions=cfg["opposite_interactions"],
            run_id=f"{cfg['run_id']}_own_opposite", **common,
        )
        completed.append("own_opposite")
    return {"status": "completed", "blocks": completed, "outputs": outputs}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    result = run(parser.parse_args(argv).config)
    print(f"Completed interaction-result blocks: {', '.join(result['blocks'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
