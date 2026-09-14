"""Configuration-driven entry point for standalone figure regeneration."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKING_ROOT = HERE.parents[2]
if str(WORKING_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKING_ROOT))

from Modules.Analysis_Figen_block.common.config import parse_float, parse_int, parse_on_off, read_key_value_config, require_values  # noqa: E402
from Modules.Analysis_Figen_block.common.paths import resolve_config_path, validate_results_root  # noqa: E402
from Modules.Analysis_Figen_block.interaction_result_analysis.blocks.config_helpers import parse_int_assignments, parse_path_assignments  # noqa: E402
from Modules.Analysis_Figen_block.figures.blocks.ahc_physchem_trajectory_workflow import run_workflow  # noqa: E402
from Modules.Analysis_Figen_block.figures.blocks.chemical_space_workflow import run_chemical_space_figures  # noqa: E402
from Modules.Analysis_Figen_block.figures.blocks.physchem_last_window_workflow import run_last_window  # noqa: E402
from Modules.Analysis_Figen_block.figures.blocks.hotspot_variants_workflow import run_full_last_figure, run_top_candidate_figure  # noqa: E402
from Modules.Analysis_Figen_block.figures.blocks.interaction_hotspot_workflow import run_interaction_hotspot  # noqa: E402
from Modules.Analysis_Figen_block.figures.blocks.own_opposite_hotspot_workflow import run_figure as run_own_opposite_figure  # noqa: E402

TOGGLES = ("ahc_physchem", "chemical_space", "physchem_last_window", "interaction_hotspot", "hotspot_full_last", "hotspot_top_candidates", "own_opposite")
ALLOWED_KEYS = {
    "results_root", "run_id", "resume", *TOGGLES,
    "ahc_input", "ahc_job_name", "ahc_threshold", "reference", "ref_keys", "descriptors", "mw_upper_bound",
    "cache_path", "pr_input", "pps_input", "last_n", "pr_threshold", "pps_threshold",
    "interaction_residue_csv", "interaction_typed_csv", "interaction_pocket_files", "interaction_reference_csv",
    "full_csv", "last_csv", "typed_csv", "pocket_files", "reference_csv", "denominators",
    "top_own_csv", "top_opposite_csv", "top_typed_csv", "top_pocket_files", "top_reference_csv", "top_denominators",
    "ownopp_prevalence_csv", "ownopp_typed_csv", "ownopp_pocket_files", "ownopp_reference_csv", "ownopp_denominators",
}


def _int_map(value: str, key: str) -> dict[str, int]:
    return parse_int_assignments([item.strip() for item in value.split(",") if item.strip()], label=key)


def _path_map(value: str, key: str, config_path: Path) -> dict[str, Path]:
    raw = parse_path_assignments([item.strip() for item in value.split(",") if item.strip()], label=key)
    return {name: resolve_config_path(path, config_path) for name, path in raw.items()}


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
        raise ValueError("At least one figure block must be ON")
    requirements = {
        "ahc_physchem": ("ahc_input", "ahc_job_name", "ahc_threshold"),
        "chemical_space": ("cache_path",),
        "physchem_last_window": ("pr_input", "pps_input"),
        "interaction_hotspot": ("interaction_residue_csv", "interaction_typed_csv", "interaction_pocket_files", "interaction_reference_csv"),
        "hotspot_full_last": ("full_csv", "last_csv", "typed_csv", "pocket_files", "reference_csv", "denominators"),
        "hotspot_top_candidates": ("top_own_csv", "top_opposite_csv", "top_typed_csv", "top_pocket_files", "top_reference_csv", "top_denominators"),
        "own_opposite": ("ownopp_prevalence_csv", "ownopp_typed_csv", "ownopp_pocket_files", "ownopp_reference_csv", "ownopp_denominators"),
    }
    for block, keys in requirements.items():
        if cfg[block]:
            require_values(raw, keys, source=document.path)
    path_keys = ("ahc_input", "reference", "cache_path", "pr_input", "pps_input", "interaction_residue_csv", "interaction_typed_csv", "interaction_reference_csv", "full_csv", "last_csv", "typed_csv", "reference_csv", "top_own_csv", "top_opposite_csv", "top_typed_csv", "top_reference_csv", "ownopp_prevalence_csv", "ownopp_typed_csv", "ownopp_reference_csv")
    for key in path_keys:
        cfg[key] = resolve_config_path(raw[key], document.path) if raw.get(key) else None
    cfg.update({
        "ahc_job_name": raw.get("ahc_job_name"),
        "ahc_threshold": parse_float(raw["ahc_threshold"], key="ahc_threshold") if raw.get("ahc_threshold") else None,
        "ref_keys": [item.strip() for item in raw.get("ref_keys", "OM").split(",") if item.strip()],
        "descriptors": [item.strip() for item in raw.get("descriptors", "").split(",") if item.strip()] or None,
        "mw_upper_bound": parse_float(raw["mw_upper_bound"], key="mw_upper_bound") if raw.get("mw_upper_bound") else None,
        "last_n": parse_int(raw.get("last_n", "50"), key="last_n", minimum=1),
        "pr_threshold": parse_float(raw["pr_threshold"], key="pr_threshold") if raw.get("pr_threshold") else None,
        "pps_threshold": parse_float(raw["pps_threshold"], key="pps_threshold") if raw.get("pps_threshold") else None,
    })
    for key in ("denominators", "top_denominators", "ownopp_denominators"):
        cfg[key] = _int_map(raw[key], key) if raw.get(key) else None
    for key in ("interaction_pocket_files", "pocket_files", "top_pocket_files", "ownopp_pocket_files"):
        cfg[key] = _path_map(raw[key], key, document.path) if raw.get(key) else None
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
    if cfg["ahc_physchem"]:
        outputs["ahc_physchem"] = run_workflow(input_csv=cfg["ahc_input"], job_name=cfg["ahc_job_name"], threshold=cfg["ahc_threshold"], run_id=f"{cfg['run_id']}_ahc_physchem", ref_csv=cfg["reference"], ref_keys=cfg["ref_keys"], descriptors=cfg["descriptors"], mw_upper_bound=cfg["mw_upper_bound"], **common)
        completed.append("ahc_physchem")
    if cfg["chemical_space"]:
        outputs["chemical_space"] = run_chemical_space_figures(cache_path=cfg["cache_path"], run_id=f"{cfg['run_id']}_chemical_space", **common)
        completed.append("chemical_space")
    if cfg["physchem_last_window"]:
        outputs["physchem_last_window"] = run_last_window(pr_csv=cfg["pr_input"], pps_csv=cfg["pps_input"], last_n=cfg["last_n"], pr_threshold=cfg["pr_threshold"], pps_threshold=cfg["pps_threshold"], run_id=f"{cfg['run_id']}_last_window", **common)
        completed.append("physchem_last_window")
    if cfg["interaction_hotspot"]:
        outputs["interaction_hotspot"] = run_interaction_hotspot(
            residue_csv=cfg["interaction_residue_csv"],
            typed_csv=cfg["interaction_typed_csv"],
            pocket_files=cfg["interaction_pocket_files"],
            reference_csv=cfg["interaction_reference_csv"],
            run_id=f"{cfg['run_id']}_interaction_hotspot", **common,
        )
        completed.append("interaction_hotspot")
    if cfg["hotspot_full_last"]:
        outputs["hotspot_full_last"] = run_full_last_figure(full_csv=cfg["full_csv"], last_csv=cfg["last_csv"], typed_csv=cfg["typed_csv"], pocket_files=cfg["pocket_files"], reference_csv=cfg["reference_csv"], denominators=cfg["denominators"], run_id=f"{cfg['run_id']}_hotspot_full_last", **common)
        completed.append("hotspot_full_last")
    if cfg["hotspot_top_candidates"]:
        outputs["hotspot_top_candidates"] = run_top_candidate_figure(own_csv=cfg["top_own_csv"], opposite_csv=cfg["top_opposite_csv"], typed_csv=cfg["top_typed_csv"], pocket_files=cfg["top_pocket_files"], reference_csv=cfg["top_reference_csv"], denominators=cfg["top_denominators"], run_id=f"{cfg['run_id']}_hotspot_top", **common)
        completed.append("hotspot_top_candidates")
    if cfg["own_opposite"]:
        outputs["own_opposite"] = run_own_opposite_figure(prevalence_csv=cfg["ownopp_prevalence_csv"], typed_csv=cfg["ownopp_typed_csv"], pocket_files=cfg["ownopp_pocket_files"], reference_csv=cfg["ownopp_reference_csv"], denominators=cfg["ownopp_denominators"], run_id=f"{cfg['run_id']}_own_opposite", **common)
        completed.append("own_opposite")
    return {"status": "completed", "blocks": completed, "outputs": outputs}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    result = run(parser.parse_args(argv).config)
    print(f"Completed figure blocks: {', '.join(result['blocks'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
