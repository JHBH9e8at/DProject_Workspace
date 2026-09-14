"""Configuration-driven entry point for structural result analyses."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKING_ROOT = HERE.parents[2]
if str(WORKING_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKING_ROOT))

from Modules.Analysis_Figen_block.common.config import (  # noqa: E402
    parse_float, parse_int, parse_on_off, read_key_value_config, require_values,
)
from Modules.Analysis_Figen_block.common.paths import resolve_config_path, validate_results_root  # noqa: E402
from Modules.Analysis_Figen_block.structural_analysis.blocks.pocket_volume import parse_targets  # noqa: E402
from Modules.Analysis_Figen_block.structural_analysis.blocks.pocket_volume_runner import run_pocket_volume  # noqa: E402
from Modules.Analysis_Figen_block.structural_analysis.blocks.redocking_workflow import run_redocking_validation  # noqa: E402

ALLOWED_KEYS = {
    "results_root", "run_id", "resume", "pocket_volume", "redocking",
    "fpocket_dir", "targets", "ligand", "ligand_resname", "distance_threshold",
    "all_spheres", "iterations", "seed", "chunk_size", "radius_offset",
    "output_prefix", "pps_docking_csv", "pps_rmsd_csv", "pr_docking_csv",
    "pr_rmsd_csv", "score_ylim", "rmsd_ylim",
}


def _pair(value: str, *, key: str) -> tuple[float, float]:
    parts = [item.strip() for item in value.split(",")]
    if len(parts) != 2:
        raise ValueError(f"{key} must contain two comma-separated numbers")
    result = (parse_float(parts[0], key=key), parse_float(parts[1], key=key))
    if result[0] >= result[1]:
        raise ValueError(f"{key} lower bound must be less than upper bound")
    return result


def parse_config(config_path: str | Path) -> dict[str, object]:
    document = read_key_value_config(config_path, allowed_keys=ALLOWED_KEYS)
    raw = document.values
    require_values(raw, ("results_root", "run_id"), source=document.path)
    pocket = parse_on_off(raw.get("pocket_volume", "off"), key="pocket_volume")
    redocking = parse_on_off(raw.get("redocking", "off"), key="redocking")
    if not pocket and not redocking:
        raise ValueError("At least one of pocket_volume or redocking must be ON")
    if pocket:
        require_values(raw, ("fpocket_dir", "targets", "ligand"), source=document.path)
    if redocking:
        require_values(raw, ("pps_docking_csv", "pps_rmsd_csv", "pr_docking_csv", "pr_rmsd_csv"), source=document.path)

    cfg: dict[str, object] = {
        "config_path": document.path,
        "results_root": validate_results_root(resolve_config_path(raw["results_root"], document.path)),
        "run_id": raw["run_id"],
        "resume": parse_on_off(raw.get("resume", "off"), key="resume"),
        "pocket_volume": pocket,
        "redocking": redocking,
    }
    for key in ("fpocket_dir", "ligand", "pps_docking_csv", "pps_rmsd_csv", "pr_docking_csv", "pr_rmsd_csv"):
        cfg[key] = resolve_config_path(raw[key], document.path) if raw.get(key) else None
    cfg.update({
        "targets": parse_targets(raw["targets"]) if pocket else [],
        "ligand_resname": raw.get("ligand_resname") or None,
        "distance_threshold": parse_float(raw.get("distance_threshold", "3.0"), key="distance_threshold"),
        "all_spheres": parse_on_off(raw.get("all_spheres", "off"), key="all_spheres"),
        "iterations": parse_int(raw.get("iterations", "1000000"), key="iterations", minimum=1),
        "seed": parse_int(raw.get("seed", "42"), key="seed"),
        "chunk_size": parse_int(raw.get("chunk_size", "100000"), key="chunk_size", minimum=1),
        "radius_offset": parse_float(raw.get("radius_offset", "-1.6"), key="radius_offset"),
        "output_prefix": raw.get("output_prefix", "ligand_local_pocket_volume"),
        "score_ylim": _pair(raw.get("score_ylim", "-12,-1"), key="score_ylim"),
        "rmsd_ylim": _pair(raw.get("rmsd_ylim", "0,2.5"), key="rmsd_ylim"),
    })
    return cfg


def run(config_path: str | Path, *, run_id: str | None = None,
        results_root: str | Path | None = None) -> dict[str, object]:
    cfg = parse_config(config_path)
    if run_id is not None:
        cfg["run_id"] = run_id
    if results_root is not None:
        cfg["results_root"] = validate_results_root(results_root)
    completed: list[str] = []
    outputs: dict[str, object] = {}
    if cfg["pocket_volume"]:
        outputs["pocket_volume"] = run_pocket_volume(
            fpocket_dir=cfg["fpocket_dir"], targets=cfg["targets"], ligand=cfg["ligand"],
            run_id=f"{cfg['run_id']}_pocket_volume", results_root=cfg["results_root"],
            resume=cfg["resume"], output_prefix=cfg["output_prefix"],
            ligand_resname=cfg["ligand_resname"], distance_threshold=cfg["distance_threshold"],
            all_spheres=cfg["all_spheres"], iterations=cfg["iterations"], seed=cfg["seed"],
            chunk_size=cfg["chunk_size"], radius_offset=cfg["radius_offset"],
        )
        completed.append("pocket_volume")
    if cfg["redocking"]:
        outputs["redocking"] = run_redocking_validation(
            pps_docking_csv=cfg["pps_docking_csv"], pps_rmsd_csv=cfg["pps_rmsd_csv"],
            pr_docking_csv=cfg["pr_docking_csv"], pr_rmsd_csv=cfg["pr_rmsd_csv"],
            run_id=f"{cfg['run_id']}_redocking", results_root=cfg["results_root"],
            resume=cfg["resume"], score_ylim=cfg["score_ylim"], rmsd_ylim=cfg["rmsd_ylim"],
        )
        completed.append("redocking")
    return {"status": "completed", "blocks": completed, "outputs": outputs}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    result = run(parser.parse_args(argv).config)
    print(f"Completed structural blocks: {', '.join(result['blocks'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
