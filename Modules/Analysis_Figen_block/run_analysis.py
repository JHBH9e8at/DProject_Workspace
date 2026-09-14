"""Run selected Analysis_Figen modules from one configuration file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKING_ROOT = HERE.parents[1]
if str(WORKING_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKING_ROOT))

from Modules.Analysis_Figen_block.common.config import parse_on_off, read_key_value_config, require_values  # noqa: E402
from Modules.Analysis_Figen_block.common.paths import resolve_config_path, validate_results_root  # noqa: E402
from Modules.Analysis_Figen_block.ahc_analysis.run_ahc_analysis import run as run_ahc  # noqa: E402
from Modules.Analysis_Figen_block.structural_analysis.run_structural_analysis import run as run_structural  # noqa: E402
from Modules.Analysis_Figen_block.interaction_result_analysis.run_interaction_result_analysis import run as run_interaction_results  # noqa: E402
from Modules.Analysis_Figen_block.figures.run_figures import run as run_figures  # noqa: E402

MODULES = ("ahc_analysis", "structural_analysis", "interaction_result_analysis", "figures")
RUNNERS = {
    "ahc_analysis": run_ahc,
    "structural_analysis": run_structural,
    "interaction_result_analysis": run_interaction_results,
    "figures": run_figures,
}
CONFIG_KEYS = {
    "ahc_analysis": "ahc_config",
    "structural_analysis": "structural_config",
    "interaction_result_analysis": "interaction_result_config",
    "figures": "figures_config",
}
ALLOWED_KEYS = {"results_root", "run_id", "dry_run", *MODULES, *CONFIG_KEYS.values()}


def parse_config(config_path: str | Path) -> dict[str, object]:
    document = read_key_value_config(config_path, allowed_keys=ALLOWED_KEYS)
    raw = document.values
    require_values(raw, ("results_root", "run_id"), source=document.path)
    cfg: dict[str, object] = {
        "config_path": document.path,
        "results_root": validate_results_root(resolve_config_path(raw["results_root"], document.path)),
        "run_id": raw["run_id"],
        "dry_run": parse_on_off(raw.get("dry_run", "off"), key="dry_run"),
    }
    for module in MODULES:
        enabled = parse_on_off(raw.get(module, "off"), key=module)
        cfg[module] = enabled
        key = CONFIG_KEYS[module]
        if enabled:
            require_values(raw, (key,), source=document.path)
            cfg[key] = resolve_config_path(raw[key], document.path)
        else:
            cfg[key] = resolve_config_path(raw[key], document.path) if raw.get(key) else None
    if not any(cfg[module] for module in MODULES):
        raise ValueError("At least one analysis module must be ON")
    return cfg


def build_plan(cfg: dict[str, object]) -> list[dict[str, object]]:
    return [
        {"module": module, "config": cfg[CONFIG_KEYS[module]]}
        for module in MODULES if cfg[module]
    ]


def run(config_path: str | Path) -> dict[str, object]:
    cfg = parse_config(config_path)
    plan = build_plan(cfg)
    print(f"Analysis run_id: {cfg['run_id']}")
    print(f"Results root: {cfg['results_root']}")
    for index, item in enumerate(plan, start=1):
        print(f"[{index}/{len(plan)}] {item['module']}: {item['config']}")
    if cfg["dry_run"]:
        return {"status": "dry_run", "plan": plan, "completed": []}
    completed: list[str] = []
    results: dict[str, object] = {}
    for item in plan:
        module = str(item["module"])
        results[module] = RUNNERS[module](
            item["config"], run_id=str(cfg["run_id"]), results_root=cfg["results_root"]
        )
        completed.append(module)
    return {"status": "completed", "plan": plan, "completed": completed, "results": results}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    run(parser.parse_args(argv).config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
