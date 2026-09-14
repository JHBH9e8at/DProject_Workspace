"""Single configuration-driven entry point for completed AHC analysis."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKING_ROOT = HERE.parents[2]
if str(WORKING_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKING_ROOT))

from Modules.ACHresultspackage.CompletedAnalysis.blocks import (  # noqa: E402
    chemical_space, last_window, physchem, scaffold,
)
from Modules.Analysis_Figen_block.common.paths import validate_results_root  # noqa: E402


DEFAULT_REFERENCE = WORKING_ROOT / "AHC_Related" / "ref_structures" / "ref_desc.csv"
REQUIRED = ("pr_input", "pps_input", "results_root", "run_id")
TOGGLES = ("physchem", "last_window", "scaffold", "chemical_space")


def _path(value: str, base: Path) -> Path:
    candidate = Path(value).expanduser()
    return (base / candidate).resolve(strict=False) if not candidate.is_absolute() else candidate.resolve(strict=False)


def _on_off(value: str, key: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"on", "off"}:
        raise ValueError(f"{key} must be on or off")
    return normalized == "on"


def parse_config(path: str | Path) -> dict[str, object]:
    config_path = Path(path).expanduser().resolve(strict=False)
    raw: dict[str, str] = {}
    with config_path.open("r", encoding="utf-8") as handle:
        for lineno, source in enumerate(handle, start=1):
            line = source.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(f"{config_path}:{lineno}: expected key=value")
            key, value = line.split("=", 1)
            key = key.strip().lower()
            if key in raw:
                raise ValueError(f"{config_path}:{lineno}: duplicate key {key!r}")
            raw[key] = value.strip()
    missing = [key for key in REQUIRED if not raw.get(key)]
    if missing:
        raise ValueError(f"{config_path}: missing required keys: {missing}")

    base = config_path.parent
    cfg: dict[str, object] = dict(raw)
    for key in ("pr_input", "pps_input", "results_root"):
        cfg[key] = _path(raw[key], base)
    cfg["results_root"] = validate_results_root(cfg["results_root"])
    cfg["reference"] = _path(raw.get("reference", str(DEFAULT_REFERENCE)), base)
    cfg["last_n"] = int(raw.get("last_n", "50"))
    cfg["random_state"] = int(raw.get("random_state", "42"))
    cfg["resume"] = _on_off(raw.get("resume", "off"), "resume")
    cfg["usecache"] = _on_off(raw.get("usecache", "on"), "usecache")
    defaults = {"physchem": "on", "last_window": "on", "scaffold": "off", "chemical_space": "off"}
    for key in TOGGLES:
        cfg[key] = _on_off(raw.get(key, defaults[key]), key)
    cfg["methods"] = tuple(item.strip().lower() for item in raw.get("methods", "umap,tsne").split(",") if item.strip())
    cfg["features"] = tuple(item.strip().lower() for item in raw.get("features", "fp,all_desc,select_desc").split(",") if item.strip())
    if not cfg["methods"] or set(cfg["methods"]) - {"umap", "tsne"}:
        raise ValueError("methods must contain umap and/or tsne")
    if not cfg["features"] or set(cfg["features"]) - {"fp", "all_desc", "select_desc"}:
        raise ValueError("features contains an unsupported value")
    if cfg["last_n"] < 1:
        raise ValueError("last_n must be at least 1")
    cfg["config_path"] = config_path
    return cfg


def run(
    config_path: str | Path,
    *,
    run_id: str | None = None,
    results_root: str | Path | None = None,
) -> dict[str, object]:
    cfg = parse_config(config_path)
    if run_id is not None:
        cfg["run_id"] = run_id
    if results_root is not None:
        cfg["results_root"] = validate_results_root(results_root)
    common = dict(
        pr_input=Path(cfg["pr_input"]), pps_input=Path(cfg["pps_input"]),
        results_root=Path(cfg["results_root"]), run_id=str(cfg["run_id"]),
        resume=bool(cfg["resume"]),
    )
    completed: list[str] = []
    internal_dependencies: list[str] = []
    if cfg["physchem"]:
        physchem.run(**common, reference=Path(cfg["reference"]), last_n=int(cfg["last_n"]))
        completed.append("physchem")
    if cfg["last_window"]:
        last_window.run(**common, last_n=int(cfg["last_n"]))
        completed.append("last_window")
    membership = None
    if cfg["chemical_space"] or cfg["scaffold"]:
        methods = tuple(cfg["methods"]) if cfg["chemical_space"] else ("umap",)
        features = tuple(cfg["features"]) if cfg["chemical_space"] else ("fp",)
        if cfg["scaffold"]:
            methods = ("umap",) + tuple(item for item in methods if item != "umap")
            features = ("fp",) + tuple(item for item in features if item != "fp")
        chemical_result = chemical_space.run(
            **common, reference=Path(cfg["reference"]), methods=methods,
            features=features, random_state=int(cfg["random_state"]),
            render_figures=bool(cfg["chemical_space"]),
            usecache=bool(cfg["usecache"]),
        )
        membership = chemical_result["membership"]
        if cfg["chemical_space"]:
            completed.append("chemical_space")
        else:
            internal_dependencies.append("chemical_space:umap/fp membership")
    if cfg["scaffold"]:
        if membership is None:
            raise RuntimeError("chemical-space calculation did not produce joint UMAP/fingerprint membership")
        scaffold.run(**common, membership=Path(membership), last_n=int(cfg["last_n"]))
        completed.append("scaffold")

    summary_dir = Path(cfg["results_root"]) / "analysis" / "completed_ahc" / str(cfg["run_id"])
    summary_dir.mkdir(parents=True, exist_ok=bool(cfg["resume"]))
    summary = {
        "status": "completed", "run_id": cfg["run_id"], "blocks": completed,
        "internal_dependencies": internal_dependencies,
        "chemical_space_cache_hit": (
            bool(chemical_result["cache_hit"])
            if cfg["chemical_space"] or cfg["scaffold"] else None
        ),
        "config": str(cfg["config_path"]),
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    (summary_dir / "completed_analysis.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8",
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run all enabled completed AHC analysis blocks")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args(argv)
    run(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
