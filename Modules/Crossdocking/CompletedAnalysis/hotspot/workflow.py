"""Route interaction analysis into an isolated work_results run."""

from __future__ import annotations

import argparse
from pathlib import Path

from Modules.Analysis_Figen_block.common.output_naming import next_available_path
from Modules.Analysis_Figen_block.common.paths import validate_results_root
from Modules.Analysis_Figen_block.common.provenance import update_manifest
from Modules.Analysis_Figen_block.common.run_context import create_run_context, generate_run_id
from .config import parse_interaction_config


SCRIPT_ID = "PY-120"


def run_interaction_workflow(*, config_path: str | Path, run_id: str,
                             results_root: str | Path | None = None,
                             resume: bool = False, runner=None,
                             upstream_run_ids: tuple[str, ...] = ()):
    cfg = parse_interaction_config(config_path)
    context = create_run_context(("analysis", "interactions"), run_id,
                                 results_root=results_root, resume=resume)
    products = next_available_path(context.run_dir, "products"); products.mkdir()
    execution_mode = cfg["execution_mode"]
    if runner is None:
        if execution_mode == "standard":
            from Modules.Crossdocking.CompletedAnalysis.interaction.standard import run_batch_interactions
        else:
            from Modules.Crossdocking.CompletedAnalysis.interaction.multiprocessing import run_batch_interactions
        runner = run_batch_interactions
    call = {key: value for key, value in cfg.items() if key != "execution_mode"}
    if execution_mode == "standard":
        for key in ("posecheck_workers", "posecheck_chunk_size", "prolif_workers"):
            call.pop(key, None)
    call.update(output_dir=products, resume=False)
    inputs = [Path(config_path), Path(cfg["pps_results"]), Path(cfg["pr_results"]),
              Path(cfg["pps_receptor"]), Path(cfg["pr_receptor"])]
    try:
        result = runner(**call)
        outputs = sorted(path for path in products.rglob("*") if path.is_file())
        update_manifest(context, script_id=SCRIPT_ID, status="completed",
                        config={key: value for key, value in cfg.items() if key not in {"pps_results", "pr_results"}},
                        inputs=inputs, outputs=outputs,
                        upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config=cfg,
                        upstream_run_ids=upstream_run_ids); raise
    return context, products, result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run standard or MP interaction analysis.")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--run-id")
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--upstream-run-id", action="append", default=[])
    args = parser.parse_args(argv)
    root = validate_results_root(args.results_root) if args.results_root else None
    context, products, _ = run_interaction_workflow(
        config_path=args.config,
        run_id=args.run_id or generate_run_id("PR_PPS", "interactions"),
        results_root=root,
        resume=args.resume,
        upstream_run_ids=tuple(args.upstream_run_id),
    )
    print(f"Interaction run: {context.run_dir}")
    print(f"Products: {products}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
