"""Run AHC iteration merging and population cleaning under work_results."""

from __future__ import annotations

import argparse
from pathlib import Path

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import RunContext, create_run_context, generate_run_id
from .population import clean_population, format_cleaning_log, merge_iteration_scores


SCRIPT_ID = "PY-103"


def run_population_preparation(*, iterations_dir: str | Path, run_name: str,
                               run_id: str, results_root: str | Path | None = None,
                               resume: bool = False, upstream_run_ids=()) -> tuple[RunContext, Path, Path, Path]:
    context = create_run_context(("analysis", "ahc_population"), run_id,
                                 results_root=results_root, resume=resume)
    config = {"iterations_dir": str(Path(iterations_dir).resolve(strict=False)), "run_name": run_name}
    try:
        merged, included, excluded = merge_iteration_scores(iterations_dir)
        merged_path, cleaned_path, log_path = allocate_versioned_group(
            context.tables_dir,
            (f"{run_name}_merged_scores.csv", f"{run_name}_cleaned_population.csv", f"{run_name}_cleaning_log.txt"),
        )
        cleaned, steps = clean_population(merged, run_name=run_name)
        merged.to_csv(merged_path, index=False)
        cleaned.to_csv(cleaned_path, index=False)
        log_path.write_text(format_cleaning_log(
            input_label=str(merged_path), output_label=str(cleaned_path),
            initial_count=len(merged), final_count=len(cleaned), steps=steps,
        ), encoding="utf-8")
        config["excluded_final_iteration"] = str(excluded)
        update_manifest(context, script_id=SCRIPT_ID, status="completed", config=config,
                        inputs=included, outputs=(merged_path, cleaned_path, log_path), upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config=config, upstream_run_ids=upstream_run_ids)
        raise
    return context, merged_path, cleaned_path, log_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Merge and clean one AHC population")
    parser.add_argument("--iterations-dir", required=True, type=Path)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--run-id")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--results-root", type=Path); parser.add_argument("--upstream-run-id", action="append", default=[])
    args = parser.parse_args(argv)
    run_population_preparation(iterations_dir=args.iterations_dir, run_name=args.run_name,
                               run_id=args.run_id or generate_run_id(args.run_name, "ahc_population"),
                               results_root=validate_results_root(args.results_root) if args.results_root else None,
                               resume=args.resume, upstream_run_ids=tuple(args.upstream_run_id))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
