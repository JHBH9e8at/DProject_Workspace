"""Run consolidated scaffold-diversity analysis under work_results."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from .scaffold_diversity import analyze_source, load_population


SCRIPT_ID = "PY-110"
TABLE_NAMES = ("fp_by_step", "scaffold_by_step", "aggregate", "top_scaffolds_last_window")


def run_scaffold_diversity(*, populations: dict[str, str | Path], membership_csv: str | Path,
                           run_id: str, results_root: str | Path | None = None,
                           resume: bool = False, last_n: int = 50, fp_radius: int = 2,
                           fp_size: int = 16384, diversity_sample_size: int = 3000,
                           novelty_reference_size: int = 5000, random_state: int = 42,
                           upstream_run_ids=()):
    context = create_run_context(("analysis", "scaffold_diversity"), run_id,
                                 results_root=results_root, resume=resume)
    config = {"last_n": last_n, "fp_radius": fp_radius, "fp_size": fp_size,
              "diversity_sample_size": diversity_sample_size,
              "novelty_reference_size": novelty_reference_size, "random_state": random_state,
              "sources": sorted(populations)}
    inputs = [Path(membership_csv), *(Path(path) for path in populations.values())]
    try:
        combined: dict[str, list[pd.DataFrame]] = {}
        for source, path in populations.items():
            frame = load_population(path, source, membership_csv)
            results = analyze_source(frame, source, last_n=last_n, fp_radius=fp_radius,
                fp_size=fp_size, diversity_sample_size=diversity_sample_size,
                novelty_reference_size=novelty_reference_size, random_state=random_state)
            for name, table in results.items(): combined.setdefault(name, []).append(table)
        tables = {name: pd.concat(parts, ignore_index=True) for name, parts in combined.items()}
        paths = allocate_versioned_group(context.tables_dir, tuple(f"{name}.csv" for name in TABLE_NAMES))
        for name, path in zip(TABLE_NAMES, paths): tables[name].to_csv(path, index=False)
        update_manifest(context, script_id=SCRIPT_ID, status="completed", config=config,
                        inputs=inputs, outputs=paths, upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config=config, upstream_run_ids=upstream_run_ids); raise
    return context, tables, dict(zip(TABLE_NAMES, paths))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Analyze AHC scaffold and fingerprint diversity")
    parser.add_argument("--pr-input", required=True, type=Path); parser.add_argument("--pps-input", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path); parser.add_argument("--run-id")
    parser.add_argument("--last-n", type=int, default=50); parser.add_argument("--resume", action="store_true")
    parser.add_argument("--results-root", type=Path); parser.add_argument("--upstream-run-id", action="append", default=[])
    args = parser.parse_args(argv)
    run_scaffold_diversity(populations={"PR": args.pr_input, "PPS": args.pps_input},
        membership_csv=args.membership, run_id=args.run_id or generate_run_id("PR_PPS", "scaffold_diversity"),
        last_n=args.last_n, results_root=validate_results_root(args.results_root) if args.results_root else None,
        resume=args.resume, upstream_run_ids=tuple(args.upstream_run_id))
    return 0


if __name__ == "__main__": raise SystemExit(main())
