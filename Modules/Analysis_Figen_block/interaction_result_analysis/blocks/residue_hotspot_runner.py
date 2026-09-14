"""Write residue-hotspot analysis products below work_results."""

from __future__ import annotations

import argparse
from pathlib import Path

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from .config_helpers import parse_int_assignments
from .residue_hotspot import summarize_interactions


SCRIPT_ID = "PY-123"


def run_residue_hotspot(*, interactions_csv: str | Path, denominators: dict[str, int],
                        run_id: str, results_root: str | Path | None = None,
                        resume: bool = False,
                        upstream_run_ids: tuple[str, ...] = ()):
    context = create_run_context(("analysis", "residue_hotspot"), run_id,
                                 results_root=results_root, resume=resume)
    config = {"pose_denominators": denominators}
    try:
        residue, typed, dominant = summarize_interactions(interactions_csv, denominators)
        paths = allocate_versioned_group(context.tables_dir,
            ("residue_prevalence.csv", "interaction_type_prevalence.csv", "dominant_interaction_type.csv"))
        for table, path in zip((residue, typed, dominant), paths): table.to_csv(path, index=False)
        update_manifest(context, script_id=SCRIPT_ID, status="completed", config=config,
                        inputs=(Path(interactions_csv),), outputs=paths,
                        upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config=config,
                        upstream_run_ids=upstream_run_ids); raise
    return context, (residue, typed, dominant), paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build residue-hotspot tables.")
    parser.add_argument("--interactions", required=True, type=Path)
    parser.add_argument("--denominator", action="append", required=True)
    parser.add_argument("--run-id"); parser.add_argument("--results-root", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--upstream-run-id", action="append", default=[])
    args = parser.parse_args(argv)
    root = validate_results_root(args.results_root) if args.results_root else None
    context, _, paths = run_residue_hotspot(
        interactions_csv=args.interactions,
        denominators=parse_int_assignments(args.denominator, label="denominator"),
        run_id=args.run_id or generate_run_id("PR_PPS", "residue_hotspot"),
        results_root=root, resume=args.resume,
        upstream_run_ids=tuple(args.upstream_run_id),
    )
    print(f"Hotspot run: {context.run_dir}")
    for path in paths: print(f"Table: {path}")
    return 0


if __name__ == "__main__": raise SystemExit(main())
