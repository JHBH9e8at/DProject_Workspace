"""Create redocking analysis tables and figures in isolated result runs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import RunContext, create_run_context, generate_run_id
from .redocking_validation import build_redocking_tables


WORKFLOW_SCRIPT_ID = "PY-099"
ANALYSIS_SCRIPT_ID = "PY-095"
FIGURE_SCRIPT_ID = "PY-097"
FIGURE_STEMS = (
    "redocking_docking_score_comparison",
    "redocking_rmsd_comparison",
    "redocking_PR_mirrored",
    "redocking_PPS_mirrored",
)


def _figure_path_map(paths: tuple[Path, ...]) -> dict[str, tuple[Path, Path]]:
    if len(paths) != len(FIGURE_STEMS) * 2:
        raise ValueError("Expected one PNG/SVG pair per redocking figure")
    return {
        stem: (paths[index * 2], paths[index * 2 + 1])
        for index, stem in enumerate(FIGURE_STEMS)
    }


def run_redocking_validation(
    *,
    pps_docking_csv: str | Path,
    pps_rmsd_csv: str | Path,
    pr_docking_csv: str | Path,
    pr_rmsd_csv: str | Path,
    run_id: str,
    results_root: str | Path | None = None,
    resume: bool = False,
    score_ylim: tuple[float, float] = (-12, -1),
    rmsd_ylim: tuple[float, float] = (0, 2.5),
    upstream_run_ids: tuple[str, ...] = (),
) -> tuple[RunContext, RunContext, list[Path], list[Path]]:
    """Write analysis tables, render figures, and record both manifests."""

    inputs = tuple(
        Path(path).expanduser().resolve(strict=False)
        for path in (
            pps_docking_csv,
            pps_rmsd_csv,
            pr_docking_csv,
            pr_rmsd_csv,
        )
    )
    config = {
        "workflow_script_id": WORKFLOW_SCRIPT_ID,
        "pps_docking_csv": str(inputs[0]),
        "pps_rmsd_csv": str(inputs[1]),
        "pr_docking_csv": str(inputs[2]),
        "pr_rmsd_csv": str(inputs[3]),
        "score_ylim": list(score_ylim),
        "rmsd_ylim": list(rmsd_ylim),
    }
    analysis_context = create_run_context(
        ("analysis", "redocking"),
        run_id,
        results_root=results_root,
        resume=resume,
    )
    try:
        poses, summary = build_redocking_tables(
            pps_docking_csv=inputs[0],
            pps_rmsd_csv=inputs[1],
            pr_docking_csv=inputs[2],
            pr_rmsd_csv=inputs[3],
        )
        pose_path, summary_path = allocate_versioned_group(
            analysis_context.tables_dir,
            ("redocking_poses.csv", "redocking_summary.csv"),
        )
        poses.to_csv(pose_path, index=False)
        summary.to_csv(summary_path, index=False)
        update_manifest(
            analysis_context,
            script_id=ANALYSIS_SCRIPT_ID,
            status="completed",
            config=config,
            inputs=inputs,
            outputs=(pose_path, summary_path),
            upstream_run_ids=upstream_run_ids,
        )
    except Exception:
        update_manifest(
            analysis_context,
            script_id=ANALYSIS_SCRIPT_ID,
            status="failed",
            config=config,
            upstream_run_ids=upstream_run_ids,
        )
        raise

    figure_context = create_run_context(
        ("Figs", "redocking_validation"),
        run_id,
        results_root=analysis_context.results_root,
        resume=resume,
    )
    try:
        matplotlib_config = figure_context.intermediate_dir / ".matplotlib"
        matplotlib_config.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_config))
        from ...figures.blocks.redocking_validation import render_figure_set

        requested_names = tuple(
            name
            for stem in FIGURE_STEMS
            for name in (f"{stem}.png", f"{stem}.svg")
        )
        figure_paths = allocate_versioned_group(
            figure_context.run_dir,
            requested_names,
        )
        pose_product = pd.read_csv(pose_path)
        written = render_figure_set(
            pose_product,
            _figure_path_map(figure_paths),
            score_ylim=score_ylim,
            rmsd_ylim=rmsd_ylim,
        )
        update_manifest(
            figure_context,
            script_id=FIGURE_SCRIPT_ID,
            status="completed",
            config=config,
            inputs=(pose_path,),
            outputs=written,
            upstream_run_ids=(*upstream_run_ids, analysis_context.run_id),
        )
    except Exception:
        update_manifest(
            figure_context,
            script_id=FIGURE_SCRIPT_ID,
            status="failed",
            config=config,
            upstream_run_ids=(*upstream_run_ids, analysis_context.run_id),
        )
        raise
    return (
        analysis_context,
        figure_context,
        [pose_path, summary_path],
        written,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate redocking analysis tables and validation figures."
    )
    parser.add_argument("--pps-docking", required=True, type=Path)
    parser.add_argument("--pps-rmsd", required=True, type=Path)
    parser.add_argument("--pr-docking", required=True, type=Path)
    parser.add_argument("--pr-rmsd", required=True, type=Path)
    parser.add_argument("--run-id")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--upstream-run-id", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    results_root = (
        validate_results_root(args.results_root)
        if args.results_root is not None
        else None
    )
    run_id = args.run_id or generate_run_id("PR_PPS", "redocking")
    analysis, figures, tables, figure_paths = run_redocking_validation(
        pps_docking_csv=args.pps_docking,
        pps_rmsd_csv=args.pps_rmsd,
        pr_docking_csv=args.pr_docking,
        pr_rmsd_csv=args.pr_rmsd,
        run_id=run_id,
        results_root=results_root,
        resume=args.resume,
        upstream_run_ids=tuple(args.upstream_run_id),
    )
    print(f"Analysis run: {analysis.run_dir}")
    for path in tables:
        print(f"Table: {path}")
    print(f"Figure run: {figures.run_dir}")
    for path in figure_paths:
        print(f"Figure: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
