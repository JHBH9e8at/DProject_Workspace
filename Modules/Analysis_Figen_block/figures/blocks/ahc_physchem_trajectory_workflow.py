"""AHC physicochemical/trajectory analysis and Figure workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from ...ahc_analysis.blocks.physchem_trajectory import build_physchem_trajectory_tables
from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from .ahc_physchem_trajectory import render_figures
from ...ahc_analysis.blocks.physchem_trajectory import SUPPORTED_DESCRIPTORS


SCRIPT_ID = "PY-107"
DEFAULT_REF_CSV = (
    Path(__file__).resolve().parents[4]
    / "AHC_Related" / "ref_structures" / "ref_desc.csv"
)


def run_workflow(*, input_csv: str | Path, job_name: str, threshold: float,
                 run_id: str, results_root: str | Path | None = None,
                 resume: bool = False, ref_csv: str | Path | None = None,
                 ref_keys: list[str] | None = None,
                 descriptors: list[str] | None = None,
                 mw_upper_bound: float | None = None,
                 upstream_run_ids: tuple[str, ...] = ()):
    selected_refs = ["OM"] if ref_keys is None else ref_keys
    selected_ref_csv = DEFAULT_REF_CSV if selected_refs and ref_csv is None else ref_csv
    analysis = create_run_context(("analysis", "ahc_physchem_trajectory"), run_id,
                                  results_root=results_root, resume=resume)
    figures = create_run_context(("Figs", "ahc_physchem_trajectory", job_name), run_id,
                                 results_root=results_root, resume=resume)
    config = {"job_name": job_name, "threshold": threshold,
              "input_csv": str(Path(input_csv).resolve(strict=False)),
              "ref_keys": selected_refs, "descriptors": descriptors or [],
              "ref_csv": str(selected_ref_csv) if selected_ref_csv is not None else None,
              "mw_upper_bound": mw_upper_bound}
    try:
        rows, trajectory, summary = build_physchem_trajectory_tables(
            input_csv, job_name=job_name, threshold=threshold, descriptors=descriptors)
        row_path, trajectory_path, summary_path = allocate_versioned_group(
            analysis.tables_dir, (f"{job_name}_physchem_rows.csv",
                                  f"{job_name}_docking_trajectory.csv",
                                  f"{job_name}_physchem_summary.csv"))
        rows.to_csv(row_path, index=False); trajectory.to_csv(trajectory_path, index=False)
        summary.to_csv(summary_path, index=False)
        update_manifest(analysis, script_id=SCRIPT_ID, status="completed", config=config,
                        inputs=(Path(input_csv),), outputs=(row_path, trajectory_path, summary_path),
                        upstream_run_ids=upstream_run_ids)

        names = (
            f"{job_name}_physchem_distribution.png", f"{job_name}_physchem_distribution.svg",
            f"{job_name}_physchem_docking_scatter.png", f"{job_name}_physchem_docking_scatter.svg",
            f"{job_name}_docking_score_trajectory.png", f"{job_name}_docking_score_trajectory.svg",
        )
        paths = allocate_versioned_group(figures.run_dir, names)
        keys = ("distribution_png", "distribution_svg", "scatter_png", "scatter_svg",
                "trajectory_png", "trajectory_svg")
        output_paths = dict(zip(keys, paths))
        render_figures(rows_csv=row_path, trajectory_csv=trajectory_path,
                       output_paths=output_paths, job_name=job_name, threshold=threshold,
                       ref_csv=selected_ref_csv, ref_keys=selected_refs,
                       descriptors=descriptors, mw_upper_bound=mw_upper_bound)

        figure_inputs = [row_path, trajectory_path]
        if selected_ref_csv is not None and selected_refs:
            figure_inputs.append(Path(selected_ref_csv))
        figure_config = {**config, "analysis_run": str(analysis.run_dir)}
        update_manifest(figures, script_id="PY-106", status="completed", config=figure_config,
                        inputs=figure_inputs, outputs=paths,
                        upstream_run_ids=(*upstream_run_ids, analysis.run_id))
    except Exception:
        update_manifest(analysis, script_id=SCRIPT_ID, status="failed", config=config,
                        upstream_run_ids=upstream_run_ids)
        update_manifest(figures, script_id="PY-106", status="failed", config=config,
                        upstream_run_ids=(*upstream_run_ids, analysis.run_id))
        raise
    return analysis, figures, (row_path, trajectory_path, summary_path), output_paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate AHC physicochemical and trajectory products")
    parser.add_argument("--input", required=True, type=Path); parser.add_argument("--job-name", required=True)
    parser.add_argument("--threshold", required=True, type=float); parser.add_argument("--run-id")
    parser.add_argument("--resume", action="store_true"); parser.add_argument("--ref-csv", type=Path)
    parser.add_argument("--ref", nargs="*", choices=("OM", "MAV", "AFI"), default=None,
                        help="Reference molecules (default: OM); pass --ref alone to disable")
    parser.add_argument("--property", dest="descriptors", action="append",
                        choices=SUPPORTED_DESCRIPTORS,
                        help="Physicochemical descriptor to plot; repeat to select multiple")
    parser.add_argument("--mw-upper-bound", type=float,
                        help="Shared molecular-weight x-axis upper bound for comparisons")
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--upstream-run-id", action="append", default=[])
    args = parser.parse_args(argv)
    results_root = (
        validate_results_root(args.results_root)
        if args.results_root is not None
        else None
    )
    run_workflow(input_csv=args.input, job_name=args.job_name, threshold=args.threshold,
                 run_id=args.run_id or generate_run_id(args.job_name, "ahc_physchem"),
                 results_root=results_root, resume=args.resume,
                 ref_csv=args.ref_csv, ref_keys=args.ref, descriptors=args.descriptors,
                 mw_upper_bound=args.mw_upper_bound,
                 upstream_run_ids=tuple(args.upstream_run_id))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
