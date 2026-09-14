"""Isolated adapter for the existing cross-docking score-analysis chain."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..common.output_naming import next_available_path
from ..common.paths import validate_results_root
from ..common.provenance import update_manifest
from ..common.run_context import create_run_context, generate_run_id
from .population_overlap import run_population_overlap
from .run_analysis import run_analysis


SCRIPT_ID = "PY-116"


def _resolve_double_arm_dir(double_arm_dir, calculation_manifest):
    if (double_arm_dir is None) == (calculation_manifest is None):
        raise ValueError("Provide exactly one of double_arm_dir or calculation_manifest")
    if calculation_manifest is None:
        return Path(double_arm_dir), None, None
    path = Path(calculation_manifest).resolve(strict=True)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or payload.get("status") != "completed":
        raise ValueError(f"Invalid or incomplete calculation manifest: {path}")
    directory = Path(payload["double_arm_dir"]).resolve(strict=True)
    declared = {record["path"] for record in payload.get("files", [])}
    if not any(value.endswith("double_arm_manifest.csv") for value in declared) or not any(
            value.endswith("double_arm_best_docking_scores.csv") for value in declared):
        raise ValueError("Calculation manifest lacks required double-arm products")
    upstream_run_id = payload.get("run_id")
    if not upstream_run_id:
        output_root = payload.get("output_root")
        upstream_run_id = Path(output_root).name if output_root else None
    return directory, path, upstream_run_id


def _merge_upstream_run_ids(*groups):
    """Preserve upstream order while removing empty and duplicate IDs."""

    merged = []
    for group in groups:
        for value in group:
            if value and value not in merged:
                merged.append(value)
    return tuple(merged)


def run_score_analysis(*, double_arm_dir: str | Path | None = None,
                       calculation_manifest: str | Path | None = None, run_id: str,
                       results_root: str | Path | None = None, resume: bool = False,
                       pps_scores: str | Path | None = None, pr_scores: str | Path | None = None,
                       selectivity_threshold: float = 2.0,
                       strong_score_threshold: float = -8.0,
                       upstream_run_ids: tuple[str, ...] = ()):
    double_arm_dir, manifest_path, calculation_run_id = _resolve_double_arm_dir(
        double_arm_dir, calculation_manifest
    )
    resolved_upstream_ids = _merge_upstream_run_ids(
        upstream_run_ids, (calculation_run_id,)
    )
    if (pps_scores is None) != (pr_scores is None):
        raise ValueError("pps_scores and pr_scores must be provided together")
    context = create_run_context(("analysis", "crossdocking_scores"), run_id,
                                 results_root=results_root, resume=resume)
    product_dir = next_available_path(context.run_dir, "products")
    product_dir.mkdir()
    selectivity_dir = product_dir / "selectivity"
    config = {"selectivity_threshold": selectivity_threshold,
              "strong_score_threshold": strong_score_threshold,
              "population_overlap": pps_scores is not None}
    inputs = [Path(double_arm_dir) / "double_arm_manifest.csv",
              Path(double_arm_dir) / "double_arm_best_docking_scores.csv"]
    if manifest_path is not None:
        inputs.insert(0, manifest_path)
    try:
        paired = run_analysis(double_arm_dir, selectivity_dir,
                              selectivity_threshold, strong_score_threshold)
        overlap = None
        if pps_scores is not None:
            overlap = run_population_overlap(pps_scores, pr_scores,
                product_dir / "population_overlap", selectivity_threshold=selectivity_threshold)
            inputs.extend([Path(pps_scores), Path(pr_scores)])
        outputs = sorted(path for path in product_dir.rglob("*.csv"))
        update_manifest(context, script_id=SCRIPT_ID, status="completed", config=config,
                        inputs=inputs, outputs=outputs,
                        upstream_run_ids=resolved_upstream_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config=config,
                        upstream_run_ids=resolved_upstream_ids); raise
    return context, product_dir, paired, overlap


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run cross-docking score QC, selectivity, and optional overlap."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--calculation-manifest", type=Path)
    source.add_argument("--double-arm-dir", type=Path)
    parser.add_argument("--pps-scores", type=Path)
    parser.add_argument("--pr-scores", type=Path)
    parser.add_argument("--selectivity-threshold", type=float, default=2.0)
    parser.add_argument("--strong-score-threshold", type=float, default=-8.0)
    parser.add_argument("--run-id")
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--upstream-run-id", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    results_root = (
        validate_results_root(args.results_root)
        if args.results_root is not None
        else None
    )
    context, products, paired, overlap = run_score_analysis(
        double_arm_dir=args.double_arm_dir,
        calculation_manifest=args.calculation_manifest,
        run_id=args.run_id or generate_run_id("PR_PPS", "crossdocking_scores"),
        results_root=results_root,
        resume=args.resume,
        pps_scores=args.pps_scores,
        pr_scores=args.pr_scores,
        selectivity_threshold=args.selectivity_threshold,
        strong_score_threshold=args.strong_score_threshold,
        upstream_run_ids=tuple(args.upstream_run_id),
    )
    print(f"Analysis run: {context.run_dir}")
    print(f"Products: {products}")
    print(f"Paired molecules: {len(paired)}")
    if overlap is not None:
        print(f"Overlapping molecules: {len(overlap)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
