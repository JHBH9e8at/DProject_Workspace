"""Command-line entry point for double-arm state-selectivity analysis."""

import argparse
from pathlib import Path

try:
    from .selectivity import build_selectivity_table, rank_selectivity, summarize_selectivity
    from .validation import load_double_arm_results, qc_has_errors, validate_double_arm_results
except ImportError:
    from selectivity import build_selectivity_table, rank_selectivity, summarize_selectivity
    from validation import load_double_arm_results, qc_has_errors, validate_double_arm_results


def run_analysis(
    input_dir,
    output_dir,
    selectivity_threshold=2.0,
    strong_score_threshold=-8.0,
):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest, best_scores = load_double_arm_results(input_dir)
    qc_report = validate_double_arm_results(manifest, best_scores)
    qc_file = output_dir / "qc_report.csv"
    qc_report.to_csv(qc_file, index=False)

    if qc_has_errors(qc_report):
        raise ValueError(f"Double-arm QC failed. See {qc_file}")

    paired = build_selectivity_table(
        manifest=manifest,
        best_scores=best_scores,
        selectivity_threshold=selectivity_threshold,
        strong_score_threshold=strong_score_threshold,
    )
    ranking = rank_selectivity(paired)
    summary = summarize_selectivity(paired)

    paired_file = output_dir / "paired_state_scores.csv"
    ranking_file = output_dir / "selectivity_ranking.csv"
    summary_file = output_dir / "selectivity_summary.csv"
    paired.to_csv(paired_file, index=False)
    ranking.to_csv(ranking_file, index=False)
    summary.to_csv(summary_file, index=False)

    print("Double-arm selectivity analysis completed")
    print(f"Analyzed molecules: {len(paired)}")
    print(f"QC report: {qc_file}")
    print(f"Paired scores: {paired_file}")
    print(f"Ranking: {ranking_file}")
    print(f"Summary: {summary_file}")
    return paired


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze paired own-state and opposite-state double-arm docking scores"
    )
    parser.add_argument(
        "--input-dir", required=True, help="double_arm_docking.py output directory"
    )
    parser.add_argument("--output-dir", required=True, help="Analysis output directory")
    parser.add_argument(
        "--selectivity-threshold",
        type=float,
        default=2.0,
        help="Minimum absolute score difference for a selective classification",
    )
    parser.add_argument(
        "--strong-score-threshold",
        type=float,
        default=-8.0,
        help="Docking score cutoff for a strong-binding classification",
    )
    args = parser.parse_args()

    run_analysis(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        selectivity_threshold=args.selectivity_threshold,
        strong_score_threshold=args.strong_score_threshold,
    )
