"""Run the complete score-analysis block for double-arm docking results."""

import argparse
from pathlib import Path

def run_score_block(
    input_dir,
    output_dir,
    selectivity_threshold=2.0,
    strong_score_threshold=-8.0,
    pps_scores=None,
    pr_scores=None,
    pps_score_col="PPS_r_i_docking_score",
    pr_score_col="PR_r_i_docking_score",
    smiles_col="smiles",
    valid_col="valid",
):
    """Run double-arm QC/selectivity and optional raw-population overlap."""
    from .population_overlap import run_population_overlap
    from .run_analysis import run_analysis

    if bool(pps_scores) != bool(pr_scores):
        raise ValueError("--pps-scores and --pr-scores must be provided together")

    output_dir = Path(output_dir).resolve()
    selectivity_dir = output_dir / "selectivity"
    paired = run_analysis(
        input_dir=input_dir,
        output_dir=selectivity_dir,
        selectivity_threshold=selectivity_threshold,
        strong_score_threshold=strong_score_threshold,
    )

    overlap = None
    if pps_scores and pr_scores:
        overlap = run_population_overlap(
            pps_scores=pps_scores,
            pr_scores=pr_scores,
            output_dir=output_dir / "population_overlap",
            pps_score_col=pps_score_col,
            pr_score_col=pr_score_col,
            smiles_col=smiles_col,
            valid_col=valid_col,
            selectivity_threshold=selectivity_threshold,
        )

    print("Score-analysis block completed")
    print(f"Output directory: {output_dir}")
    return {"selectivity": paired, "population_overlap": overlap}


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="Double-arm result directory")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--selectivity-threshold", type=float, default=2.0)
    parser.add_argument("--strong-score-threshold", type=float, default=-8.0)
    parser.add_argument("--pps-scores", help="Optional raw PPS score CSV")
    parser.add_argument("--pr-scores", help="Optional raw PR score CSV")
    parser.add_argument("--pps-score-col", default="PPS_r_i_docking_score")
    parser.add_argument("--pr-score-col", default="PR_r_i_docking_score")
    parser.add_argument("--smiles-col", default="smiles")
    parser.add_argument("--valid-col", default="valid")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_score_block(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        selectivity_threshold=args.selectivity_threshold,
        strong_score_threshold=args.strong_score_threshold,
        pps_scores=args.pps_scores,
        pr_scores=args.pr_scores,
        pps_score_col=args.pps_score_col,
        pr_score_col=args.pr_score_col,
        smiles_col=args.smiles_col,
        valid_col=args.valid_col,
    )


if __name__ == "__main__":
    main()
