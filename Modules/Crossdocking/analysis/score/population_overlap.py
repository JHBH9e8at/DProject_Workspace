"""Find canonical molecules shared by raw PPS and PR AHC populations."""

import argparse
from pathlib import Path

import pandas as pd

try:
    from .identity import aggregate_population, prepare_population_instances
except ImportError:
    from identity import aggregate_population, prepare_population_instances


def score_preference(margin, threshold):
    if margin >= threshold:
        return "PPS"
    if margin <= -threshold:
        return "PR"
    return "inconclusive"


def classify_overlap(row, threshold):
    best_margin = row["best_margin_PR_minus_PPS"]
    median_margin = row["median_margin_PR_minus_PPS"]
    best_preference = score_preference(best_margin, threshold)
    median_preference = score_preference(median_margin, threshold)

    if best_preference == median_preference and best_preference == "PPS":
        return "robust_PPS_preference"
    if best_preference == median_preference and best_preference == "PR":
        return "robust_PR_preference"
    if best_preference == median_preference == "inconclusive":
        return "inconclusive"
    if best_margin > 0 and median_margin > 0:
        return "weak_PPS_preference"
    if best_margin < 0 and median_margin < 0:
        return "weak_PR_preference"
    return "best_median_disagreement"


def build_overlap_table(pps_aggregated, pr_aggregated, selectivity_threshold=2.0):
    """Inner-join canonical molecules and calculate paired state differences."""
    overlap = pps_aggregated.merge(
        pr_aggregated,
        on=["canonical_isomeric_smiles", "canonical_formal_charge"],
        how="inner",
        validate="one_to_one",
    )

    overlap["best_margin_PR_minus_PPS"] = (
        overlap["PR_best_score"] - overlap["PPS_best_score"]
    )
    overlap["median_margin_PR_minus_PPS"] = (
        overlap["PR_median_score"] - overlap["PPS_median_score"]
    )
    overlap["mean_margin_PR_minus_PPS"] = (
        overlap["PR_mean_score"] - overlap["PPS_mean_score"]
    )
    overlap["best_preferred_state"] = overlap["best_margin_PR_minus_PPS"].map(
        lambda margin: score_preference(margin, selectivity_threshold)
    )
    overlap["median_preferred_state"] = overlap["median_margin_PR_minus_PPS"].map(
        lambda margin: score_preference(margin, selectivity_threshold)
    )
    overlap["best_median_sign_agreement"] = (
        overlap["best_margin_PR_minus_PPS"] * overlap["median_margin_PR_minus_PPS"]
    ) > 0
    overlap["selectivity_threshold"] = selectivity_threshold
    if overlap.empty:
        overlap["overlap_class"] = pd.Series(dtype="object")
    else:
        overlap["overlap_class"] = overlap.apply(
            classify_overlap, axis=1, threshold=selectivity_threshold
        )

    class_priority = {
        "robust_PPS_preference": 0,
        "weak_PPS_preference": 1,
        "inconclusive": 2,
        "best_median_disagreement": 3,
        "weak_PR_preference": 4,
        "robust_PR_preference": 5,
    }
    overlap["_class_priority"] = overlap["overlap_class"].map(class_priority).fillna(99)
    overlap = overlap.sort_values(
        ["_class_priority", "best_margin_PR_minus_PPS"],
        ascending=[True, False],
    ).drop(columns="_class_priority")
    overlap.insert(0, "overlap_rank", range(1, len(overlap) + 1))
    return overlap.reset_index(drop=True)


def build_summary(
    pps_total_rows,
    pr_total_rows,
    pps_valid,
    pr_valid,
    pps_aggregated,
    pr_aggregated,
    overlap,
    pps_excluded,
    pr_excluded,
):
    """Create an auditable key/value summary of population overlap coverage."""
    pps_unique = len(pps_aggregated)
    pr_unique = len(pr_aggregated)
    overlap_count = len(overlap)
    union_count = pps_unique + pr_unique - overlap_count

    records = [
        ("PPS_input_instances", pps_total_rows),
        ("PR_input_instances", pr_total_rows),
        ("PPS_valid_instances", len(pps_valid)),
        ("PR_valid_instances", len(pr_valid)),
        ("PPS_excluded_instances", len(pps_excluded)),
        ("PR_excluded_instances", len(pr_excluded)),
        ("PPS_unique_canonical_molecules", pps_unique),
        ("PR_unique_canonical_molecules", pr_unique),
        ("overlapping_canonical_molecules", overlap_count),
        ("canonical_molecule_union", union_count),
        (
            "overlap_fraction_of_PPS",
            overlap_count / pps_unique if pps_unique else 0.0,
        ),
        (
            "overlap_fraction_of_PR",
            overlap_count / pr_unique if pr_unique else 0.0,
        ),
        (
            "canonical_jaccard_similarity",
            overlap_count / union_count if union_count else 0.0,
        ),
    ]

    if not overlap.empty:
        records.extend(
            [
                ("median_best_margin_PR_minus_PPS", overlap["best_margin_PR_minus_PPS"].median()),
                ("mean_best_margin_PR_minus_PPS", overlap["best_margin_PR_minus_PPS"].mean()),
                (
                    "robust_PPS_preference_count",
                    int((overlap["overlap_class"] == "robust_PPS_preference").sum()),
                ),
                (
                    "robust_PR_preference_count",
                    int((overlap["overlap_class"] == "robust_PR_preference").sum()),
                ),
                (
                    "best_median_disagreement_count",
                    int((overlap["overlap_class"] == "best_median_disagreement").sum()),
                ),
            ]
        )

    return pd.DataFrame(records, columns=["metric", "value"])


def run_population_overlap(
    pps_scores,
    pr_scores,
    output_dir,
    pps_score_col="PPS_r_i_docking_score",
    pr_score_col="PR_r_i_docking_score",
    smiles_col="smiles",
    valid_col="valid",
    selectivity_threshold=2.0,
):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    pps_valid, pps_audit, pps_excluded, pps_total = prepare_population_instances(
        pps_scores, "PPS", pps_score_col, smiles_col, valid_col
    )
    pr_valid, pr_audit, pr_excluded, pr_total = prepare_population_instances(
        pr_scores, "PR", pr_score_col, smiles_col, valid_col
    )

    pps_aggregated = aggregate_population(pps_valid, "PPS", pps_score_col)
    pr_aggregated = aggregate_population(pr_valid, "PR", pr_score_col)
    overlap = build_overlap_table(
        pps_aggregated, pr_aggregated, selectivity_threshold
    )

    pps_duplicates = pps_aggregated[
        pps_aggregated["PPS_occurrence_count"] > 1
    ].copy()
    pr_duplicates = pr_aggregated[
        pr_aggregated["PR_occurrence_count"] > 1
    ].copy()
    excluded = pd.concat([pps_excluded, pr_excluded], ignore_index=True)
    summary = build_summary(
        pps_total,
        pr_total,
        pps_audit,
        pr_audit,
        pps_aggregated,
        pr_aggregated,
        overlap,
        pps_excluded,
        pr_excluded,
    )

    outputs = {
        "population_overlap_summary.csv": summary,
        "overlapping_molecules.csv": overlap,
        "PPS_population_aggregated.csv": pps_aggregated,
        "PR_population_aggregated.csv": pr_aggregated,
        "PPS_internal_duplicates.csv": pps_duplicates,
        "PR_internal_duplicates.csv": pr_duplicates,
        "PPS_canonicalized_instances.csv": pps_audit,
        "PR_canonicalized_instances.csv": pr_audit,
        "excluded_instances.csv": excluded,
    }
    for filename, table in outputs.items():
        table.to_csv(output_dir / filename, index=False)

    print("Population-overlap analysis completed")
    print(f"PPS valid instances: {len(pps_audit)}")
    print(f"PR valid instances: {len(pr_audit)}")
    print(f"PPS unique canonical molecules: {len(pps_aggregated)}")
    print(f"PR unique canonical molecules: {len(pr_aggregated)}")
    print(f"Overlapping canonical molecules: {len(overlap)}")
    print(f"Output directory: {output_dir}")
    return overlap


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare canonical molecules shared by raw PPS and PR AHC populations"
    )
    parser.add_argument("--pps-scores", required=True, help="<Raw> PPS scores.csv")
    parser.add_argument("--pr-scores", required=True, help="<Raw> PR scores.csv")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--pps-score-col", default="PPS_r_i_docking_score", help="colNameforDscore")
    parser.add_argument("--pr-score-col", default="PR_r_i_docking_score", help="colNameforDscore")
    parser.add_argument("--smiles-col", default="smiles")
    parser.add_argument("--valid-col", default="valid")
    parser.add_argument(
        "--selectivity-threshold",
        type=float,
        default=2.0,
        help="Minimum absolute score difference for robust state preference",
    )
    args = parser.parse_args()

    run_population_overlap(
        pps_scores=args.pps_scores,
        pr_scores=args.pr_scores,
        output_dir=args.output_dir,
        pps_score_col=args.pps_score_col,
        pr_score_col=args.pr_score_col,
        smiles_col=args.smiles_col,
        valid_col=args.valid_col,
        selectivity_threshold=args.selectivity_threshold,
    )


# exmp

# python -m analysis.score.population_overlap \
#   --pps-scores /path/to/PPS/scores.csv \
#   --pr-scores /path/to/PR/scores.csv \
#   --output-dir /path/to/population_overlap \
#   --selectivity-threshold 2.0
