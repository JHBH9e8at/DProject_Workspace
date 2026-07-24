"""Build paired own-state/opposite-state selectivity metrics."""

import pandas as pd


def classify_selectivity(
    own_score,
    opposite_score,
    selectivity_margin,
    selectivity_threshold=2.0,
    strong_score_threshold=-8.0,
):
    """Classify one paired docking result using explicit score thresholds."""
    if selectivity_margin >= selectivity_threshold:
        if own_score <= strong_score_threshold:
            return "selective_own_state"
        return "own_state_preferred_weak"

    if selectivity_margin <= -selectivity_threshold:
        if opposite_score <= strong_score_threshold:
            return "reverse_selective"
        return "reverse_preferred_weak"

    if own_score <= strong_score_threshold and opposite_score <= strong_score_threshold:
        return "nonselective_strong"
    return "inconclusive"


def build_selectivity_table(
    manifest,
    best_scores,
    selectivity_threshold=2.0,
    strong_score_threshold=-8.0,
):
    """Create one paired own/opposite-state score row per completed molecule."""
    manifest_by_id = manifest.set_index("crossdock_id", drop=False)
    rows = []

    for _, best_row in best_scores.iterrows():
        crossdock_id = best_row["crossdock_id"]
        source = manifest_by_id.loc[crossdock_id]
        if isinstance(source, pd.DataFrame):
            raise ValueError(f"Duplicate manifest crossdock_id: {crossdock_id}")

        own_state = source["ligand_state"]
        opposite_state = source["receptor_state"]
        own_score_column = f"{own_state}_r_i_docking_score"
        own_score = float(source[own_score_column])
        opposite_score = float(best_row["r_i_docking_score"])
        margin = opposite_score - own_score

        row = {
            "crossdock_id": crossdock_id,
            "crossdock_arm": source["crossdock_arm"],
            "molecule_id": source["molecule_id"],
            "own_state": own_state,
            "opposite_state": opposite_state,
            "source_step": source.get("step", source.get("source_step", pd.NA)),
            "ligand_prep_mode": best_row.get(
                "ligand_prep_mode", source.get("ligand_prep_mode", pd.NA)
            ),
            "prepared_variant_count": best_row.get("prepared_variant_count", pd.NA),
            "docked_pose_count": best_row.get("docked_pose_count", pd.NA),
            "own_state_score": own_score,
            "opposite_state_score": opposite_score,
            "selectivity_margin": margin,
            "absolute_score_difference": abs(margin),
            "preferred_state": own_state if margin > 0 else (opposite_state if margin < 0 else "tie"),
            "own_state_preferred": margin > 0,
            "selectivity_threshold": selectivity_threshold,
            "strong_score_threshold": strong_score_threshold,
            "selectivity_class": classify_selectivity(
                own_score,
                opposite_score,
                margin,
                selectivity_threshold,
                strong_score_threshold,
            ),
            "best_crossdock_variant": best_row.get("s_lp_Variant", pd.NA),
            "best_crossdock_smiles": best_row.get("smiles", pd.NA),
        }

        own_le_column = f"{own_state}_r_i_glide_ligand_efficiency"
        if own_le_column in source.index and "r_i_glide_ligand_efficiency" in best_row.index:
            own_le = pd.to_numeric(source[own_le_column], errors="coerce")
            opposite_le = pd.to_numeric(
                best_row["r_i_glide_ligand_efficiency"], errors="coerce"
            )
            row["own_state_ligand_efficiency"] = own_le
            row["opposite_state_ligand_efficiency"] = opposite_le
            row["ligand_efficiency_margin"] = opposite_le - own_le

        rows.append(row)

    return pd.DataFrame(rows)


def rank_selectivity(paired_table):
    """Rank candidates by class, margin, and own-state docking strength."""
    priority = {
        "selective_own_state": 0,
        "own_state_preferred_weak": 1,
        "nonselective_strong": 2,
        "inconclusive": 3,
        "reverse_preferred_weak": 4,
        "reverse_selective": 5,
    }
    ranked = paired_table.copy()
    ranked["class_priority"] = ranked["selectivity_class"].map(priority).fillna(99)
    ranked = ranked.sort_values(
        ["class_priority", "selectivity_margin", "own_state_score"],
        ascending=[True, False, True],
    ).reset_index(drop=True)
    ranked.insert(0, "selectivity_rank", range(1, len(ranked) + 1))
    return ranked.drop(columns="class_priority")


def summarize_selectivity(paired_table):
    """Summarize candidate classes separately for PPS- and PR-origin arms."""
    return (
        paired_table.groupby(["crossdock_arm", "own_state", "selectivity_class"], dropna=False)
        .agg(
            molecule_count=("crossdock_id", "size"),
            mean_own_score=("own_state_score", "mean"),
            mean_opposite_score=("opposite_state_score", "mean"),
            mean_selectivity_margin=("selectivity_margin", "mean"),
            median_selectivity_margin=("selectivity_margin", "median"),
        )
        .reset_index()
    )
