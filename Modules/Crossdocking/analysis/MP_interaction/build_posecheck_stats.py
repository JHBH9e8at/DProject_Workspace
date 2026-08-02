import math
from pathlib import Path

import pandas as pd


SOURCE = Path(r"Q:\coding_dir\701_Project\Resultsbin\MP_Interrun\aggregate\pose_quality.csv")
OUTPUT = Path(r"Q:\coding_dir\Workspace1\posecheck_analysis")
OUTPUT.mkdir(parents=True, exist_ok=True)

quality = pd.read_csv(SOURCE, low_memory=False)
required = {
    "population_state",
    "pose_index",
    "molecule_id",
    "docking_score",
    "heavy_atom_count",
    "rotatable_bond_count",
    "clash_count",
    "clashes_per_heavy_atom",
    "strain_energy",
    "strain_per_rotatable_bond",
}
missing = sorted(required - set(quality.columns))
if missing:
    raise ValueError(f"Missing required columns: {missing}")

states = ["PPS", "PR"]
if set(quality["population_state"].dropna().unique()) != set(states):
    raise ValueError("Expected exactly PPS and PR population states")

top_frames = []
selection_rows = []
for state in states:
    frame = quality.loc[quality["population_state"] == state].copy()
    frame = frame.sort_values(
        ["docking_score", "pose_index", "molecule_id"],
        ascending=[True, True, True],
        kind="mergesort",
    )
    selected_n = math.ceil(len(frame) * 0.01)
    selected = frame.head(selected_n).copy()
    selected.insert(0, "top1pct_rank", range(1, selected_n + 1))
    top_frames.append(selected)
    selection_rows.append(
        {
            "population": state,
            "full_population_n": len(frame),
            "selected_n": selected_n,
            "selected_fraction": selected_n / len(frame),
            "best_docking_score": selected["docking_score"].min(),
            "top1pct_boundary_docking_score": selected["docking_score"].max(),
            "next_unselected_docking_score": (
                frame.iloc[selected_n]["docking_score"]
                if selected_n < len(frame)
                else pd.NA
            ),
            "selection_rule": "ceil(N*0.01), ascending docking_score; deterministic tie-break by pose_index then molecule_id",
        }
    )

top1 = pd.concat(top_frames, ignore_index=True)
cohorts = {
    "full_population": quality,
    "docking_top1pct": top1,
}


def summarize(series, statistic):
    values = pd.to_numeric(series, errors="coerce")
    if statistic == "count_nonmissing":
        return int(values.notna().sum())
    if statistic == "count_missing":
        return int(values.isna().sum())
    if statistic == "mean":
        return values.mean()
    if statistic == "median":
        return values.median()
    if statistic == "q05":
        return values.quantile(0.05)
    if statistic == "q25":
        return values.quantile(0.25)
    if statistic == "q75":
        return values.quantile(0.75)
    if statistic == "q95":
        return values.quantile(0.95)
    if statistic == "min":
        return values.min()
    if statistic == "max":
        return values.max()
    raise ValueError(statistic)


metric_specs = [
    ("pose_count", None, "count_rows", "count"),
    ("docking_score_mean", "docking_score", "mean", "Glide score"),
    ("docking_score_median", "docking_score", "median", "Glide score"),
    ("docking_score_q05", "docking_score", "q05", "Glide score"),
    ("docking_score_q95", "docking_score", "q95", "Glide score"),
    ("heavy_atom_count_mean", "heavy_atom_count", "mean", "atoms"),
    ("heavy_atom_count_median", "heavy_atom_count", "median", "atoms"),
    ("rotatable_bond_count_mean", "rotatable_bond_count", "mean", "bonds"),
    ("rotatable_bond_count_median", "rotatable_bond_count", "median", "bonds"),
    ("clash_count_mean", "clash_count", "mean", "count"),
    ("clash_count_median", "clash_count", "median", "count"),
    ("clash_count_q05", "clash_count", "q05", "count"),
    ("clash_count_q25", "clash_count", "q25", "count"),
    ("clash_count_q75", "clash_count", "q75", "count"),
    ("clash_count_q95", "clash_count", "q95", "count"),
    ("clash_count_max", "clash_count", "max", "count"),
    ("clashes_per_heavy_atom_mean", "clashes_per_heavy_atom", "mean", "ratio"),
    ("clashes_per_heavy_atom_median", "clashes_per_heavy_atom", "median", "ratio"),
    ("clashes_per_heavy_atom_q05", "clashes_per_heavy_atom", "q05", "ratio"),
    ("clashes_per_heavy_atom_q25", "clashes_per_heavy_atom", "q25", "ratio"),
    ("clashes_per_heavy_atom_q75", "clashes_per_heavy_atom", "q75", "ratio"),
    ("clashes_per_heavy_atom_q95", "clashes_per_heavy_atom", "q95", "ratio"),
    ("clashes_per_heavy_atom_max", "clashes_per_heavy_atom", "max", "ratio"),
    ("strain_energy_mean", "strain_energy", "mean", "PoseCheck strain energy"),
    ("strain_energy_median", "strain_energy", "median", "PoseCheck strain energy"),
    ("strain_energy_q05", "strain_energy", "q05", "PoseCheck strain energy"),
    ("strain_energy_q25", "strain_energy", "q25", "PoseCheck strain energy"),
    ("strain_energy_q75", "strain_energy", "q75", "PoseCheck strain energy"),
    ("strain_energy_q95", "strain_energy", "q95", "PoseCheck strain energy"),
    ("strain_energy_max", "strain_energy", "max", "PoseCheck strain energy"),
    ("strain_per_rotatable_bond_mean", "strain_per_rotatable_bond", "mean", "ratio"),
    ("strain_per_rotatable_bond_median", "strain_per_rotatable_bond", "median", "ratio"),
    ("strain_per_rotatable_bond_q05", "strain_per_rotatable_bond", "q05", "ratio"),
    ("strain_per_rotatable_bond_q25", "strain_per_rotatable_bond", "q25", "ratio"),
    ("strain_per_rotatable_bond_q75", "strain_per_rotatable_bond", "q75", "ratio"),
    ("strain_per_rotatable_bond_q95", "strain_per_rotatable_bond", "q95", "ratio"),
    ("strain_per_rotatable_bond_max", "strain_per_rotatable_bond", "max", "ratio"),
    ("strain_per_rotatable_bond_missing", "strain_per_rotatable_bond", "count_missing", "count"),
]

summary_rows = []
for scope, cohort in cohorts.items():
    grouped = {state: cohort.loc[cohort["population_state"] == state] for state in states}
    for metric, column, statistic, unit in metric_specs:
        row = {"scope": scope, "metric": metric, "unit": unit}
        for state in states:
            frame = grouped[state]
            row[state] = len(frame) if statistic == "count_rows" else summarize(frame[column], statistic)
        summary_rows.append(row)

relationship_rows = []
quality_metrics = [
    "clash_count",
    "clashes_per_heavy_atom",
    "strain_energy",
    "strain_per_rotatable_bond",
]
for scope, cohort in cohorts.items():
    for method in ["spearman", "pearson"]:
        for metric in quality_metrics:
            row = {
                "scope": scope,
                "correlation_method": method,
                "pose_quality_metric": metric,
            }
            for state in states:
                pair = cohort.loc[
                    cohort["population_state"] == state,
                    ["docking_score", metric],
                ].apply(pd.to_numeric, errors="coerce").dropna()
                if method == "spearman":
                    row[state] = pair["docking_score"].rank().corr(
                        pair[metric].rank(), method="pearson"
                    )
                else:
                    row[state] = pair["docking_score"].corr(
                        pair[metric], method="pearson"
                    )
                row[f"{state}_n"] = len(pair)
            relationship_rows.append(row)

summary_long = pd.DataFrame(summary_rows)
summary = pd.DataFrame(
    [
        {
            "metric": metric,
            "PPS_full": summary_long.loc[
                (summary_long["scope"] == "full_population")
                & (summary_long["metric"] == metric),
                "PPS",
            ].iloc[0],
            "PPS_T1": summary_long.loc[
                (summary_long["scope"] == "docking_top1pct")
                & (summary_long["metric"] == metric),
                "PPS",
            ].iloc[0],
            "PR_full": summary_long.loc[
                (summary_long["scope"] == "full_population")
                & (summary_long["metric"] == metric),
                "PR",
            ].iloc[0],
            "PR_T1": summary_long.loc[
                (summary_long["scope"] == "docking_top1pct")
                & (summary_long["metric"] == metric),
                "PR",
            ].iloc[0],
        }
        for metric, _column, _statistic, _unit in metric_specs
    ]
)
relationships = pd.DataFrame(relationship_rows)
selection = pd.DataFrame(selection_rows)

summary.to_csv(OUTPUT / "pose_quality_summary_full_vs_top1pct.csv", index=False)
relationships.to_csv(OUTPUT / "docking_pose_quality_relationship_full_vs_top1pct.csv", index=False)
selection.to_csv(OUTPUT / "top1pct_selection_summary.csv", index=False)
top1.to_csv(OUTPUT / "pose_quality_top1pct_selected_poses.csv", index=False)

# Reconciliation checks.
assert len(quality) == 40108
assert selection.set_index("population").loc["PPS", "selected_n"] == 205
assert selection.set_index("population").loc["PR", "selected_n"] == 197
assert len(top1) == 402
assert not top1.duplicated(["population_state", "pose_index"]).any()
assert list(summary.columns) == ["metric", "PPS_full", "PPS_T1", "PR_full", "PR_T1"]
assert len(summary) == len(metric_specs)
assert len(relationships) == 2 * 2 * len(quality_metrics)

print(selection.to_string(index=False))
print("\nPose quality summary preview:")
print(summary.head(12).to_string(index=False))
print("\nRelationship table:")
print(relationships.to_string(index=False))
print(f"\nOutput directory: {OUTPUT}")
