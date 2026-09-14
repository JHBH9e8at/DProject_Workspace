"""Stable data products for AHC physicochemical and docking-trajectory figures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


SCRIPT_ID = "PY-105"
DESCRIPTORS = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_TPSA", "desc_NumHAcceptors", "desc_NumHDonors",
]
OPTIONAL_DESCRIPTORS = ["desc_NumRotatableBonds"]
SUPPORTED_DESCRIPTORS = [*DESCRIPTORS, *OPTIONAL_DESCRIPTORS]


def build_physchem_trajectory_tables(
    input_csv: str | Path, *, job_name: str, threshold: float,
    descriptors: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source = Path(input_csv)
    frame = pd.read_csv(source)
    score_col = f"{job_name}_r_i_docking_score"
    selected = descriptors or DESCRIPTORS
    unknown = [column for column in selected if column not in SUPPORTED_DESCRIPTORS]
    if unknown:
        raise ValueError(f"Unsupported physicochemical descriptors: {unknown}")
    required = ["step", score_col, *selected]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{source}: missing required columns: {missing}")

    rows = frame[required].copy()
    for column in required:
        rows[column] = pd.to_numeric(rows[column], errors="coerce")
    if rows["step"].isna().any() or rows[score_col].isna().any():
        raise ValueError(f"{source}: step or docking score contains non-numeric values")

    max_step = rows["step"].max()
    first_last50 = max_step - 50
    rows["last50"] = rows["step"] >= first_last50
    rows["docking_class"] = rows[score_col].le(threshold).map({True: "good", False: "weak"})

    failures = rows.loc[rows[score_col] == 0].groupby("step").size().rename("failure_count")
    valid = rows.loc[rows[score_col] != 0]
    trajectory = valid.groupby("step")[score_col].agg(["count", "mean", "std"]).reset_index()
    trajectory = trajectory.merge(failures, on="step", how="outer").sort_values("step")
    trajectory["failure_count"] = trajectory["failure_count"].fillna(0).astype(int)
    trajectory["job_name"] = job_name
    trajectory["score_column"] = score_col

    summaries = []
    valid_rows = rows.loc[rows[score_col] != 0]
    for descriptor in selected:
        values = valid_rows[descriptor].dropna()
        summaries.append({
            "job_name": job_name, "descriptor": descriptor, "count": len(values),
            "mean": values.mean(), "std": values.std(), "min": values.min(),
            "median": values.median(), "max": values.max(),
        })
    return rows, trajectory, pd.DataFrame(summaries)
