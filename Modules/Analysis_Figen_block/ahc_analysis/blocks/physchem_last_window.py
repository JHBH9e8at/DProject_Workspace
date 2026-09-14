"""Stable PR/PPS last-N physicochemical scatter input tables."""

from pathlib import Path
import pandas as pd


SCRIPT_ID = "PY-147"
DESCRIPTORS = ["desc_MolWt", "desc_NumRotatableBonds", "desc_CLogP", "desc_TPSA"]


def prepare_last_window(input_csv, *, job_name, last_n):
    if int(last_n) < 1:
        raise ValueError("last_n must be at least 1")
    score_col = f"{job_name}_r_i_docking_score"
    frame = pd.read_csv(input_csv)
    columns = ["step", score_col, *DESCRIPTORS]
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f"{Path(input_csv)}: missing required columns: {missing}")
    rows = frame[columns].copy()
    for column in columns:
        rows[column] = pd.to_numeric(rows[column], errors="coerce")
    loaded_rows = len(rows)
    rows = rows.dropna(subset=["step", score_col])
    rows = rows[rows[score_col].ne(0)].copy()
    if rows.empty:
        raise ValueError(f"No plottable {job_name} rows remain")
    max_step = rows.step.max()
    first_step = max(rows.step.min(), max_step - int(last_n) + 1)
    rows.insert(0, "population_state", job_name)
    rows["window"] = rows.step.ge(first_step).map({True: "recent", False: "earlier"})
    summary = {"population_state": job_name, "score_column": score_col,
               "loaded_rows": loaded_rows, "plottable_rows": len(rows),
               "recent_rows": int(rows.window.eq("recent").sum()),
               "first_last_step": first_step, "max_step": max_step,
               "last_n": int(last_n)}
    return rows, summary


def prepare_pr_pps(pr_csv, pps_csv, *, last_n):
    pr, pr_summary = prepare_last_window(pr_csv, job_name="PR", last_n=last_n)
    pps, pps_summary = prepare_last_window(pps_csv, job_name="PPS", last_n=last_n)
    return pd.concat((pr, pps), ignore_index=True), pd.DataFrame((pr_summary, pps_summary))
