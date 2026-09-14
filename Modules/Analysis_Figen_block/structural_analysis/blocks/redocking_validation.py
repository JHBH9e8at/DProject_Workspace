"""Load, validate, and summarize PR/PPS redocking validation data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


POSE_COLUMNS = ["state", "pose_index", "docking_score", "rmsd", "within_2A"]
SUMMARY_COLUMNS = [
    "state",
    "pose_count",
    "docking_score_min",
    "docking_score_max",
    "docking_score_mean",
    "docking_score_median",
    "rmsd_min",
    "rmsd_max",
    "rmsd_mean",
    "rmsd_median",
    "rmsd_within_2A_count",
    "rmsd_within_2A_fraction",
]


def load_state_dataset(
    state: str,
    docking_csv: str | Path,
    rmsd_csv: str | Path,
) -> pd.DataFrame:
    """Load one state's score/RMSD pair into a validated pose-level table."""

    normalized_state = state.strip().upper()
    if normalized_state not in {"PR", "PPS"}:
        raise ValueError(f"state must be PR or PPS, received {state!r}")
    docking_path = Path(docking_csv).expanduser().resolve(strict=True)
    rmsd_path = Path(rmsd_csv).expanduser().resolve(strict=True)
    docking = pd.read_csv(docking_path)
    rmsd = pd.read_csv(rmsd_path)
    docking.columns = docking.columns.str.strip()
    rmsd.columns = (
        rmsd.columns.str.strip().str.replace('"', "", regex=False).str.strip()
    )
    if "docking score" not in docking.columns:
        raise ValueError(f"Missing 'docking score' column: {docking_path}")
    if "RMS" not in rmsd.columns:
        raise ValueError(f"Missing 'RMS' column: {rmsd_path}")
    scores = pd.to_numeric(docking["docking score"], errors="raise")
    rmsd_values = pd.to_numeric(
        rmsd["RMS"].astype(str).str.replace('"', "", regex=False).str.strip(),
        errors="raise",
    )
    if len(scores) != len(rmsd_values):
        raise ValueError(
            f"{normalized_state}: docking ({len(scores)}) and RMSD "
            f"({len(rmsd_values)}) row counts differ"
        )
    if len(scores) == 0:
        raise ValueError(f"{normalized_state}: redocking dataset is empty")
    if not np.isfinite(scores.to_numpy(dtype=float)).all():
        raise ValueError(f"{normalized_state}: docking scores contain nonfinite values")
    if not np.isfinite(rmsd_values.to_numpy(dtype=float)).all():
        raise ValueError(f"{normalized_state}: RMSD values contain nonfinite values")
    # The source tables are paired by row position; both must describe the same
    # pose sequence in the same order.
    return pd.DataFrame(
        {
            "state": normalized_state,
            "pose_index": np.arange(1, len(scores) + 1, dtype=int),
            "docking_score": scores.to_numpy(dtype=float),
            "rmsd": rmsd_values.to_numpy(dtype=float),
            # Redocking-success indicator: I_i = 1 when RMSD_i < 2.0 Å.
            "within_2A": rmsd_values.to_numpy(dtype=float) < 2.0,
        },
        columns=POSE_COLUMNS,
    )


def build_redocking_tables(
    *,
    pps_docking_csv: str | Path,
    pps_rmsd_csv: str | Path,
    pr_docking_csv: str | Path,
    pr_rmsd_csv: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return a combined pose table and state-level summary."""

    pps = load_state_dataset("PPS", pps_docking_csv, pps_rmsd_csv)
    pr = load_state_dataset("PR", pr_docking_csv, pr_rmsd_csv)
    if len(pps) != len(pr):
        raise ValueError(f"PPS ({len(pps)}) and PR ({len(pr)}) row counts differ")
    poses = pd.concat((pps, pr), ignore_index=True)
    summary_rows: list[dict[str, float | int | str]] = []
    for state in ("PPS", "PR"):
        data = poses.loc[poses["state"] == state]
        summary_rows.append(
            {
                "state": state,
                "pose_count": len(data),
                "docking_score_min": float(data["docking_score"].min()),
                "docking_score_max": float(data["docking_score"].max()),
                "docking_score_mean": float(data["docking_score"].mean()),
                "docking_score_median": float(data["docking_score"].median()),
                "rmsd_min": float(data["rmsd"].min()),
                "rmsd_max": float(data["rmsd"].max()),
                "rmsd_mean": float(data["rmsd"].mean()),
                "rmsd_median": float(data["rmsd"].median()),
                "rmsd_within_2A_count": int(data["within_2A"].sum()),
                # Boolean mean = N(RMSD < 2 Å) / N(total poses).
                "rmsd_within_2A_fraction": float(data["within_2A"].mean()),
            }
        )
    return poses, pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS)
