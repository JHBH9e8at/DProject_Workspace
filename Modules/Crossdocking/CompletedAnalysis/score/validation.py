"""Load and validate double-arm cross-docking result tables."""

from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_ARMS = {
    "L_PPS_to_R_PR": ("PPS", "PR"),
    "L_PR_to_R_PPS": ("PR", "PPS"),
}

MANIFEST_REQUIRED = {
    "crossdock_id",
    "crossdock_arm",
    "ligand_state",
    "receptor_state",
    "molecule_id",
    "pipeline_status",
}

BEST_REQUIRED = {
    "crossdock_id",
    "crossdock_arm",
    "ligand_state",
    "receptor_state",
    "molecule_id",
    "r_i_docking_score",
}


def load_double_arm_results(input_dir):
    """Load the combined manifest and per-instance best docking scores."""
    input_dir = Path(input_dir).resolve()
    manifest_file = input_dir / "double_arm_manifest.csv"
    best_file = input_dir / "double_arm_best_docking_scores.csv"

    missing = [str(path) for path in (manifest_file, best_file) if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing double-arm result file(s): {missing}")

    return pd.read_csv(manifest_file), pd.read_csv(best_file)


def _issue(rows, severity, check, message, affected_count=0):
    rows.append(
        {
            "severity": severity,
            "check": check,
            "affected_count": int(affected_count),
            "message": message,
        }
    )


def _check_required_columns(table, required, table_name, rows):
    missing = sorted(required - set(table.columns))
    if missing:
        _issue(
            rows,
            "ERROR",
            f"{table_name}_required_columns",
            f"Missing required columns: {missing}",
            len(missing),
        )
        return False
    _issue(rows, "PASS", f"{table_name}_required_columns", "All required columns present")
    return True


def validate_double_arm_results(manifest, best_scores):
    """Return a machine-readable QC report for paired selectivity analysis."""
    rows = []
    manifest_ok = _check_required_columns(
        manifest, MANIFEST_REQUIRED, "manifest", rows
    )
    best_ok = _check_required_columns(best_scores, BEST_REQUIRED, "best_scores", rows)

    if not (manifest_ok and best_ok):
        return pd.DataFrame(rows)

    for table, name in ((manifest, "manifest"), (best_scores, "best_scores")):
        duplicates = int(table["crossdock_id"].duplicated().sum())
        _issue(
            rows,
            "ERROR" if duplicates else "PASS",
            f"{name}_unique_crossdock_id",
            "Duplicate crossdock_id values found" if duplicates else "crossdock_id values are unique",
            duplicates,
        )

    bad_arm_rows = 0
    for _, row in manifest.iterrows():
        expected = EXPECTED_ARMS.get(row["crossdock_arm"])
        if expected is None or expected != (row["ligand_state"], row["receptor_state"]):
            bad_arm_rows += 1
    _issue(
        rows,
        "ERROR" if bad_arm_rows else "PASS",
        "arm_state_mapping",
        "Arm/state labels are inconsistent" if bad_arm_rows else "Both arm/state mappings are consistent",
        bad_arm_rows,
    )

    observed_arms = set(manifest["crossdock_arm"].dropna())
    missing_arms = sorted(set(EXPECTED_ARMS) - observed_arms)
    _issue(
        rows,
        "ERROR" if missing_arms else "PASS",
        "both_arms_present",
        f"Missing expected arms: {missing_arms}" if missing_arms else "Both expected arms are present",
        len(missing_arms),
    )

    completed_ids = set(
        manifest.loc[manifest["pipeline_status"] == "completed", "crossdock_id"]
    )
    best_ids = set(best_scores["crossdock_id"])
    missing_best = completed_ids - best_ids
    unexpected_best = best_ids - completed_ids
    mismatch_count = len(missing_best) + len(unexpected_best)
    _issue(
        rows,
        "ERROR" if mismatch_count else "PASS",
        "completed_to_best_score_reconciliation",
        (
            f"Missing best rows: {len(missing_best)}; unexpected best rows: {len(unexpected_best)}"
            if mismatch_count
            else "Every completed manifest row has exactly one best-score row"
        ),
        mismatch_count,
    )

    numeric_scores = pd.to_numeric(best_scores["r_i_docking_score"], errors="coerce")
    invalid_scores = int((~np.isfinite(numeric_scores)).sum())
    nonnegative_scores = int((numeric_scores >= 0).sum())
    _issue(
        rows,
        "ERROR" if invalid_scores else "PASS",
        "finite_opposite_state_scores",
        "Invalid or missing docking scores found" if invalid_scores else "All opposite-state scores are finite",
        invalid_scores,
    )
    _issue(
        rows,
        "ERROR" if nonnegative_scores else "PASS",
        "negative_opposite_state_scores",
        "Zero or positive docking scores found" if nonnegative_scores else "All opposite-state scores are negative",
        nonnegative_scores,
    )

    missing_own_columns = []
    invalid_own_scores = 0
    for state in ("PPS", "PR"):
        column = f"{state}_r_i_docking_score"
        state_rows = manifest["ligand_state"] == state
        if column not in manifest.columns:
            missing_own_columns.append(column)
            continue
        values = pd.to_numeric(manifest.loc[state_rows, column], errors="coerce")
        invalid_own_scores += int((~np.isfinite(values)).sum())

    _issue(
        rows,
        "ERROR" if missing_own_columns else "PASS",
        "own_state_score_columns",
        (
            f"Missing own-state score columns: {missing_own_columns}"
            if missing_own_columns
            else "PPS and PR own-state score columns are present"
        ),
        len(missing_own_columns),
    )
    _issue(
        rows,
        "ERROR" if invalid_own_scores else "PASS",
        "finite_own_state_scores",
        "Invalid own-state scores found" if invalid_own_scores else "All own-state scores are finite",
        invalid_own_scores,
    )

    if "best_docking_score" in manifest.columns:
        joined = manifest[["crossdock_id", "best_docking_score"]].merge(
            best_scores[["crossdock_id", "r_i_docking_score"]],
            on="crossdock_id",
            how="inner",
        )
        left = pd.to_numeric(joined["best_docking_score"], errors="coerce")
        right = pd.to_numeric(joined["r_i_docking_score"], errors="coerce")
        mismatches = int((left.sub(right).abs() > 1e-8).sum())
        _issue(
            rows,
            "ERROR" if mismatches else "PASS",
            "manifest_best_score_match",
            "Manifest and best-score values differ" if mismatches else "Manifest and best-score values match",
            mismatches,
        )

    return pd.DataFrame(rows)


def qc_has_errors(qc_report):
    return bool((qc_report["severity"] == "ERROR").any())
