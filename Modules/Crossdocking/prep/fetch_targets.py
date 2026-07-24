"""Locate selected AHC Glide pose files and write a cross-docking manifest."""

import argparse
from pathlib import Path

import pandas as pd


def build_target_records(
    input_csv,
    run_dir,
    run_name,
    step_col="step",
    variant_col=None,
    limit=None,
):
    df = pd.read_csv(input_csv)

    if variant_col is None:
        variant_col = f"{run_name}_best_variant"

    required_cols = [step_col, variant_col]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required column(s): {missing_cols}")

    if limit is not None:
        df = df.head(limit).copy()
    else:
        df = df.copy()

    if df.empty:
        raise ValueError("Input CSV contains no selected molecules")

    if df[variant_col].isna().any():
        raise ValueError(f"Missing values found in {variant_col}")
    if df[variant_col].duplicated().any():
        duplicates = df.loc[df[variant_col].duplicated(), variant_col].tolist()
        raise ValueError(f"Duplicate best-variant IDs found: {duplicates[:5]}")

    run_dir = Path(run_dir).resolve()
    glide_dir = run_dir / f"{run_name}_GlideDock"

    records = []
    for _, row in df.iterrows():
        step = int(row[step_col])
        variant = str(row[variant_col])
        target_file = glide_dir / str(step) / f"{variant}_lib.sdfgz"

        record = row.to_dict()
        record["molecule_id"] = variant
        record["source_run"] = run_name
        record["target_sdfgz"] = str(target_file)
        records.append(record)

    return records


def inspect_target(record):
    """Check one target path. Kept separate so Dask can reuse it."""
    target_file = Path(record["target_sdfgz"])
    exists = target_file.is_file()
    size = target_file.stat().st_size if exists else 0

    checked = record.copy()
    checked["target_exists"] = exists
    checked["target_size_bytes"] = size
    checked["target_ready"] = exists and size > 0
    return checked


def write_manifest(records, output_csv, allow_unavailable=False):
    manifest = pd.DataFrame(records)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_csv, index=False)

    ready_count = int(manifest["target_ready"].sum())
    unavailable_count = len(manifest) - ready_count

    print(f"Manifest: {output_csv.resolve()}")
    print(f"Selected targets: {len(manifest)}")
    print(f"Ready targets: {ready_count}")
    print(f"Unavailable targets: {unavailable_count}")

    if unavailable_count and not allow_unavailable:
        raise FileNotFoundError(
            f"{unavailable_count} target file(s) are missing or empty. "
            "The manifest was written for inspection."
        )

    return manifest


def fetch_targets(
    input_csv,
    run_dir,
    run_name,
    output_csv,
    step_col="step",
    variant_col=None,
    limit=None,
    allow_unavailable=False,
):
    records = build_target_records(
        input_csv=input_csv,
        run_dir=run_dir,
        run_name=run_name,
        step_col=step_col,
        variant_col=variant_col,
        limit=limit,
    )
    checked_records = [inspect_target(record) for record in records]
    return write_manifest(checked_records, output_csv, allow_unavailable)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Locate selected *_lib.sdfgz files and write a manifest"
    )
    parser.add_argument("--input", required=True, help="Filtered/selected CSV")
    parser.add_argument("--run-dir", required=True, help="AHC result directory")
    parser.add_argument("--run-name", required=True, help="Run prefix, e.g. PPS or PR")
    parser.add_argument("--output", required=True, help="Output manifest CSV")
    parser.add_argument("--step-col", default="step")
    parser.add_argument(
        "--variant-col",
        default=None,
        help="Best-variant column; defaults to <run-name>_best_variant",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="FOR TEST ONLY!:::normal filtering output should already be selected",
    )
    parser.add_argument(
        "--allow-unavailable",
        action="store_true",
        help="Do not fail when target files are missing or empty Should not really be an issue",
    )
    args = parser.parse_args()

    fetch_targets(
        input_csv=args.input,
        run_dir=args.run_dir,
        run_name=args.run_name,
        output_csv=args.output,
        step_col=args.step_col,
        variant_col=args.variant_col,
        limit=args.limit,
        allow_unavailable=args.allow_unavailable,
    )




#### exampel

# For linux
# python fetch_target.py \
# --input Q:\coding_dir\701_Project\workspace\Modules\Crossdocking\Batchdocking\T1p_PPS.csv \
# --run-dir "Q:\coding_dir\Workspace1\Docking\2026_06_02_SMILES-RNN_PPS_tester_run_2" \
# --run-name PPS \
# --output Q:\coding_dir\701_Project\workspace\Modules\Crossdocking\Batchdocking\target.csv \
# --allow-unavailable

#For PS [ oinly for local testing]
# python fetch_target.py `
#   --input "Q:\coding_dir\701_Project\workspace\Modules\Crossdocking\patching\testingbin\T1p_PPS.csv" `
#   --run-dir "Q:\coding_dir\Workspace1\Docking\2026_06_02_SMILES-RNN_PPS_tester_run_2" `
#   --run-name PPS `
#   --output "Q:\coding_dir\701_Project\workspace\Modules\Crossdocking\patching\testingbin\target.csv" `
#   --allow-unavailable


#   --limit 100