"""Prepare every selected AHC ligand required by a cross-docking arm."""

import argparse
import re
from pathlib import Path

def _safe_name(value):
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return name or "molecule"


def run_preparation(
    input_csv,
    run_dir,
    run_name,
    output_dir,
    step_col="step",
    variant_col=None,
    best_dscore_top_per=1.0,
    testmode_top_n=None,
    allow_unavailable=False,
    resume=False,
    fail_fast=False,
):
    """Build a target manifest and materialize every available ligand as SDF."""
    from .fetch_targets import fetch_targets
    from .filter_population import filter_population
    from .unpack_sdfgz import unpack_sdfgz

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_file = output_dir / f"{run_name}_filtered.csv"
    manifest_file = output_dir / "preparation_manifest.csv"
    ligands_dir = output_dir / "ligands"
    ligands_dir.mkdir(exist_ok=True)

    filter_population(
        input_csv=input_csv,
        output_csv=filtered_file,
        run_name=run_name,
        best_dscore_top_per=best_dscore_top_per,
        testmode_top_n=testmode_top_n,
    )

    manifest = fetch_targets(
        input_csv=filtered_file,
        run_dir=run_dir,
        run_name=run_name,
        output_csv=manifest_file,
        step_col=step_col,
        variant_col=variant_col,
        limit=None,
        allow_unavailable=allow_unavailable,
    ).copy()
    manifest["prepared_sdf"] = ""
    manifest["preparation_status"] = "pending"
    manifest["preparation_error"] = ""

    for index, row in manifest.iterrows():
        if not bool(row["target_ready"]):
            manifest.at[index, "preparation_status"] = "unavailable"
            continue

        ligand_dir = ligands_dir / f"{index + 1:04d}_{_safe_name(row['molecule_id'])}"
        source = Path(row["target_sdfgz"])
        expected_sdf = ligand_dir / source.with_suffix(".sdf").name

        try:
            if resume and expected_sdf.is_file() and expected_sdf.stat().st_size > 0:
                prepared_sdf = expected_sdf
                status = "reused"
            else:
                prepared_sdf = unpack_sdfgz(
                    input_file=source,
                    output_dir=ligand_dir,
                    log_dir=output_dir / "logs",
                )
                status = "completed"
            manifest.at[index, "prepared_sdf"] = str(prepared_sdf)
            manifest.at[index, "preparation_status"] = status
        except Exception as error:
            manifest.at[index, "preparation_status"] = "failed"
            manifest.at[index, "preparation_error"] = str(error)
            manifest.to_csv(manifest_file, index=False)
            if fail_fast:
                raise

        manifest.to_csv(manifest_file, index=False)

    completed = manifest["preparation_status"].isin(["completed", "reused"])
    failed = manifest["preparation_status"].eq("failed")
    print("Preparation block completed")
    print(f"Prepared ligands: {int(completed.sum())}")
    print(f"Failed ligands: {int(failed.sum())}")
    print(f"Manifest: {manifest_file}")

    if failed.any():
        raise RuntimeError(
            f"{int(failed.sum())} ligand(s) failed preparation; see {manifest_file}"
        )
    return manifest


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Raw AHC scores CSV")
    parser.add_argument("--run-dir", required=True, help="AHC result directory")
    parser.add_argument("--run-name", required=True, help="Run prefix, e.g. PPS or PR")
    parser.add_argument("--output-dir", required=True, help="Preparation result directory")
    parser.add_argument("--step-col", default="step")
    parser.add_argument("--variant-col")
    parser.add_argument(
        "--best-dscore-top-per",
        type=float,
        default=1.0,
        help="Top own-state docking-score percentage selected after filtering",
    )
    parser.add_argument(
        "--testmode-top-n",
        type=int,
        help="Testing override: select exactly the top N instead of a percentage",
    )
    parser.add_argument("--allow-unavailable", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_preparation(
        input_csv=args.input,
        run_dir=args.run_dir,
        run_name=args.run_name,
        output_dir=args.output_dir,
        step_col=args.step_col,
        variant_col=args.variant_col,
        best_dscore_top_per=args.best_dscore_top_per,
        testmode_top_n=args.testmode_top_n,
        allow_unavailable=args.allow_unavailable,
        resume=args.resume,
        fail_fast=args.fail_fast,
    )


if __name__ == "__main__":
    main()
