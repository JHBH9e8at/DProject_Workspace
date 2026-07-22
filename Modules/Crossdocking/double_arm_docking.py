"""Run both directions of the myosin-state cross-docking experiment."""

import argparse
from pathlib import Path

import pandas as pd

from single_arm_docking import run_pipeline


ARMS = (
    {
        "arm": "L_PPS_to_R_PR",
        "ligand_state": "PPS",
        "receptor_state": "PR",
        "input_key": "pps_input",
        "run_dir_key": "pps_run_dir",
        "run_name_key": "pps_run_name",
        "variant_col_key": "pps_variant_col",
        "grid_key": "pr_grid",
    },
    {
        "arm": "L_PR_to_R_PPS",
        "ligand_state": "PR",
        "receptor_state": "PPS",
        "input_key": "pr_input",
        "run_dir_key": "pr_run_dir",
        "run_name_key": "pr_run_name",
        "variant_col_key": "pr_variant_col",
        "grid_key": "pps_grid",
    },
)


def label_arm_table(table, arm):
    """Add collision-safe state and arm identifiers to an output table."""
    table = table.copy()
    table.insert(0, "crossdock_arm", arm["arm"])
    table.insert(1, "ligand_state", arm["ligand_state"])
    table.insert(2, "receptor_state", arm["receptor_state"])
    if "molecule_id" in table.columns:
        table.insert(
            3,
            "crossdock_id",
            arm["arm"] + "__" + table["molecule_id"].astype(str),
        )
    return table


def combine_arm_outputs(output_dir):
    """Combine both arm manifests and score tables without changing source IDs."""
    manifests = []
    scores = []
    best_scores = []

    for arm in ARMS:
        arm_dir = output_dir / arm["arm"]
        manifest_file = arm_dir / "target_manifest.csv"
        score_file = arm_dir / "combined_docking_scores.csv"
        best_score_file = arm_dir / "best_docking_scores.csv"

        if manifest_file.is_file():
            manifests.append(label_arm_table(pd.read_csv(manifest_file), arm))
        if score_file.is_file():
            scores.append(label_arm_table(pd.read_csv(score_file), arm))
        if best_score_file.is_file():
            best_scores.append(label_arm_table(pd.read_csv(best_score_file), arm))

    if manifests:
        combined_manifest = pd.concat(manifests, ignore_index=True)
        combined_manifest.to_csv(output_dir / "double_arm_manifest.csv", index=False)

        summary = (
            combined_manifest.groupby(
                ["crossdock_arm", "ligand_state", "receptor_state", "pipeline_status"],
                dropna=False,
            )
            .size()
            .rename("molecule_count")
            .reset_index()
        )
        summary.to_csv(output_dir / "double_arm_summary.csv", index=False)

    if scores:
        pd.concat(scores, ignore_index=True).to_csv(
            output_dir / "double_arm_docking_scores.csv", index=False
        )
    if best_scores:
        pd.concat(best_scores, ignore_index=True).to_csv(
            output_dir / "double_arm_best_docking_scores.csv", index=False
        )


def run_double_arm(
    pps_input,
    pps_run_dir,
    pr_input,
    pr_run_dir,
    pps_grid,
    pr_grid,
    output_dir,
    pps_run_name="PPS",
    pr_run_name="PR",
    step_col="step",
    pps_variant_col=None,
    pr_variant_col=None,
    limit=None,
    precision="SP",
    allow_unavailable=False,
    resume=False,
    fail_fast=False,
    ligand_prep_mode="off",
    smiles_col="smiles",
):
    """Run PPS-to-PR and PR-to-PPS arms sequentially and combine outputs."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    values = {
        "pps_input": pps_input,
        "pps_run_dir": pps_run_dir,
        "pps_run_name": pps_run_name,
        "pps_variant_col": pps_variant_col,
        "pps_grid": pps_grid,
        "pr_input": pr_input,
        "pr_run_dir": pr_run_dir,
        "pr_run_name": pr_run_name,
        "pr_variant_col": pr_variant_col,
        "pr_grid": pr_grid,
    }
    arm_errors = []

    for number, arm in enumerate(ARMS, start=1):
        print("\n" + "=" * 72)
        print(f"ARM {number}/2: {arm['arm']}")
        print("=" * 72)
        try:
            run_pipeline(
                input_csv=values[arm["input_key"]],
                run_dir=values[arm["run_dir_key"]],
                run_name=values[arm["run_name_key"]],
                grid_file=values[arm["grid_key"]],
                output_dir=output_dir / arm["arm"],
                step_col=step_col,
                variant_col=values[arm["variant_col_key"]],
                limit=limit,
                precision=precision,
                allow_unavailable=allow_unavailable,
                resume=resume,
                fail_fast=fail_fast,
                ligand_prep_mode=ligand_prep_mode,
                smiles_col=smiles_col,
            )
        except Exception as error:
            arm_errors.append((arm["arm"], error))
            print(f"ARM FAILED: {arm['arm']}: {type(error).__name__}: {error}")
            if fail_fast:
                break
        finally:
            combine_arm_outputs(output_dir)

    combine_arm_outputs(output_dir)

    print("\nDouble-arm cross-docking finished")
    print(f"Output directory: {output_dir}")
    print(f"Arm failures: {len(arm_errors)}")

    if arm_errors:
        details = "; ".join(f"{arm}: {error}" for arm, error in arm_errors)
        raise RuntimeError(f"One or more cross-docking arms failed. {details}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run PPS-to-PR and PR-to-PPS cross-docking sequentially"
    )
    parser.add_argument("--pps-input", required=True, help="Filtered PPS >>>CSV<<<")
    parser.add_argument("--pps-run-dir", required=True, help="PPS AHC result >>directory<<")
    parser.add_argument("--pr-input", required=True, help="Filtered PR >>CSV<<")
    parser.add_argument("--pr-run-dir", required=True, help="PR AHC result >>directory<<")
    parser.add_argument("--pps-grid", required=True, help="PPS receptor grid ZIP")
    parser.add_argument("--pr-grid", required=True, help="PR receptor grid ZIP")
    parser.add_argument("--output-dir", required=True, help="Double-arm result >>directory<<")
    parser.add_argument("--pps-run-name", default="PPS")
    parser.add_argument("--pr-run-name", default="PR")
    parser.add_argument("--step-col", default="step")
    parser.add_argument("--pps-variant-col", default=None)
    parser.add_argument("--pr-variant-col", default=None)
    parser.add_argument("--smiles-col", default="smiles")
    parser.add_argument("--limit", type=int, default=None, help="Per-arm testing limit")
    parser.add_argument(
        "--precision", default="SP", choices=["HTVS", "SP", "XP"] help = " use SP in most cases !")
    parser.add_argument("--allow-unavailable", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument(
        "--ligand-prep-mode",
        default="off",
        choices=["off", "on"],
        help=(
            "off: reuse existing AHC *_lib.sdfgz poses; "
            "on: run LigPrep from each filtered CSV SMILES"
        ),
    )
    args = parser.parse_args()

    run_double_arm(
        pps_input=args.pps_input,
        pps_run_dir=args.pps_run_dir,
        pr_input=args.pr_input,
        pr_run_dir=args.pr_run_dir,
        pps_grid=args.pps_grid,
        pr_grid=args.pr_grid,
        output_dir=args.output_dir,
        pps_run_name=args.pps_run_name,
        pr_run_name=args.pr_run_name,
        step_col=args.step_col,
        pps_variant_col=args.pps_variant_col,
        pr_variant_col=args.pr_variant_col,
        limit=args.limit,
        precision=args.precision,
        allow_unavailable=args.allow_unavailable,
        resume=args.resume,
        fail_fast=args.fail_fast,
        ligand_prep_mode=args.ligand_prep_mode,
        smiles_col=args.smiles_col,
    )
