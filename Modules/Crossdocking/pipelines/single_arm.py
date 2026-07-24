"""Run one arm of the terminal, non-Dask cross-docking pipeline."""

import argparse
import re
from pathlib import Path

import pandas as pd

from docking.glide_input import make_glide_input
from docking.glide_execution import run_glide
from docking.scores import extract_scores
from prep.fetch_targets import fetch_targets
from prep.prepare_ligands import prepare_ligand
from prep.unpack_sdfgz import unpack_sdfgz


class LigandPreparationError(RuntimeError):
    """Raised when the LigPrep stage fails before docking starts."""


def safe_name(value):
    """Convert a molecule ID into a portable filename component."""
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return name or "molecule"


def count_sdf_records(sdf_file):
    """Count complete records in an uncompressed SDF."""
    with Path(sdf_file).open("r", errors="replace") as handle:
        return sum(line.strip() == "$$$$" for line in handle)


def write_pipeline_manifest(manifest, manifest_file):
    manifest.to_csv(manifest_file, index=False)


def add_score_metadata(scores, row, selection_index, step_col):
    scores = scores.copy()
    scores.insert(0, "selection_index", selection_index)
    scores.insert(1, "molecule_id", row["molecule_id"])
    scores.insert(2, "source_run", row["source_run"])
    scores.insert(3, "source_step", row[step_col])
    scores.insert(4, "ligand_prep_mode", row["ligand_prep_mode"])
    scores.insert(5, "prepared_variant_count", row["prepared_variant_count"])
    scores.insert(6, "docked_pose_count", len(scores))
    return scores


def collect_completed_scores(manifest, step_col):
    """Rebuild the all-pose score table from completed per-job CSV files."""
    tables = []
    for index, row in manifest.iterrows():
        if row.get("pipeline_status") != "completed":
            continue
        score_file = Path(str(row["score_csv"]))
        if not score_file.is_file():
            continue
        scores = pd.read_csv(score_file)
        tables.append(add_score_metadata(scores, row, index + 1, step_col))
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def write_combined_scores(manifest, all_pose_file, best_score_file, step_col):
    combined = collect_completed_scores(manifest, step_col)
    if combined.empty:
        return combined

    combined.to_csv(all_pose_file, index=False)
    best = (
        combined.sort_values("r_i_docking_score", ascending=True)
        .groupby("selection_index", as_index=False, sort=False)
        .head(1)
        .sort_values("selection_index")
    )
    best.to_csv(best_score_file, index=False)
    return combined


def prepare_existing_pose(row, job_dir, resume):
    source_sdfgz = Path(row["target_sdfgz"])
    input_dir = job_dir / "input"
    expected_sdf = input_dir / source_sdfgz.with_suffix(".sdf").name
    if resume and expected_sdf.is_file():
        print(f"Reusing unpacked ligand: {expected_sdf}")
        return source_sdfgz, expected_sdf
    prepared_sdf = unpack_sdfgz(
        input_file=source_sdfgz,
        output_dir=input_dir,
        log_dir=job_dir / "logs",
    )
    return source_sdfgz, prepared_sdf


def prepare_from_smiles(row, job_dir, smiles_col, resume):
    molecule_id = safe_name(row["molecule_id"])
    smiles = str(row[smiles_col]).strip()
    input_smi = job_dir / "input" / f"{molecule_id}.smi"
    prepared_sdf = job_dir / "ligprep" / f"{molecule_id}_prepared.sdf"

    input_smi.parent.mkdir(parents=True, exist_ok=True)
    if not (resume and input_smi.is_file()):
        input_smi.write_text(f"{smiles}\t{molecule_id}\n", encoding="utf-8")

    if resume and prepared_sdf.is_file() and prepared_sdf.stat().st_size > 0:
        print(f"Reusing LigPrep output: {prepared_sdf}")
        return input_smi, prepared_sdf

    return input_smi, prepare_ligand(input_smi, prepared_sdf)


def write_job_best_score(scores, best_score_file):
    best = scores.nsmallest(1, "r_i_docking_score")
    best.to_csv(best_score_file, index=False)
    return float(best.iloc[0]["r_i_docking_score"])


def run_one_target(
    row,
    selection_index,
    grid_file,
    jobs_dir,
    precision,
    resume,
    ligand_prep_mode,
    smiles_col,
):
    molecule_id = str(row["molecule_id"])
    job_dir = jobs_dir / f"{selection_index:04d}_{safe_name(molecule_id)}"
    score_file = job_dir / "docking_scores.csv"
    best_score_file = job_dir / "best_docking_score.csv"

    if resume and score_file.is_file():
        print(f"Already docked, reusing: {score_file}")
        scores = pd.read_csv(score_file)
        best_score = write_job_best_score(scores, best_score_file)
        pose_files = sorted(job_dir.glob("*_lib.sdfgz"))
        prepared_files = sorted((job_dir / "ligprep").glob("*_prepared.sdf"))
        if ligand_prep_mode == "off":
            prepared_files = sorted((job_dir / "input").glob("*.sdf"))
        prepared_sdf = prepared_files[0] if prepared_files else ""
        variant_count = count_sdf_records(prepared_sdf) if prepared_sdf else 0
        if ligand_prep_mode == "on":
            source_ligand = job_dir / "input" / f"{safe_name(molecule_id)}.smi"
        else:
            source_ligand = Path(row["target_sdfgz"])
        return {
            "job_dir": job_dir,
            "pose_file": pose_files[0] if pose_files else "",
            "score_file": score_file,
            "best_score_file": best_score_file,
            "best_score": best_score,
            "source_ligand_file": source_ligand,
            "prepared_ligand_file": prepared_sdf,
            "prepared_variant_count": variant_count,
        }

    if job_dir.exists() and any(job_dir.iterdir()) and not resume:
        raise FileExistsError(
            f"Job directory is not empty: {job_dir}. Use --resume or a new "
            "--output-dir."
        )

    job_dir.mkdir(parents=True, exist_ok=True)
    if ligand_prep_mode == "on":
        try:
            source_ligand, prepared_sdf = prepare_from_smiles(
                row, job_dir, smiles_col, resume
            )
        except Exception as error:
            raise LigandPreparationError(str(error)) from error
    else:
        source_ligand, prepared_sdf = prepare_existing_pose(row, job_dir, resume)

    variant_count = count_sdf_records(prepared_sdf)
    if variant_count == 0:
        message = f"Prepared ligand contains no complete SDF records: {prepared_sdf}"
        if ligand_prep_mode == "on":
            raise LigandPreparationError(message)
        raise RuntimeError(message)

    pose_files = sorted(job_dir.glob("*_lib.sdfgz"))
    if resume and pose_files:
        pose_file = pose_files[0]
        print(f"Reusing Glide pose file: {pose_file}")
    else:
        glide_input = make_glide_input(
            prepared_ligand=prepared_sdf,
            grid_file=grid_file,
            output_in=job_dir / "docking.in",
            precision=precision,
        )
        pose_file = run_glide(glide_input=glide_input, output_dir=job_dir)

    scores = extract_scores(pose_file=pose_file, output_csv=score_file)
    best_score = write_job_best_score(scores, best_score_file)
    return {
        "job_dir": job_dir,
        "pose_file": pose_file,
        "score_file": score_file,
        "best_score_file": best_score_file,
        "best_score": best_score,
        "source_ligand_file": source_ligand,
        "prepared_ligand_file": prepared_sdf,
        "prepared_variant_count": variant_count,
    }


def run_pipeline(
    input_csv,
    run_dir,
    run_name,
    grid_file,
    output_dir,
    step_col="step",
    variant_col=None,
    limit=None,
    precision="SP",
    allow_unavailable=False,
    resume=False,
    fail_fast=False,
    ligand_prep_mode="off",
    smiles_col="smiles",
):
    ligand_prep_mode = ligand_prep_mode.lower()
    if ligand_prep_mode not in {"on", "off"}:
        raise ValueError("ligand_prep_mode must be 'on' or 'off'")

    output_dir = Path(output_dir).resolve()
    grid_file = Path(grid_file).resolve()
    if not grid_file.is_file():
        raise FileNotFoundError(f"Receptor grid not found: {grid_file}")

    output_dir.mkdir(parents=True, exist_ok=True)
    jobs_dir = output_dir / "jobs"
    jobs_dir.mkdir(exist_ok=True)
    manifest_file = output_dir / "target_manifest.csv"
    combined_file = output_dir / "combined_docking_scores.csv"
    best_file = output_dir / "best_docking_scores.csv"

    manifest = fetch_targets(
        input_csv=input_csv,
        run_dir=run_dir,
        run_name=run_name,
        output_csv=manifest_file,
        step_col=step_col,
        variant_col=variant_col,
        limit=limit,
        # Existing pose availability is irrelevant when LigPrep starts from SMILES.
        allow_unavailable=allow_unavailable or ligand_prep_mode == "on",
    )

    if ligand_prep_mode == "on":
        if smiles_col not in manifest.columns:
            raise ValueError(f"Missing SMILES column: {smiles_col}")
        smiles_ready = manifest[smiles_col].notna() & manifest[smiles_col].astype(str).str.strip().ne("")
        manifest["pipeline_ready"] = smiles_ready
        missing_input_count = int((~smiles_ready).sum())
        if missing_input_count and not allow_unavailable:
            raise ValueError(
                f"{missing_input_count} selected molecule(s) have no SMILES in {smiles_col}"
            )
    else:
        manifest["pipeline_ready"] = manifest["target_ready"]

    manifest["pipeline_status"] = manifest["pipeline_ready"].map(
        {True: "pending", False: "unavailable"}
    )
    manifest["pipeline_error"] = ""
    if ligand_prep_mode == "on":
        manifest.loc[~manifest["pipeline_ready"], "pipeline_error"] = "Missing SMILES"
    manifest["ligand_prep_mode"] = ligand_prep_mode
    manifest["input_smiles"] = manifest[smiles_col] if ligand_prep_mode == "on" else ""
    manifest["ligprep_status"] = "pending" if ligand_prep_mode == "on" else "not_requested"
    manifest.loc[~manifest["pipeline_ready"], "ligprep_status"] = "unavailable"
    manifest["source_ligand_file"] = ""
    manifest["prepared_ligand_file"] = ""
    manifest["prepared_variant_count"] = 0
    manifest["job_dir"] = ""
    manifest["pose_file"] = ""
    manifest["score_csv"] = ""
    manifest["best_score_csv"] = ""
    manifest["best_docking_score"] = pd.NA
    write_pipeline_manifest(manifest, manifest_file)

    failures = []
    ready_indices = manifest.index[manifest["pipeline_ready"]].tolist()
    total_ready = len(ready_indices)

    for count, index in enumerate(ready_indices, start=1):
        row = manifest.loc[index]
        molecule_id = row["molecule_id"]
        expected_job_dir = jobs_dir / f"{index + 1:04d}_{safe_name(molecule_id)}"
        manifest.at[index, "job_dir"] = str(expected_job_dir)
        manifest.at[index, "score_csv"] = str(expected_job_dir / "docking_scores.csv")
        manifest.at[index, "best_score_csv"] = str(expected_job_dir / "best_docking_score.csv")
        print(f"\n[{count}/{total_ready}] Docking {molecule_id} (LigPrep: {ligand_prep_mode})")
        try:
            result = run_one_target(
                row=row,
                selection_index=index + 1,
                grid_file=grid_file,
                jobs_dir=jobs_dir,
                precision=precision,
                resume=resume,
                ligand_prep_mode=ligand_prep_mode,
                smiles_col=smiles_col,
            )
            manifest.at[index, "pipeline_status"] = "completed"
            manifest.at[index, "ligprep_status"] = (
                "completed" if ligand_prep_mode == "on" else "not_requested"
            )
            manifest.at[index, "job_dir"] = str(result["job_dir"])
            manifest.at[index, "pose_file"] = str(result["pose_file"])
            manifest.at[index, "score_csv"] = str(result["score_file"])
            manifest.at[index, "best_score_csv"] = str(result["best_score_file"])
            manifest.at[index, "best_docking_score"] = result["best_score"]
            manifest.at[index, "source_ligand_file"] = str(result["source_ligand_file"])
            manifest.at[index, "prepared_ligand_file"] = str(result["prepared_ligand_file"])
            manifest.at[index, "prepared_variant_count"] = result["prepared_variant_count"]
        except Exception as error:
            manifest.at[index, "pipeline_status"] = "failed"
            if ligand_prep_mode == "on":
                manifest.at[index, "ligprep_status"] = (
                    "failed" if isinstance(error, LigandPreparationError) else "completed"
                )
            manifest.at[index, "pipeline_error"] = f"{type(error).__name__}: {error}"
            failures.append((molecule_id, error))
            print(f"FAILED {molecule_id}: {type(error).__name__}: {error}")
        finally:
            write_pipeline_manifest(manifest, manifest_file)
            write_combined_scores(manifest, combined_file, best_file, step_col)

        if failures and fail_fast:
            break

    completed = int((manifest["pipeline_status"] == "completed").sum())
    unavailable = int((manifest["pipeline_status"] == "unavailable").sum())
    failed = int((manifest["pipeline_status"] == "failed").sum())
    write_pipeline_manifest(manifest, manifest_file)
    write_combined_scores(manifest, combined_file, best_file, step_col)

    print("\nCross-docking pipeline finished")
    print(f"LigPrep mode: {ligand_prep_mode}")
    print(f"Completed: {completed}")
    print(f"Failed: {failed}")
    print(f"Unavailable: {unavailable}")
    print(f"Manifest: {manifest_file}")
    if combined_file.is_file():
        print(f"All docking poses: {combined_file}")
        print(f"Best score per instance: {best_file}")

    if failures:
        raise RuntimeError(
            f"{len(failures)} docking job(s) failed. See {manifest_file} for details."
        )
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run one sequential AHC cross-docking arm"
    )
    parser.add_argument("--input", required=True, help="Filtered/selected CSV")
    parser.add_argument("--run-dir", required=True, help="AHC result directory")
    parser.add_argument("--run-name", required=True, help="Run prefix, e.g. PPS or PR")
    parser.add_argument("--grid", required=True, help="Opposite-state Glide grid ZIP")
    parser.add_argument("--output-dir", required=True, help="Pipeline result directory")
    parser.add_argument("--step-col", default="step")
    parser.add_argument(
        "--variant-col",
        default=None,
        help="Best-variant column; defaults to <run-name>_best_variant",
    )
    parser.add_argument("--smiles-col", default="smiles")
    parser.add_argument("--limit", type=int, default=None, help="Testing limit")
    parser.add_argument("--precision", default="SP", choices=["HTVS", "SP", "XP"])
    parser.add_argument(
        "--ligand-prep-mode",
        default="off",
        choices=["off", "on"],
        help=(
            "off: unpack and dock the existing AHC *_lib.sdfgz pose; "
            "on: run LigPrep from the filtered CSV SMILES before docking"
        ),
    )
    parser.add_argument(
        "--allow-unavailable",
        action="store_true",
        help="Skip missing source poses (off) or missing SMILES (on)",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    args = parser.parse_args()

    run_pipeline(
        input_csv=args.input,
        run_dir=args.run_dir,
        run_name=args.run_name,
        grid_file=args.grid,
        output_dir=args.output_dir,
        step_col=args.step_col,
        variant_col=args.variant_col,
        limit=args.limit,
        precision=args.precision,
        allow_unavailable=args.allow_unavailable,
        resume=args.resume,
        fail_fast=args.fail_fast,
        ligand_prep_mode=args.ligand_prep_mode,
        smiles_col=args.smiles_col,
    )
# python -m pipelines.single_arm \
# --input /home/andy/proj/Results/temptesting/target.csv \
# --run-dir /home/andy/proj/Results/AHC_Resultsbin/2026_06_02_SMILES-RNN_PPS_tester_run_2 \
# --run-name PPS \
# --grid /home/andy/proj/DProject_Workspace/AHC_Related/ref_structures/grids/PR/PR_Auto_Grid.zip \
# --output-dir /home/andy/proj/Results/temptesting/PPStoPRtest \
# --ligand-prep-mode on \
# --precision SP \
# --limit 5
