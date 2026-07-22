"""Run the terminal, non-Dask cross-docking pipeline."""

import argparse
import re
from pathlib import Path

import pandas as pd

from s1_Prep.fetch_target import fetch_targets
from s1_Prep.unpack_p201 import unpack_sdfgz
from s2_Dock.extract_scores import extract_scores
from s2_Dock.make_glide_input import make_glide_input
from s2_Dock.run_glide import run_glide


def safe_name(value):
    """Convert a molecule ID into a portable directory name."""
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return name or "molecule"


def write_pipeline_manifest(manifest, manifest_file):
    manifest.to_csv(manifest_file, index=False)


def add_score_metadata(scores, row, selection_index, step_col):
    scores = scores.copy()
    scores.insert(0, "selection_index", selection_index)
    scores.insert(1, "molecule_id", row["molecule_id"])
    scores.insert(2, "source_run", row["source_run"])
    scores.insert(3, "source_step", row[step_col])
    return scores


def collect_completed_scores(manifest, step_col):
    """Rebuild the combined score table from completed per-job CSV files."""
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


def write_combined_scores(manifest, output_file, step_col):
    combined = collect_completed_scores(manifest, step_col)
    if not combined.empty:
        combined.to_csv(output_file, index=False)
    return combined


def prepare_job_ligand(source_sdfgz, job_dir, resume):
    input_dir = job_dir / "input"
    expected_sdf = input_dir / Path(source_sdfgz).with_suffix(".sdf").name
    if resume and expected_sdf.is_file():
        print(f"Reusing unpacked ligand: {expected_sdf}")
        return expected_sdf
    return unpack_sdfgz(
        input_file=source_sdfgz,
        output_dir=input_dir,
        log_dir=job_dir / "logs",
    )


def run_one_target(row, selection_index, grid_file, jobs_dir, precision, resume):
    molecule_id = str(row["molecule_id"])
    job_dir = jobs_dir / f"{selection_index:04d}_{safe_name(molecule_id)}"
    score_file = job_dir / "docking_scores.csv"

    if resume and score_file.is_file():
        print(f"Already completed, reusing: {score_file}")
        pose_files = sorted(job_dir.glob("*_lib.sdfgz"))
        return job_dir, pose_files[0] if pose_files else "", score_file

    if job_dir.exists() and any(job_dir.iterdir()) and not resume:
        raise FileExistsError(
            f"Job directory is not empty: {job_dir}. Use --resume or a new "
            "--output-dir."
        )

    job_dir.mkdir(parents=True, exist_ok=True)
    prepared_sdf = prepare_job_ligand(row["target_sdfgz"], job_dir, resume)

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

    extract_scores(pose_file=pose_file, output_csv=score_file)
    return job_dir, pose_file, score_file


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
):
    output_dir = Path(output_dir).resolve()
    grid_file = Path(grid_file).resolve()
    if not grid_file.is_file():
        raise FileNotFoundError(f"Receptor grid not found: {grid_file}")

    output_dir.mkdir(parents=True, exist_ok=True)
    jobs_dir = output_dir / "jobs"
    jobs_dir.mkdir(exist_ok=True)
    manifest_file = output_dir / "target_manifest.csv"
    combined_file = output_dir / "combined_docking_scores.csv"

    manifest = fetch_targets(
        input_csv=input_csv,
        run_dir=run_dir,
        run_name=run_name,
        output_csv=manifest_file,
        step_col=step_col,
        variant_col=variant_col,
        limit=limit,
        allow_unavailable=allow_unavailable,
    )
    manifest["pipeline_status"] = manifest["target_ready"].map(
        {True: "pending", False: "unavailable"}
    )
    manifest["pipeline_error"] = ""
    manifest["job_dir"] = ""
    manifest["pose_file"] = ""
    manifest["score_csv"] = ""
    write_pipeline_manifest(manifest, manifest_file)

    failures = []
    ready_indices = manifest.index[manifest["target_ready"]].tolist()
    total_ready = len(ready_indices)

    for count, index in enumerate(ready_indices, start=1):
        row = manifest.loc[index]
        molecule_id = row["molecule_id"]
        expected_job_dir = jobs_dir / f"{index + 1:04d}_{safe_name(molecule_id)}"
        manifest.at[index, "job_dir"] = str(expected_job_dir)
        manifest.at[index, "score_csv"] = str(expected_job_dir / "docking_scores.csv")
        print(f"\n[{count}/{total_ready}] Docking {molecule_id}")
        try:
            job_dir, pose_file, score_file = run_one_target(
                row=row,
                selection_index=index + 1,
                grid_file=grid_file,
                jobs_dir=jobs_dir,
                precision=precision,
                resume=resume,
            )
            manifest.at[index, "pipeline_status"] = "completed"
            manifest.at[index, "job_dir"] = str(job_dir)
            manifest.at[index, "pose_file"] = str(pose_file)
            manifest.at[index, "score_csv"] = str(score_file)
        except Exception as error:
            manifest.at[index, "pipeline_status"] = "failed"
            manifest.at[index, "pipeline_error"] = f"{type(error).__name__}: {error}"
            failures.append((molecule_id, error))
            print(f"FAILED {molecule_id}: {type(error).__name__}: {error}")
        finally:
            write_pipeline_manifest(manifest, manifest_file)
            write_combined_scores(manifest, combined_file, step_col)

        if failures and fail_fast:
            break

    completed = int((manifest["pipeline_status"] == "completed").sum())
    unavailable = int((manifest["pipeline_status"] == "unavailable").sum())
    failed = int((manifest["pipeline_status"] == "failed").sum())
    write_pipeline_manifest(manifest, manifest_file)
    write_combined_scores(manifest, combined_file, step_col)

    print("\nCross-docking pipeline finished")
    print(f"Completed: {completed}")
    print(f"Failed: {failed}")
    print(f"Unavailable: {unavailable}")
    print(f"Manifest: {manifest_file}")
    if combined_file.is_file():
        print(f"Combined scores: {combined_file}")

    if failures:
        raise RuntimeError(
            f"{len(failures)} docking job(s) failed. See {manifest_file} for details."
        )
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the sequential AHC pose cross-docking pipeline"
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
    parser.add_argument("--limit", type=int, default=None, help="Testing limit")
    parser.add_argument(
        "--precision", default="SP", choices=["HTVS", "SP", "XP"]
    )
    parser.add_argument(
        "--allow-unavailable",
        action="store_true",
        help="Skip source files that are missing or empty",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse completed or partially completed job outputs",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop after the first failed docking job",
    )
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
    )
