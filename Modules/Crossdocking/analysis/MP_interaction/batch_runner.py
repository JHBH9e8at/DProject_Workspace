"""Batch interaction analysis for original PPS and PR AHC docking results."""

import json
import math
import re
import gzip
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


STATE_SPECS = {
    "PPS": {
        "results_key": "pps_results",
        "run_dir_key": "pps_run_dir",
        "run_name_key": "pps_run_name",
        "variant_col_key": "pps_variant_col",
        "receptor_key": "pps_receptor",
    },
    "PR": {
        "results_key": "pr_results",
        "run_dir_key": "pr_run_dir",
        "run_name_key": "pr_run_name",
        "variant_col_key": "pr_variant_col",
        "receptor_key": "pr_receptor",
    },
}


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _safe_name(value):
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return name or "molecule"


def _true_mask(series):
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def select_population(
    results_csv,
    run_name,
    step_col="step",
    variant_col=None,
    smiles_col="smiles",
    testmode_top_n=None,
):
    """Select usable AHC results; optionally restrict each state to its top N."""
    run_name = str(run_name).strip()
    variant_col = variant_col or f"{run_name}_best_variant"
    score_col = f"{run_name}_r_i_docking_score"
    table = pd.read_csv(results_csv)
    required = {step_col, variant_col, smiles_col, score_col}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{results_csv}: missing required column(s): {missing}")

    working = table.copy()
    if "valid" in working.columns:
        working = working.loc[_true_mask(working["valid"])].copy()
    if "unique" in working.columns:
        working = working.loc[_true_mask(working["unique"])].copy()

    working[score_col] = pd.to_numeric(working[score_col], errors="coerce")
    finite_score = working[score_col].map(
        lambda value: pd.notna(value) and math.isfinite(value) and value != 0
    )
    usable_variant = (
        working[variant_col].notna()
        & working[variant_col].astype(str).str.strip().ne("")
        & working[variant_col].astype(str).str.strip().ne("0.0")
    )
    usable_smiles = (
        working[smiles_col].notna()
        & working[smiles_col].astype(str).str.strip().ne("")
    )
    working = working.loc[finite_score & usable_variant & usable_smiles].copy()
    if working.empty:
        raise ValueError(f"{results_csv}: no usable docking results remain")

    working[step_col] = pd.to_numeric(working[step_col], errors="raise").astype(int)
    working[variant_col] = working[variant_col].astype(str).str.strip()
    working = working.sort_values(score_col, ascending=True)
    working = working.drop_duplicates(variant_col, keep="first").reset_index(drop=True)
    working.insert(0, "interaction_selection_rank", range(1, len(working) + 1))

    if testmode_top_n is not None:
        if testmode_top_n <= 0:
            raise ValueError("testmode_top_n must be a positive integer")
        working = working.head(testmode_top_n).copy()
    return working


def build_population_manifest(
    results_csv,
    run_dir,
    run_name,
    receptor_file,
    state,
    step_col="step",
    variant_col=None,
    smiles_col="smiles",
    testmode_top_n=None,
):
    """Resolve selected result rows to their original AHC Glide pose files."""
    selected = select_population(
        results_csv=results_csv,
        run_name=run_name,
        step_col=step_col,
        variant_col=variant_col,
        smiles_col=smiles_col,
        testmode_top_n=testmode_top_n,
    )
    variant_col = variant_col or f"{run_name}_best_variant"
    glide_dir = Path(run_dir).resolve() / f"{run_name}_GlideDock"
    receptor_file = Path(receptor_file).resolve()

    records = []
    for _, row in selected.iterrows():
        step = int(row[step_col])
        variant = str(row[variant_col])
        pose_file = glide_dir / str(step) / f"{variant}_lib.sdfgz"
        record = row.to_dict()
        record.update(
            {
                "population_state": state,
                "run_name": run_name,
                "molecule_id": variant,
                "source_step": step,
                "input_smiles": str(row[smiles_col]),
                "receptor_file": str(receptor_file),
                "pose_file": str(pose_file),
                "pose_file_exists": pose_file.is_file(),
                "pose_file_size_bytes": (
                    pose_file.stat().st_size if pose_file.is_file() else 0
                ),
                "interaction_status": "pending",
                "interaction_error": "",
                "job_dir": "",
            }
        )
        records.append(record)
    return pd.DataFrame(records)


def _job_is_complete(job_dir, analysis):
    required = [job_dir / "pose_metadata.csv", job_dir / "run_summary.csv"]
    if analysis in {"posecheck", "both"}:
        required.append(job_dir / "pose_quality.csv")
    if analysis in {"prolif", "both"}:
        required.append(job_dir / "interactions_long.csv")
    return all(path.is_file() for path in required)


def _materialize_population_sdf(manifest, output_sdf, source_map_csv):
    """Concatenate per-molecule SDFGZ files and map global pose indices."""
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    source_rows = []
    global_pose_index = 0
    with output_sdf.open("wb") as destination:
        for _, row in manifest.iterrows():
            if not bool(row["pose_file_exists"]) or int(row["pose_file_size_bytes"]) <= 0:
                continue
            pose_file = Path(row["pose_file"])
            with gzip.open(pose_file, "rb") as source:
                payload = source.read()
            record_count = payload.count(b"$$$$")
            if record_count <= 0:
                raise ValueError(f"No SDF record terminator found: {pose_file}")
            destination.write(payload)
            if not payload.endswith(b"\n"):
                destination.write(b"\n")
            for source_record in range(1, record_count + 1):
                source_rows.append(
                    {
                        "pose_index": global_pose_index,
                        "population_state": row["population_state"],
                        "run_name": row["run_name"],
                        "molecule_id": row["molecule_id"],
                        "source_step": row["source_step"],
                        "input_smiles": row["input_smiles"],
                        "source_pose_file": row["pose_file"],
                        "source_sdf_record": source_record,
                    }
                )
                global_pose_index += 1
    source_map = pd.DataFrame(source_rows)
    if source_map.empty:
        raise ValueError("No available SDF records were materialized")
    source_map.to_csv(source_map_csv, index=False)
    return source_map


def _collect_population_outputs(job_dir, source_map, state, aggregate):
    output_names = {
        "pose_metadata": "pose_metadata.csv",
        "pose_quality": "pose_quality.csv",
        "interactions_long": "interactions_long.csv",
        "run_summary": "run_summary.csv",
    }
    for key, filename in output_names.items():
        path = job_dir / filename
        if path.is_file():
            table = pd.read_csv(path)
            if "pose_index" in table.columns:
                table = table.merge(
                    source_map,
                    on="pose_index",
                    how="left",
                    validate="many_to_one",
                )
            else:
                table.insert(0, "population_state", state)
            aggregate[key].append(table)


def run_batch_interactions(
    pps_results,
    pr_results,
    pps_run_dir,
    pr_run_dir,
    pps_receptor,
    pr_receptor,
    output_dir,
    pps_run_name="PPS",
    pr_run_name="PR",
    step_col="step",
    pps_variant_col="PPS_best_variant",
    pr_variant_col="PR_best_variant",
    smiles_col="smiles",
    testmode_top_n=None,
    analysis="both",
    residue_map=None,
    include_secondary_interactions=False,
    max_clashes=None,
    max_clashes_per_heavy_atom=None,
    max_strain_energy=None,
    posecheck_workers=1,
    posecheck_chunk_size=500,
    prolif_workers=None,
    allow_unavailable=False,
    resume=False,
    fail_fast=False,
    analyzer=None,
):
    """Run the existing single-complex analysis over PPS and PR populations."""
    if int(posecheck_workers) > 1 or (
        prolif_workers is not None and int(prolif_workers) > 1
    ):
        for variable in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            os.environ[variable] = "1"

    if analyzer is None:
        try:
            from .single_complex import run_single_complex
        except ImportError:
            from single_complex import run_single_complex

        analyzer = run_single_complex

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir = output_dir / "manifests"
    jobs_dir = output_dir / "populations"
    aggregate_dir = output_dir / "aggregate"
    for directory in (manifests_dir, jobs_dir, aggregate_dir):
        directory.mkdir(parents=True, exist_ok=True)

    config = {
        "pps_results": pps_results,
        "pr_results": pr_results,
        "pps_run_dir": pps_run_dir,
        "pr_run_dir": pr_run_dir,
        "pps_receptor": pps_receptor,
        "pr_receptor": pr_receptor,
        "pps_run_name": pps_run_name,
        "pr_run_name": pr_run_name,
        "pps_variant_col": pps_variant_col,
        "pr_variant_col": pr_variant_col,
    }
    manifests = []
    for state, spec in STATE_SPECS.items():
        manifest = build_population_manifest(
            results_csv=config[spec["results_key"]],
            run_dir=config[spec["run_dir_key"]],
            run_name=config[spec["run_name_key"]],
            receptor_file=config[spec["receptor_key"]],
            state=state,
            step_col=step_col,
            variant_col=config[spec["variant_col_key"]],
            smiles_col=smiles_col,
            testmode_top_n=testmode_top_n,
        )
        manifests.append(manifest)

    manifest = pd.concat(manifests, ignore_index=True)
    manifest_file = manifests_dir / "interaction_manifest.csv"
    manifest.to_csv(manifest_file, index=False)
    unavailable = ~manifest["pose_file_exists"] | manifest["pose_file_size_bytes"].eq(0)
    if unavailable.any() and not allow_unavailable:
        raise FileNotFoundError(
            f"{int(unavailable.sum())} pose file(s) are missing or empty; "
            f"see {manifest_file}"
        )

    aggregate = {
        "pose_metadata": [],
        "pose_quality": [],
        "interactions_long": [],
        "run_summary": [],
    }
    started_at = _utc_now()
    for state in ("PPS", "PR"):
        state_mask = manifest["population_state"].eq(state)
        ready_mask = (
            state_mask
            & manifest["pose_file_exists"]
            & manifest["pose_file_size_bytes"].gt(0)
        )
        unavailable_mask = state_mask & ~ready_mask
        manifest.loc[unavailable_mask, "interaction_status"] = "unavailable"
        state_manifest = manifest.loc[ready_mask].copy()
        if state_manifest.empty:
            continue

        job_dir = jobs_dir / state
        combined_sdf = job_dir / f"{state}_selected_poses.sdf"
        source_map_csv = job_dir / "pose_source_map.csv"
        manifest.loc[state_mask, "job_dir"] = str(job_dir)
        try:
            source_map = _materialize_population_sdf(
                state_manifest, combined_sdf, source_map_csv
            )
            print(
                f"[{state}] Materialized {len(source_map)} pose record(s): "
                f"{combined_sdf}",
                flush=True,
            )
            if resume and _job_is_complete(job_dir, analysis):
                status = "reused"
                print(f"[{state}] Reusing completed analysis", flush=True)
            else:
                print(
                    f"[{state}] Starting {analysis} analysis with receptor "
                    f"{state_manifest.iloc[0]['receptor_file']}",
                    flush=True,
                )
                analyzer(
                    protein_file=state_manifest.iloc[0]["receptor_file"],
                    pose_file=combined_sdf,
                    receptor_state=state,
                    output_dir=job_dir,
                    analysis=analysis,
                    residue_map=residue_map,
                    include_secondary_interactions=include_secondary_interactions,
                    max_clashes=max_clashes,
                    max_clashes_per_heavy_atom=max_clashes_per_heavy_atom,
                    max_strain_energy=max_strain_energy,
                    posecheck_workers=posecheck_workers,
                    posecheck_chunk_size=posecheck_chunk_size,
                    prolif_workers=prolif_workers,
                    resume=resume,
                )
                status = "completed"
                print(f"[{state}] Analysis completed", flush=True)
            manifest.loc[ready_mask, "interaction_status"] = status
            _collect_population_outputs(
                job_dir, source_map, state, aggregate
            )
        except Exception as error:
            manifest.loc[ready_mask, "interaction_status"] = "failed"
            manifest.loc[ready_mask, "interaction_error"] = str(error)
            manifest.to_csv(manifest_file, index=False)
            if fail_fast:
                raise
        manifest.to_csv(manifest_file, index=False)

    output_files = {}
    for key, tables in aggregate.items():
        if tables:
            output_file = aggregate_dir / f"{key}.csv"
            pd.concat(tables, ignore_index=True).to_csv(output_file, index=False)
            output_files[key] = str(output_file)

    counts = manifest["interaction_status"].value_counts().to_dict()
    failed_count = int(counts.get("failed", 0))
    summary = {
        "status": "failed" if failed_count else "completed",
        "started_at_utc": started_at,
        "finished_at_utc": _utc_now(),
        "selection_mode": "testmode_top_n" if testmode_top_n else "full_run",
        "testmode_top_n": testmode_top_n,
        "analysis": analysis,
        "posecheck_workers": posecheck_workers,
        "posecheck_chunk_size": posecheck_chunk_size,
        "prolif_workers": prolif_workers,
        "selected_molecule_count": len(manifest),
        "status_counts": counts,
        "manifest_file": str(manifest_file),
        "aggregate_outputs": output_files,
    }
    summary_file = output_dir / "interaction_run_summary.json"
    summary_file.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    if failed_count:
        raise RuntimeError(
            f"{failed_count} interaction job(s) failed; see {manifest_file}"
        )
    return manifest, summary
