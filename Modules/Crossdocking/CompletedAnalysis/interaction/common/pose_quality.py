"""Shared serial/parallel PoseCheck adapter and quality thresholding."""

from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import multiprocessing
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd


_POSECHECK_WORKER_CLASS = None
_POSECHECK_WORKER_PROTEIN = None


def _limit_native_threads(threads=1):
    """Prevent each worker process from creating its own large native thread pool."""
    value = str(int(threads))
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = value


def _require_posecheck():
    try:
        from posecheck import PoseCheck
    except ImportError as error:
        raise RuntimeError(
            "PoseCheck is required for pose-quality analysis but is not installed. "
            "Install it in a dedicated compatible environment before running this step."
        ) from error
    return PoseCheck


def _pdb_has_explicit_hydrogens(protein_pdb):
    with Path(protein_pdb).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            atom_name = line[12:16].strip().upper() if len(line) >= 16 else ""
            if element == "H" or (not element and atom_name.startswith("H")):
                return True
    return False


def _load_posecheck_protein(checker, protein_pdb, verbose=True):
    """Preserve prepared receptor hydrogens; use Reduce only when they are absent."""
    if _pdb_has_explicit_hydrogens(protein_pdb):
        from rdkit import Chem

        if verbose:
            print(
                "PoseCheck: explicit receptor hydrogens detected; "
                "loading prepared PDB directly and bypassing Reduce",
                flush=True,
            )
        protein = Chem.MolFromPDBFile(
            str(protein_pdb),
            sanitize=False,
            removeHs=False,
        )
        if protein is None:
            raise ValueError(
                f"RDKit could not load the prepared receptor PDB: {protein_pdb}"
            )
        checker.protein = protein
        if verbose:
            print(
                f"PoseCheck: prepared receptor loaded ({protein.GetNumAtoms()} atoms)",
                flush=True,
            )
        return

    if verbose:
        print(
            "PoseCheck: receptor has no explicit hydrogens; running Reduce",
            flush=True,
        )
    checker.load_protein_from_pdb(str(protein_pdb))


def _raw_posecheck_metrics(protein_pdb, ligand_sdf, verbose=True):
    """Return raw PoseCheck metrics for one SDF using local zero-based indices."""
    PoseCheck = _require_posecheck()
    checker = PoseCheck()
    if verbose:
        print("PoseCheck: loading receptor", flush=True)
    _load_posecheck_protein(checker, protein_pdb, verbose=verbose)
    if verbose:
        print("PoseCheck: loading ligand poses", flush=True)
    checker.load_ligands_from_sdf(str(ligand_sdf))

    if verbose:
        print("PoseCheck: calculating clashes", flush=True)
    clashes = list(checker.calculate_clashes())
    if verbose:
        print("PoseCheck: calculating strain energy", flush=True)
    strain = list(checker.calculate_strain_energy())
    if len(clashes) != len(strain):
        raise RuntimeError(
            "PoseCheck returned different numbers of clash and strain results"
        )
    return pd.DataFrame(
        {
            "pose_index": range(len(clashes)),
            "clash_count": pd.to_numeric(clashes, errors="coerce"),
            "strain_energy": pd.to_numeric(strain, errors="coerce"),
        }
    )


def _initialize_posecheck_worker(protein_pdb, native_threads=1):
    """Load the prepared receptor once into each spawned worker process."""
    global _POSECHECK_WORKER_CLASS, _POSECHECK_WORKER_PROTEIN
    _limit_native_threads(native_threads)
    PoseCheck = _require_posecheck()
    checker = PoseCheck()
    _load_posecheck_protein(checker, Path(protein_pdb).resolve(), verbose=False)
    _POSECHECK_WORKER_CLASS = PoseCheck
    _POSECHECK_WORKER_PROTEIN = checker.protein


def _raw_posecheck_metrics_from_worker(ligand_sdf):
    if _POSECHECK_WORKER_CLASS is None or _POSECHECK_WORKER_PROTEIN is None:
        raise RuntimeError("PoseCheck worker receptor was not initialized")
    checker = _POSECHECK_WORKER_CLASS()
    checker.protein = _POSECHECK_WORKER_PROTEIN
    checker.load_ligands_from_sdf(str(ligand_sdf))
    clashes = list(checker.calculate_clashes())
    strain = list(checker.calculate_strain_energy())
    if len(clashes) != len(strain):
        raise RuntimeError(
            "PoseCheck returned different numbers of clash and strain results"
        )
    return pd.DataFrame(
        {
            "pose_index": range(len(clashes)),
            "clash_count": pd.to_numeric(clashes, errors="coerce"),
            "strain_energy": pd.to_numeric(strain, errors="coerce"),
        }
    )


def _attach_posecheck_status_and_metadata(quality, pose_metadata=None):
    quality = quality.copy()
    quality["posecheck_status"] = np.where(
        np.isfinite(quality["clash_count"])
        & np.isfinite(quality["strain_energy"]),
        "completed",
        "invalid_metric",
    )

    if pose_metadata is not None:
        valid_count = int((pose_metadata["pose_read_status"] == "valid").sum())
        if valid_count != len(quality):
            raise RuntimeError(
                "PoseCheck result count does not match the number of valid SDF poses: "
                f"{len(quality)} versus {valid_count}"
            )
        quality = pose_metadata.merge(
            quality, on="pose_index", how="left", validate="one_to_one"
        )

    if "heavy_atom_count" in quality.columns:
        heavy_atoms = pd.to_numeric(quality["heavy_atom_count"], errors="coerce")
        quality["clashes_per_heavy_atom"] = quality["clash_count"] / heavy_atoms
        quality.loc[heavy_atoms <= 0, "clashes_per_heavy_atom"] = np.nan
    if "rotatable_bond_count" in quality.columns:
        rotors = pd.to_numeric(quality["rotatable_bond_count"], errors="coerce")
        quality["strain_per_rotatable_bond"] = quality["strain_energy"] / rotors
        quality.loc[rotors <= 0, "strain_per_rotatable_bond"] = np.nan
    return quality


def run_posecheck(protein_pdb, ligand_sdf, pose_metadata=None):
    """Calculate raw clash and strain metrics without constructing a composite PQI."""
    protein_pdb = Path(protein_pdb).resolve()
    ligand_sdf = Path(ligand_sdf).resolve()
    if protein_pdb.suffix.lower() != ".pdb":
        raise ValueError("PoseCheck currently requires the protein input as .pdb")

    quality = _raw_posecheck_metrics(protein_pdb, ligand_sdf, verbose=True)
    return _attach_posecheck_status_and_metadata(quality, pose_metadata)


def _write_sdf_chunk(records, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        for record in records:
            handle.write(record)
            if not record.endswith(b"\n"):
                handle.write(b"\n")


def split_sdf_chunks(ligand_sdf, chunk_dir, chunk_size):
    """Split an SDF on record boundaries and return auditable chunk specifications."""
    if chunk_size <= 0:
        raise ValueError("posecheck_chunk_size must be a positive integer")
    ligand_sdf = Path(ligand_sdf).resolve()
    chunk_dir = Path(chunk_dir).resolve()
    chunk_dir.mkdir(parents=True, exist_ok=True)

    for old_chunk in chunk_dir.glob("chunk_*.sdf"):
        old_chunk.unlink()

    specs = []
    records = []
    record = bytearray()
    start_index = 0
    with ligand_sdf.open("rb") as source:
        for line in source:
            record.extend(line)
            if line.strip() != b"$$$$":
                continue
            records.append(bytes(record))
            record.clear()
            if len(records) == chunk_size:
                chunk_index = len(specs)
                path = chunk_dir / f"chunk_{chunk_index:05d}.sdf"
                _write_sdf_chunk(records, path)
                specs.append(
                    {
                        "chunk_index": chunk_index,
                        "start_index": start_index,
                        "pose_count": len(records),
                        "sdf_file": str(path),
                        "chunk_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    }
                )
                start_index += len(records)
                records = []
    if record.strip():
        raise ValueError(f"Unterminated SDF record found in {ligand_sdf}")
    if records:
        chunk_index = len(specs)
        path = chunk_dir / f"chunk_{chunk_index:05d}.sdf"
        _write_sdf_chunk(records, path)
        specs.append(
            {
                "chunk_index": chunk_index,
                "start_index": start_index,
                "pose_count": len(records),
                "sdf_file": str(path),
                "chunk_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    if not specs:
        raise ValueError(f"No SDF records found: {ligand_sdf}")
    return specs


def _run_posecheck_chunk(protein_pdb, chunk_spec, result_file):
    """Process one SDF chunk in a worker process and checkpoint its raw metrics."""
    started = time.monotonic()
    if _POSECHECK_WORKER_PROTEIN is None:
        quality = _raw_posecheck_metrics(
            Path(protein_pdb).resolve(),
            Path(chunk_spec["sdf_file"]).resolve(),
            verbose=False,
        )
    else:
        quality = _raw_posecheck_metrics_from_worker(
            Path(chunk_spec["sdf_file"]).resolve()
        )
    expected = int(chunk_spec["pose_count"])
    if len(quality) != expected:
        raise RuntimeError(
            f"Chunk {chunk_spec['chunk_index']} returned {len(quality)} poses; "
            f"expected {expected}"
        )
    quality["pose_index"] += int(chunk_spec["start_index"])
    quality.insert(0, "chunk_index", int(chunk_spec["chunk_index"]))
    quality["chunk_sha256"] = chunk_spec["chunk_sha256"]
    quality["worker_pid"] = os.getpid()
    quality["chunk_elapsed_seconds"] = time.monotonic() - started
    result_file = Path(result_file)
    result_file.parent.mkdir(parents=True, exist_ok=True)
    temporary = result_file.with_suffix(".csv.tmp")
    quality.to_csv(temporary, index=False)
    temporary.replace(result_file)
    return {
        "chunk_index": int(chunk_spec["chunk_index"]),
        "pose_count": len(quality),
        "result_file": str(result_file),
        "elapsed_seconds": time.monotonic() - started,
        "worker_pid": os.getpid(),
    }


def _chunk_result_is_valid(result_file, chunk_spec):
    result_file = Path(result_file)
    if not result_file.is_file() or result_file.stat().st_size <= 0:
        return False
    try:
        result = pd.read_csv(
            result_file,
            usecols=["pose_index", "chunk_sha256"],
        )
    except (OSError, ValueError, pd.errors.ParserError):
        return False
    expected_indices = list(
        range(
            int(chunk_spec["start_index"]),
            int(chunk_spec["start_index"]) + int(chunk_spec["pose_count"]),
        )
    )
    return (
        len(result) == int(chunk_spec["pose_count"])
        and result["pose_index"].tolist() == expected_indices
        and result["chunk_sha256"].eq(chunk_spec["chunk_sha256"]).all()
    )


def run_posecheck_parallel(
    protein_pdb,
    ligand_sdf,
    pose_metadata,
    work_dir,
    workers,
    chunk_size=500,
    resume=False,
    executor_factory=None,
    worker_function=None,
    native_threads_per_worker=1,
):
    """Run PoseCheck over SDF chunks in independent worker processes."""
    workers = int(workers)
    chunk_size = int(chunk_size)
    if workers <= 0:
        raise ValueError("posecheck_workers must be a positive integer")
    if workers == 1:
        return run_posecheck(protein_pdb, ligand_sdf, pose_metadata)

    available = (
        len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else (os.cpu_count() or 1)
    )
    workers = min(workers, available)
    _limit_native_threads(native_threads_per_worker)
    work_dir = Path(work_dir).resolve()
    chunk_dir = work_dir / "sdf_chunks"
    result_dir = work_dir / "results"
    result_dir.mkdir(parents=True, exist_ok=True)
    specs = split_sdf_chunks(ligand_sdf, chunk_dir, chunk_size)
    expected_total = sum(int(spec["pose_count"]) for spec in specs)
    if expected_total != len(pose_metadata):
        raise RuntimeError(
            f"SDF chunk count {expected_total} does not match pose metadata "
            f"count {len(pose_metadata)}"
        )

    pending = []
    reused = 0
    for spec in specs:
        result_file = result_dir / f"chunk_{spec['chunk_index']:05d}.csv"
        spec["result_file"] = str(result_file)
        if resume and _chunk_result_is_valid(result_file, spec):
            reused += 1
        else:
            pending.append(spec)

    state = (
        str(pose_metadata["receptor_state"].iloc[0])
        if "receptor_state" in pose_metadata.columns and len(pose_metadata)
        else "population"
    )
    print(
        f"PoseCheck MP [{state}]: "
        f"{expected_total} poses, {len(specs)} chunks, {workers} workers, "
        f"{reused} reused",
        flush=True,
    )
    worker_function = worker_function or _run_posecheck_chunk
    completed = reused
    reused_poses = sum(
        int(spec["pose_count"])
        for spec in specs
        if resume and _chunk_result_is_valid(spec["result_file"], spec)
    )
    completed_poses = reused_poses
    run_started = time.monotonic()
    if pending:
        if executor_factory is None:
            context = multiprocessing.get_context("spawn")
            executor = ProcessPoolExecutor(
                max_workers=workers,
                mp_context=context,
                initializer=_initialize_posecheck_worker,
                initargs=(
                    str(Path(protein_pdb).resolve()),
                    native_threads_per_worker,
                ),
            )
        else:
            executor = executor_factory(workers)
        with executor:
            futures = {
                executor.submit(
                    worker_function,
                    str(Path(protein_pdb).resolve()),
                    spec,
                    spec["result_file"],
                ): spec
                for spec in pending
            }
            for future in as_completed(futures):
                spec = futures[future]
                result = future.result()
                completed += 1
                completed_poses += int(result["pose_count"])
                elapsed = max(time.monotonic() - run_started, 1e-9)
                newly_completed = max(completed_poses - reused_poses, 0)
                throughput = newly_completed / elapsed
                remaining = expected_total - completed_poses
                eta_seconds = remaining / throughput if throughput > 0 else float("inf")
                eta_text = (
                    time.strftime("%H:%M:%S", time.gmtime(eta_seconds))
                    if np.isfinite(eta_seconds)
                    else "unknown"
                )
                print(
                    f"PoseCheck MP [{state}]: "
                    f"chunk {result['chunk_index'] + 1}/{len(specs)} completed "
                    f"({result['pose_count']} poses, "
                    f"{result['elapsed_seconds']:.1f}s); "
                    f"{completed_poses}/{expected_total} poses, "
                    f"{throughput:.2f} poses/s, ETA {eta_text}",
                    flush=True,
                )

    tables = []
    for spec in specs:
        result_file = Path(spec["result_file"])
        if not _chunk_result_is_valid(result_file, spec):
            raise RuntimeError(
                f"Missing or invalid PoseCheck chunk result: {result_file}"
            )
        tables.append(pd.read_csv(result_file))
    quality = pd.concat(tables, ignore_index=True).sort_values("pose_index")
    quality = quality.drop(
        columns=[
            "chunk_index",
            "chunk_sha256",
            "worker_pid",
            "chunk_elapsed_seconds",
        ]
    ).reset_index(drop=True)
    if quality["pose_index"].tolist() != list(range(expected_total)):
        raise RuntimeError("Parallel PoseCheck output has missing or duplicate pose indices")
    return _attach_posecheck_status_and_metadata(quality, pose_metadata)


def apply_pose_quality_thresholds(
    quality,
    max_clashes=None,
    max_clashes_per_heavy_atom=None,
    max_strain_energy=None,
):
    """Apply explicit QC thresholds while preserving every raw metric."""
    assessed = quality.copy()
    assessed["pose_quality_pass"] = assessed["posecheck_status"].eq("completed")
    failure_reasons = pd.Series("", index=assessed.index, dtype="object")

    criteria = [
        ("clash_count", max_clashes, "clash_count"),
        (
            "clashes_per_heavy_atom",
            max_clashes_per_heavy_atom,
            "clashes_per_heavy_atom",
        ),
        ("strain_energy", max_strain_energy, "strain_energy"),
    ]
    for column, threshold, label in criteria:
        if threshold is None:
            continue
        if threshold < 0:
            raise ValueError(f"{label} threshold must be nonnegative")
        if column not in assessed.columns:
            raise ValueError(f"Quality table does not contain {column}")
        failed = pd.to_numeric(assessed[column], errors="coerce") > threshold
        assessed.loc[failed, "pose_quality_pass"] = False
        failure_reasons.loc[failed] = failure_reasons.loc[failed].map(
            lambda current: f"{current};{label}" if current else label
        )

    invalid = ~assessed["posecheck_status"].eq("completed")
    failure_reasons.loc[invalid] = "invalid_posecheck_metric"
    assessed["pose_quality_failure_reason"] = failure_reasons
    return assessed
