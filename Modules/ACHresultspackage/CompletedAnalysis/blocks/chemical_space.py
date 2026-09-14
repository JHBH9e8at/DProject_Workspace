"""Run or reuse joint chemical-space calculation, then render matching figures."""

import json
from pathlib import Path
from Modules.Analysis_Figen_block.ahc_analysis.blocks.chemical_space_runner import run_chemical_space
from Modules.Analysis_Figen_block.common.provenance import sha256_file
from Modules.Analysis_Figen_block.figures.blocks.chemical_space_workflow import run_chemical_space_figures


def _find_cache(*, results_root: Path, inputs: tuple[Path, ...],
                methods: tuple[str, ...], features: tuple[str, ...],
                random_state: int):
    """Return the newest verified cache whose inputs and settings are identical."""
    cache_root = results_root / "analysis" / "chemical_space"
    if not cache_root.is_dir():
        return None
    input_hashes = [sha256_file(path) for path in inputs]
    manifests = sorted(
        cache_root.glob("*/run_manifest.json"),
        key=lambda path: path.stat().st_mtime_ns,
        reverse=True,
    )
    expected_config = {
        "methods": list(methods), "features": list(features),
        "random_state": random_state,
    }
    for manifest_path in manifests:
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest.get("script_id") != "PY-137" or manifest.get("status") != "completed":
                continue
            config = manifest.get("config", {})
            if any(config.get(key) != value for key, value in expected_config.items()):
                continue
            records = manifest.get("inputs", [])
            if [record.get("sha256") for record in records] != input_hashes:
                continue
            run_dir = manifest_path.parent
            outputs = tuple(run_dir / record["path"] for record in manifest.get("outputs", []))
            if not outputs or any(
                not path.is_file() or sha256_file(path) != record.get("sha256")
                for path, record in zip(outputs, manifest.get("outputs", []))
            ):
                continue
            return manifest, outputs
        except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
            continue
    return None


def run(*, pr_input: Path, pps_input: Path, reference: Path,
        results_root: Path, run_id: str, resume: bool,
        methods: tuple[str, ...], features: tuple[str, ...], random_state: int,
        render_figures: bool = True, usecache: bool = True):
    inputs = (pr_input, pps_input, reference)
    cached = _find_cache(
        results_root=results_root, inputs=inputs, methods=methods,
        features=features, random_state=random_state,
    ) if usecache else None
    if cached:
        manifest, outputs = cached
        context = None
        upstream_run_id = str(manifest["run"]["run_id"])
        print(f"[cache hit] {manifest['run']['run_dir']}")
    else:
        context, outputs = run_chemical_space(
            pr_path=pr_input, pps_path=pps_input, reference_path=reference,
            methods=methods, features=features, random_state=random_state,
            run_id=f"{run_id}_chemical_space", results_root=results_root, resume=resume,
        )
        upstream_run_id = context.run_id
    membership = next(
        (path for path in outputs
         if path.name.startswith("joint_umap_fp") and path.name.endswith(".csv.gz")),
        None,
    )
    figures = []
    if render_figures:
        for cache in outputs:
            if not cache.name.endswith(".csv.gz"):
                continue
            figures.append(run_chemical_space_figures(
                cache_path=cache, run_id=f"{run_id}_{cache.name.removesuffix('.csv.gz')}",
                results_root=results_root, resume=resume,
                upstream_run_ids=(upstream_run_id,),
            ))
    return {
        "context": context,
        "outputs": outputs,
        "figures": figures,
        "membership": membership,
        "cache_hit": cached is not None,
        "cache_run_id": upstream_run_id,
    }
