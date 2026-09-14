"""Run ligand-local Fpocket volume analysis in an isolated work_results run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import RunContext, create_run_context, generate_run_id
from .pocket_volume import calculate_pocket_volumes, parse_targets


SCRIPT_ID = "PY-093"
CSV_FIELDS = [
    "scope",
    "pocket_file",
    "total_alpha_spheres",
    "selected_alpha_spheres",
    "union_volume_A3",
    "standard_error_A3",
    "naive_sphere_volume_sum_A3",
    "bounding_box_volume_A3",
    "points_inside",
    "iterations",
    "bounding_box_min",
    "bounding_box_max",
]


def _write_outputs(
    rows: list[dict[str, object]],
    payload: dict[str, object],
    csv_path: Path,
    json_path: Path,
) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_pocket_volume(
    *,
    fpocket_dir: str | Path,
    targets: list[int],
    ligand: str | Path,
    run_id: str,
    results_root: str | Path | None = None,
    resume: bool = False,
    output_prefix: str = "ligand_local_pocket_volume",
    ligand_resname: str | None = None,
    distance_threshold: float = 3.0,
    all_spheres: bool = False,
    iterations: int = 1_000_000,
    seed: int = 42,
    chunk_size: int = 100_000,
    radius_offset: float = -1.6,
    upstream_run_ids=(),
) -> tuple[RunContext, Path, Path]:
    """Calculate, write, and record one pocket-volume analysis run."""

    context = create_run_context(
        ("analysis", "pocket_volume"),
        run_id,
        results_root=results_root,
        resume=resume,
    )
    config = {
        "fpocket_dir": str(Path(fpocket_dir).expanduser().resolve(strict=False)),
        "targets": list(targets),
        "ligand": str(Path(ligand).expanduser().resolve(strict=False)),
        "ligand_resname": ligand_resname,
        "distance_threshold": distance_threshold,
        "all_spheres": all_spheres,
        "iterations": iterations,
        "seed": seed,
        "chunk_size": chunk_size,
        "radius_offset": radius_offset,
        "output_prefix": output_prefix,
    }
    try:
        rows, payload, pocket_files = calculate_pocket_volumes(
            fpocket_dir=fpocket_dir,
            targets=targets,
            ligand=ligand,
            ligand_resname=ligand_resname,
            distance_threshold=distance_threshold,
            all_spheres=all_spheres,
            iterations=iterations,
            seed=seed,
            chunk_size=chunk_size,
            radius_offset=radius_offset,
        )
        csv_path, json_path = allocate_versioned_group(
            context.tables_dir,
            (f"{output_prefix}.csv", f"{output_prefix}.json"),
        )
        _write_outputs(rows, payload, csv_path, json_path)
        update_manifest(
            context,
            script_id=SCRIPT_ID,
            status="completed",
            config=config,
            inputs=(Path(ligand), *pocket_files),
            outputs=(csv_path, json_path),
            upstream_run_ids=upstream_run_ids,
        )
    except Exception:
        update_manifest(
            context,
            script_id=SCRIPT_ID,
            status="failed",
            config=config,
            upstream_run_ids=upstream_run_ids,
        )
        raise
    return context, csv_path, json_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate ligand-local Fpocket alpha-sphere union volume and "
            "write an isolated work_results run."
        )
    )
    parser.add_argument("--fpocket-dir", required=True, type=Path)
    parser.add_argument("--target", required=True, type=parse_targets)
    parser.add_argument("--ligand", required=True, type=Path)
    parser.add_argument("--ligand-resname", type=str.upper)
    parser.add_argument("--distance-threshold", type=float, default=3.0)
    parser.add_argument("--all-spheres", action="store_true")
    parser.add_argument("--iterations", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--radius-offset", type=float, default=-1.6)
    parser.add_argument("--output-prefix", default="ligand_local_pocket_volume")
    parser.add_argument("--target-label", default="pocket")
    parser.add_argument("--run-id")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--results-root", type=Path); parser.add_argument("--upstream-run-id", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_id = args.run_id or generate_run_id(args.target_label, "pocket_volume")
    context, csv_path, json_path = run_pocket_volume(
        fpocket_dir=args.fpocket_dir,
        targets=args.target,
        ligand=args.ligand,
        run_id=run_id,
        results_root=validate_results_root(args.results_root) if args.results_root else None,
        resume=args.resume,
        output_prefix=args.output_prefix,
        ligand_resname=args.ligand_resname,
        distance_threshold=args.distance_threshold,
        all_spheres=args.all_spheres,
        iterations=args.iterations,
        seed=args.seed,
        chunk_size=args.chunk_size,
        radius_offset=args.radius_offset,
        upstream_run_ids=tuple(args.upstream_run_id),
    )
    print(f"Run directory: {context.run_dir}")
    print(f"CSV:  {csv_path}")
    print(f"JSON: {json_path}")
    print(f"Manifest: {context.manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
