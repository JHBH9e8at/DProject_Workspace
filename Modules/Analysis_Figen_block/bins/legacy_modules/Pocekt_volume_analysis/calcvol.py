#!/usr/bin/env python3
"""Calculate ligand-local Fpocket alpha-sphere union volumes.

The script selects alpha spheres whose centres lie within a user-specified
distance of any ligand heavy atom, then estimates their union volume using
the same bounding-box Monte Carlo principle used by Fpocket.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Sphere:
    centre: tuple[float, float, float]
    radius: float


def parse_targets(value: str) -> list[int]:
    tokens = [token.strip() for token in value.split(",") if token.strip()]
    if not tokens or any(not re.fullmatch(r"\d+", token) for token in tokens):
        raise argparse.ArgumentTypeError("--target must look like 2,38,51")
    targets = [int(token) for token in tokens]
    if len(set(targets)) != len(targets):
        raise argparse.ArgumentTypeError("--target contains duplicate pocket numbers")
    return targets


def read_pocket_spheres(path: Path) -> list[Sphere]:
    spheres: list[Sphere] = []
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        fields = line.split()
        try:
            # Fpocket pocket*_vert.pqr: ... x y z charge radius
            x, y, z = map(float, fields[5:8])
            radius = float(fields[-1])
        except (ValueError, IndexError) as exc:
            raise ValueError(f"Could not parse alpha sphere at {path}:{line_number}") from exc
        if radius <= 0:
            raise ValueError(f"Non-positive alpha-sphere radius at {path}:{line_number}")
        spheres.append(Sphere((x, y, z), radius))
    if not spheres:
        raise ValueError(f"No alpha spheres found in {path}")
    return spheres


def infer_element(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    atom_name = line[12:16].strip()
    atom_name = re.sub(r"^[0-9]+", "", atom_name)
    return atom_name[:1].upper()


def read_ligand_heavy_atoms(path: Path, residue_name: str | None = None) -> np.ndarray:
    coordinates: list[tuple[float, float, float]] = []
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if residue_name is not None and line[17:20].strip().upper() != residue_name:
            continue
        if infer_element(line) == "H":
            continue
        try:
            coordinates.append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            )
        except ValueError as exc:
            raise ValueError(f"Could not parse ligand coordinate at {path}:{line_number}") from exc
    if not coordinates:
        raise ValueError(f"No ligand heavy atoms found in {path}")
    return np.asarray(coordinates, dtype=float)


def select_ligand_local_spheres(
    spheres: list[Sphere], ligand_coordinates: np.ndarray, threshold: float
) -> list[Sphere]:
    centres = np.asarray([sphere.centre for sphere in spheres], dtype=float)
    # Centre-to-heavy-atom distance, matching Fpocket ligand-overlap criteria.
    distance_squared = np.sum(
        (centres[:, None, :] - ligand_coordinates[None, :, :]) ** 2, axis=2
    )
    selected = np.any(distance_squared <= threshold**2, axis=1)
    return [sphere for sphere, keep in zip(spheres, selected) if keep]


def monte_carlo_union_volume(
    spheres: list[Sphere], iterations: int, seed: int, chunk_size: int, radius_offset: float
) -> dict[str, float | int | list[float]]:
    if not spheres:
        return {
            "union_volume_A3": 0.0,
            "standard_error_A3": 0.0,
            "bounding_box_volume_A3": 0.0,
            "points_inside": 0,
            "iterations": iterations,
            "naive_sphere_volume_sum_A3": 0.0,
            "bounding_box_min": [],
            "bounding_box_max": [],
        }

    centres = np.asarray([sphere.centre for sphere in spheres], dtype=float)
    raw_radii = np.asarray([sphere.radius for sphere in spheres], dtype=float)
    radii = raw_radii + radius_offset
    if np.any(radii <= 0):
        raise ValueError(
            "The radius offset produces a non-positive effective alpha-sphere radius"
        )
    lower = np.min(centres - radii[:, None], axis=0)
    upper = np.max(centres + radii[:, None], axis=0)
    box_volume = float(np.prod(upper - lower))

    rng = np.random.default_rng(seed)
    points_inside = 0
    completed = 0
    radius_squared = radii**2
    while completed < iterations:
        count = min(chunk_size, iterations - completed)
        points = rng.uniform(lower, upper, size=(count, 3))
        inside_any = np.zeros(count, dtype=bool)
        # Sphere-wise evaluation avoids a potentially very large point-by-sphere array.
        for centre, r2 in zip(centres, radius_squared):
            active = ~inside_any
            if not np.any(active):
                break
            delta = points[active] - centre
            inside_any[active] = np.einsum("ij,ij->i", delta, delta) <= r2
        points_inside += int(np.count_nonzero(inside_any))
        completed += count

    fraction = points_inside / iterations
    volume = box_volume * fraction
    standard_error = box_volume * math.sqrt(fraction * (1.0 - fraction) / iterations)
    naive_sum = float(np.sum((4.0 / 3.0) * math.pi * radii**3))
    return {
        "union_volume_A3": volume,
        "standard_error_A3": standard_error,
        "bounding_box_volume_A3": box_volume,
        "points_inside": points_inside,
        "iterations": iterations,
        "naive_sphere_volume_sum_A3": naive_sum,
        "bounding_box_min": lower.tolist(),
        "bounding_box_max": upper.tolist(),
    }


def unique_spheres(spheres: list[Sphere]) -> list[Sphere]:
    # Identical vertices can occur if multiple Fpocket pockets are combined.
    return list(dict.fromkeys(spheres))


def resolve_pocket_file(fpocket_dir: Path, target: int) -> Path:
    candidates = [
        fpocket_dir / "pockets" / f"pocket{target}_vert.pqr",
        fpocket_dir / f"pocket{target}_vert.pqr",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"Could not find pocket{target}_vert.pqr under {fpocket_dir} or its pockets directory"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate ligand-local alpha-sphere union volume from Fpocket results."
    )
    parser.add_argument("--fpocket-dir", required=True, type=Path)
    parser.add_argument("--target", required=True, type=parse_targets, help="e.g. 2,38,51")
    parser.add_argument("--ligand", required=True, type=Path, help="Ligand PDB file")
    parser.add_argument(
        "--ligand-resname",
        type=str.upper,
        help="Optional PDB residue name when --ligand is a protein-ligand complex (e.g. 2OW)",
    )
    parser.add_argument("--distance-threshold", type=float, default=3.0, metavar="ANGSTROM")
    parser.add_argument(
        "--all-spheres",
        action="store_true",
        help="Use every alpha sphere in each target pocket (no ligand-distance threshold)",
    )
    parser.add_argument("--iterations", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument(
        "--radius-offset",
        type=float,
        default=-1.6,
        metavar="ANGSTROM",
        help="Offset applied to alpha-sphere radii for volume integration (Fpocket default: -1.6)",
    )
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument("--output-prefix", default="ligand_local_pocket_volume")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not args.all_spheres and args.distance_threshold <= 0:
        raise ValueError("--distance-threshold must be positive")
    if args.iterations <= 0 or args.chunk_size <= 0:
        raise ValueError("--iterations and --chunk-size must be positive")
    if not args.fpocket_dir.is_dir():
        raise NotADirectoryError(args.fpocket_dir)
    if not args.ligand.is_file():
        raise FileNotFoundError(args.ligand)

    ligand_coordinates = read_ligand_heavy_atoms(args.ligand, args.ligand_resname)
    rows: list[dict[str, object]] = []
    all_selected: list[Sphere] = []

    for offset, target in enumerate(args.target):
        pocket_file = resolve_pocket_file(args.fpocket_dir, target)
        all_spheres = read_pocket_spheres(pocket_file)
        selected = (
            all_spheres
            if args.all_spheres
            else select_ligand_local_spheres(
                all_spheres, ligand_coordinates, args.distance_threshold
            )
        )
        all_selected.extend(selected)
        estimate = monte_carlo_union_volume(
            selected, args.iterations, args.seed + offset, args.chunk_size, args.radius_offset
        )
        rows.append(
            {
                "scope": f"pocket{target}",
                "pocket_file": str(pocket_file.resolve()),
                "total_alpha_spheres": len(all_spheres),
                "selected_alpha_spheres": len(selected),
                **estimate,
            }
        )

    combined = unique_spheres(all_selected)
    combined_estimate = monte_carlo_union_volume(
        combined,
        args.iterations,
        args.seed + len(args.target),
        args.chunk_size,
        args.radius_offset,
    )
    rows.append(
        {
            "scope": "combined",
            "pocket_file": ";".join(f"pocket{target}" for target in args.target),
            "total_alpha_spheres": "",
            "selected_alpha_spheres": len(combined),
            **combined_estimate,
        }
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"{args.output_prefix}.csv"
    json_path = args.output_dir / f"{args.output_prefix}.json"
    csv_fields = [
        "scope", "pocket_file", "total_alpha_spheres", "selected_alpha_spheres",
        "union_volume_A3", "standard_error_A3", "naive_sphere_volume_sum_A3",
        "bounding_box_volume_A3", "points_inside", "iterations",
        "bounding_box_min", "bounding_box_max",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_fields)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "method": "bounding-box Monte Carlo union of selected Fpocket alpha spheres",
        "selection_rule": (
            "all alpha spheres in each target pocket"
            if args.all_spheres
            else "alpha-sphere centre within threshold of any ligand heavy atom"
        ),
        "fpocket_dir": str(args.fpocket_dir.resolve()),
        "ligand": str(args.ligand.resolve()),
        "ligand_resname": args.ligand_resname,
        "ligand_heavy_atoms": len(ligand_coordinates),
        "targets": args.target,
        "distance_threshold_A": None if args.all_spheres else args.distance_threshold,
        "iterations": args.iterations,
        "seed": args.seed,
        "alpha_sphere_radius_offset_A": args.radius_offset,
        "results": rows,
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"Ligand heavy atoms: {len(ligand_coordinates)}")
    if args.all_spheres:
        print("Selection threshold: none (all target-pocket alpha spheres)")
    else:
        print(f"Selection threshold: {args.distance_threshold:g} A (centre-to-heavy-atom)")
    print(f"Alpha-sphere radius offset for integration: {args.radius_offset:g} A")
    for row in rows:
        print(
            f"{row['scope']}: selected={row['selected_alpha_spheres']}, "
            f"union_volume={float(row['union_volume_A3']):.3f} +/- "
            f"{float(row['standard_error_A3']):.3f} A^3"
        )
    print(f"CSV:  {csv_path.resolve()}")
    print(f"JSON: {json_path.resolve()}")


if __name__ == "__main__":
    main()
