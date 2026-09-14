#!/usr/bin/env python3
"""Calculate ligand-local Fpocket alpha-sphere union volumes.

The script selects alpha spheres whose centres lie within a user-specified
distance of any ligand heavy atom, then estimates their union volume using
the same bounding-box Monte Carlo principle used by Fpocket.
"""

from __future__ import annotations

import argparse
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


def calculate_pocket_volumes(
    *,
    fpocket_dir: str | Path,
    targets: list[int],
    ligand: str | Path,
    ligand_resname: str | None = None,
    distance_threshold: float = 3.0,
    all_spheres: bool = False,
    iterations: int = 1_000_000,
    seed: int = 42,
    chunk_size: int = 100_000,
    radius_offset: float = -1.6,
) -> tuple[list[dict[str, object]], dict[str, object], list[Path]]:
    """Calculate per-pocket and combined volumes without performing file output."""

    fpocket_path = Path(fpocket_dir).expanduser().resolve(strict=False)
    ligand_path = Path(ligand).expanduser().resolve(strict=False)
    if not targets or len(set(targets)) != len(targets):
        raise ValueError("targets must contain unique pocket numbers")
    if any(not isinstance(target, int) or target < 0 for target in targets):
        raise ValueError("targets must contain non-negative integers")
    if not all_spheres and distance_threshold <= 0:
        raise ValueError("distance_threshold must be positive")
    if iterations <= 0 or chunk_size <= 0:
        raise ValueError("iterations and chunk_size must be positive")
    if not fpocket_path.is_dir():
        raise NotADirectoryError(fpocket_path)
    if not ligand_path.is_file():
        raise FileNotFoundError(ligand_path)

    normalized_resname = ligand_resname.upper() if ligand_resname else None
    ligand_coordinates = read_ligand_heavy_atoms(ligand_path, normalized_resname)
    rows: list[dict[str, object]] = []
    all_selected: list[Sphere] = []
    pocket_files: list[Path] = []

    for offset, target in enumerate(targets):
        pocket_file = resolve_pocket_file(fpocket_path, target)
        pocket_files.append(pocket_file.resolve())
        available_spheres = read_pocket_spheres(pocket_file)
        selected = (
            available_spheres
            if all_spheres
            else select_ligand_local_spheres(
                available_spheres, ligand_coordinates, distance_threshold
            )
        )
        all_selected.extend(selected)
        estimate = monte_carlo_union_volume(
            selected, iterations, seed + offset, chunk_size, radius_offset
        )
        rows.append(
            {
                "scope": f"pocket{target}",
                "pocket_file": str(pocket_file.resolve()),
                "total_alpha_spheres": len(available_spheres),
                "selected_alpha_spheres": len(selected),
                **estimate,
            }
        )

    combined = unique_spheres(all_selected)
    combined_estimate = monte_carlo_union_volume(
        combined,
        iterations,
        seed + len(targets),
        chunk_size,
        radius_offset,
    )
    rows.append(
        {
            "scope": "combined",
            "pocket_file": ";".join(f"pocket{target}" for target in targets),
            "total_alpha_spheres": "",
            "selected_alpha_spheres": len(combined),
            **combined_estimate,
        }
    )

    payload = {
        "method": "bounding-box Monte Carlo union of selected Fpocket alpha spheres",
        "selection_rule": (
            "all alpha spheres in each target pocket"
            if all_spheres
            else "alpha-sphere centre within threshold of any ligand heavy atom"
        ),
        "fpocket_dir": str(fpocket_path),
        "ligand": str(ligand_path),
        "ligand_resname": normalized_resname,
        "ligand_heavy_atoms": len(ligand_coordinates),
        "targets": targets,
        "distance_threshold_A": None if all_spheres else distance_threshold,
        "iterations": iterations,
        "seed": seed,
        "alpha_sphere_radius_offset_A": radius_offset,
        "results": rows,
    }
    return rows, payload, pocket_files
