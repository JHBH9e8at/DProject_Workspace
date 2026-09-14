"""Resolve project paths and enforce the work-results safety boundary."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable


WORKING_DIRECTORY_NAME = "work_space"
RESULTS_DIRECTORY_NAME = "work_results"


def _as_directory(path: str | os.PathLike[str] | Path) -> Path:
    candidate = Path(path).expanduser().resolve(strict=False)
    return candidate.parent if candidate.is_file() else candidate


def discover_project_root(
    start: str | os.PathLike[str] | Path | None = None,
) -> Path:
    """Find the directory containing the sibling work_space/work_results roots."""

    current = _as_directory(start or __file__)
    for candidate in (current, *current.parents):
        if candidate.name == WORKING_DIRECTORY_NAME:
            return candidate.parent
        if (candidate / WORKING_DIRECTORY_NAME).is_dir():
            return candidate
        if (
            (candidate / "Modules" / "Analysis_Figen_block").is_dir()
            and (candidate / "requirements").is_dir()
        ):
            return candidate
    raise FileNotFoundError(
        f"Could not find a {WORKING_DIRECTORY_NAME!r} ancestor from {current}"
    )


def get_working_root(
    project_root: str | os.PathLike[str] | Path | None = None,
) -> Path:
    root = (
        Path(project_root).expanduser().resolve(strict=False)
        if project_root is not None
        else discover_project_root()
    )
    working_root = root / WORKING_DIRECTORY_NAME
    if working_root.is_dir():
        return working_root
    if (root / "Modules" / "Analysis_Figen_block").is_dir():
        return root
    raise FileNotFoundError(f"Working root does not exist: {working_root}")


def get_results_root(
    project_root: str | os.PathLike[str] | Path | None = None,
    *,
    create: bool = False,
) -> Path:
    root = (
        Path(project_root).expanduser().resolve(strict=False)
        if project_root is not None
        else discover_project_root()
    )
    results_root = root / RESULTS_DIRECTORY_NAME
    if create:
        results_root.mkdir(parents=True, exist_ok=True)
    return results_root


def validate_results_root(
    path: str | os.PathLike[str] | Path,
) -> Path:
    """Require an explicitly supplied output boundary to be `work_results`."""

    boundary = Path(path).expanduser().resolve(strict=False)
    if boundary.name.casefold() != RESULTS_DIRECTORY_NAME.casefold():
        raise ValueError(
            f"results_root must name the {RESULTS_DIRECTORY_NAME!r} directory: "
            f"{boundary}"
        )
    return boundary


def resolve_config_path(
    value: str | os.PathLike[str] | Path,
    config_path: str | os.PathLike[str] | Path,
) -> Path:
    """Resolve a path relative to the directory containing its config file."""

    candidate = Path(value).expanduser()
    if candidate.is_absolute():
        return candidate.resolve(strict=False)
    base = Path(config_path).expanduser().resolve(strict=False).parent
    return (base / candidate).resolve(strict=False)


def assert_within(
    candidate: str | os.PathLike[str] | Path,
    boundary: str | os.PathLike[str] | Path,
) -> Path:
    """Return a resolved path only when it is inside the supplied boundary."""

    resolved_candidate = Path(candidate).expanduser().resolve(strict=False)
    resolved_boundary = Path(boundary).expanduser().resolve(strict=False)
    try:
        common = Path(os.path.commonpath((resolved_candidate, resolved_boundary)))
    except ValueError as exc:
        raise ValueError(
            f"Path is on a different filesystem than boundary: {resolved_candidate}"
        ) from exc
    if common != resolved_boundary:
        raise ValueError(
            f"Path escapes allowed boundary {resolved_boundary}: {resolved_candidate}"
        )
    return resolved_candidate


def validate_component(value: str, *, label: str) -> str:
    """Validate one portable directory-name component."""

    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    cleaned = value.strip()
    if not cleaned or cleaned in {".", ".."}:
        raise ValueError(f"{label} must not be empty, '.' or '..'")
    if Path(cleaned).name != cleaned or "/" in cleaned or "\\" in cleaned:
        raise ValueError(f"{label} must be one path component: {value!r}")
    if any(ord(character) < 32 for character in cleaned):
        raise ValueError(f"{label} contains a control character")
    return cleaned


def normalize_feature_parts(feature: str | Iterable[str]) -> tuple[str, ...]:
    """Normalize a feature name or explicit feature path components."""

    raw_parts = (feature,) if isinstance(feature, str) else tuple(feature)
    if not raw_parts:
        raise ValueError("feature must contain at least one component")
    return tuple(
        validate_component(part, label=f"feature component {index}")
        for index, part in enumerate(raw_parts, start=1)
    )


def result_run_path(
    results_root: str | os.PathLike[str] | Path,
    feature: str | Iterable[str],
    run_id: str,
) -> Path:
    """Build and boundary-check work_results/<feature>/<run_id>."""

    boundary = Path(results_root).expanduser().resolve(strict=False)
    parts = normalize_feature_parts(feature)
    safe_run_id = validate_component(run_id, label="run_id")
    return assert_within(boundary.joinpath(*parts, safe_run_id), boundary)
