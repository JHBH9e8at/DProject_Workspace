"""Allocate collision-safe output names without overwriting existing files."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def _validate_filename(filename: str) -> str:
    if not isinstance(filename, str):
        raise TypeError("filename must be a string")
    value = filename.strip()
    if not value or Path(value).name != value or value in {".", ".."}:
        raise ValueError(f"filename must be one non-empty path component: {filename!r}")
    return value


def _with_index(filename: str, index: int) -> str:
    path = Path(filename)
    suffix = "".join(path.suffixes)
    stem = path.name[: -len(suffix)] if suffix else path.name
    version = "" if index == 0 else f"_{index}"
    return f"{stem}{version}{suffix}"


def allocate_versioned_group(
    output_dir: str | Path,
    filenames: Iterable[str],
) -> tuple[Path, ...]:
    """Allocate one shared suffix for a related group of output filenames."""

    directory = Path(output_dir).expanduser().resolve(strict=False)
    directory.mkdir(parents=True, exist_ok=True)
    names = tuple(_validate_filename(filename) for filename in filenames)
    if not names:
        raise ValueError("filenames must contain at least one name")
    if len(set(names)) != len(names):
        raise ValueError("filenames must not contain duplicates")
    index = 0
    while True:
        candidates = tuple(directory / _with_index(name, index) for name in names)
        if not any(candidate.exists() for candidate in candidates):
            return candidates
        index += 1


def next_available_path(output_dir: str | Path, filename: str) -> Path:
    """Return filename or filename_1, filename_2, ... when occupied."""

    return allocate_versioned_group(output_dir, (filename,))[0]

