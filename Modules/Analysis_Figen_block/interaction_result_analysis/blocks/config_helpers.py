"""Parse named integer and path mappings used by result-only analyses."""

from pathlib import Path


def parse_int_assignments(values: list[str], *, label: str) -> dict[str, int]:
    parsed: dict[str, int] = {}
    for raw in values:
        if "=" not in raw:
            raise ValueError(f"{label} must use NAME=VALUE: {raw!r}")
        key, value = (part.strip() for part in raw.split("=", 1))
        if not key or key in parsed:
            raise ValueError(f"Invalid or duplicate {label} key: {key!r}")
        parsed[key] = int(value)
        if parsed[key] <= 0:
            raise ValueError(f"{label} values must be positive: {raw!r}")
    if not parsed:
        raise ValueError(f"At least one {label} is required")
    return parsed


def parse_path_assignments(values: list[str], *, label: str) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for raw in values:
        if "=" not in raw:
            raise ValueError(f"{label} must use NAME=PATH: {raw!r}")
        key, value = (part.strip() for part in raw.split("=", 1))
        if not key or not value or key in parsed:
            raise ValueError(f"Invalid or duplicate {label}: {raw!r}")
        parsed[key] = Path(value)
    if not parsed:
        raise ValueError(f"At least one {label} is required")
    return parsed
