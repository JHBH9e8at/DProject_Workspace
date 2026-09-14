"""Strict helpers for the repository's ``key=value`` configuration files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Collection, Mapping


@dataclass(frozen=True)
class ConfigDocument:
    """Parsed text values together with the file that defines their path base."""

    path: Path
    values: dict[str, str]

    @property
    def base_dir(self) -> Path:
        return self.path.parent


def read_key_value_config(
    config_path: str | Path,
    *,
    allowed_keys: Collection[str] | None = None,
) -> ConfigDocument:
    """Read a strict UTF-8 ``key=value`` file.

    Blank lines and full-line ``#`` comments are ignored. Keys are normalized to
    lower case. Duplicate, empty, or unknown keys fail before a workflow starts.
    Inline comments are intentionally unsupported, so ``#`` remains part of a
    value when it does not begin the stripped line.
    """

    path = Path(config_path).expanduser().resolve(strict=True)
    allowed = {key.strip().lower() for key in allowed_keys} if allowed_keys is not None else None
    values: dict[str, str] = {}
    for lineno, source in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = source.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{lineno}: expected key=value")
        key, value = line.split("=", 1)
        key = key.strip().lower()
        if not key:
            raise ValueError(f"{path}:{lineno}: key must not be empty")
        if key in values:
            raise ValueError(f"{path}:{lineno}: duplicate key {key!r}")
        if allowed is not None and key not in allowed:
            raise ValueError(f"{path}:{lineno}: unknown key {key!r}")
        values[key] = value.strip()
    return ConfigDocument(path=path, values=values)


def require_values(
    values: Mapping[str, str],
    required_keys: Collection[str],
    *,
    source: str | Path | None = None,
) -> None:
    """Require non-empty values for every named key."""

    missing = [key for key in required_keys if not values.get(key, "").strip()]
    if missing:
        prefix = f"{source}: " if source is not None else ""
        raise ValueError(f"{prefix}missing required keys: {missing}")


def parse_on_off(value: str, *, key: str) -> bool:
    """Convert the public ON/OFF spelling to bool."""

    normalized = value.strip().lower()
    if normalized not in {"on", "off"}:
        raise ValueError(f"{key} must be ON or OFF")
    return normalized == "on"


def parse_int(value: str, *, key: str, minimum: int | None = None) -> int:
    """Parse an integer and optionally enforce an inclusive lower bound."""

    try:
        result = int(value.strip())
    except ValueError as exc:
        raise ValueError(f"{key} must be an integer") from exc
    if minimum is not None and result < minimum:
        raise ValueError(f"{key} must be at least {minimum}")
    return result


def parse_float(value: str, *, key: str) -> float:
    """Parse a floating-point value."""

    try:
        return float(value.strip())
    except ValueError as exc:
        raise ValueError(f"{key} must be a number") from exc
