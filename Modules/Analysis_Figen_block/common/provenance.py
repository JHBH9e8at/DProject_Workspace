"""Read and atomically update run provenance manifests."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from .paths import assert_within
from .run_context import RunContext


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    source = Path(path).expanduser().resolve(strict=True)
    digest = hashlib.sha256()
    with source.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _input_record(path: str | Path) -> dict[str, Any]:
    source = Path(path).expanduser().resolve(strict=True)
    if not source.is_file():
        raise ValueError(f"Manifest input is not a file: {source}")
    return {
        "path": str(source),
        "size_bytes": source.stat().st_size,
        "sha256": sha256_file(source),
    }


def _output_record(path: str | Path, run_dir: Path) -> dict[str, Any]:
    output = assert_within(path, run_dir)
    record: dict[str, Any] = {"path": output.relative_to(run_dir).as_posix()}
    if output.is_file():
        record.update(
            size_bytes=output.stat().st_size,
            sha256=sha256_file(output),
        )
    else:
        record["status"] = "declared"
    return record


def build_manifest(
    context: RunContext,
    *,
    script_id: str,
    status: str = "created",
    config: dict[str, Any] | None = None,
    inputs: Iterable[str | Path] = (),
    outputs: Iterable[str | Path] = (),
    upstream_run_ids: Iterable[str] = (),
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "script_id": script_id,
        "status": status,
        "created_at_utc": utc_now(),
        "updated_at_utc": utc_now(),
        "run": context.as_dict(),
        "config": dict(config or {}),
        "upstream_run_ids": list(upstream_run_ids),
        "inputs": [_input_record(path) for path in inputs],
        "outputs": [_output_record(path, context.run_dir) for path in outputs],
    }


def write_manifest(path: str | Path, manifest: dict[str, Any]) -> Path:
    """Atomically write a UTF-8 JSON manifest."""

    target = Path(path).expanduser().resolve(strict=False)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()
    return target


def read_manifest(path: str | Path) -> dict[str, Any]:
    source = Path(path).expanduser().resolve(strict=True)
    with source.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Manifest root must be a JSON object: {source}")
    return value


def update_manifest(
    context: RunContext,
    *,
    script_id: str,
    status: str,
    config: dict[str, Any] | None = None,
    inputs: Iterable[str | Path] = (),
    outputs: Iterable[str | Path] = (),
    upstream_run_ids: Iterable[str] = (),
) -> Path:
    manifest = build_manifest(
        context,
        script_id=script_id,
        status=status,
        config=config,
        inputs=inputs,
        outputs=outputs,
        upstream_run_ids=upstream_run_ids,
    )
    return write_manifest(context.manifest_path, manifest)
