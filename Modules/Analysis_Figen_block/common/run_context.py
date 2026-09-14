"""Create isolated, traceable run directories below work_results."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

from .paths import get_results_root, result_run_path, validate_component


@dataclass(frozen=True)
class RunContext:
    feature: tuple[str, ...]
    run_id: str
    results_root: Path
    run_dir: Path
    inputs_dir: Path
    intermediate_dir: Path
    tables_dir: Path
    logs_dir: Path
    manifest_path: Path

    def as_dict(self) -> dict[str, str | list[str]]:
        return {
            "feature": list(self.feature),
            "run_id": self.run_id,
            "results_root": str(self.results_root),
            "run_dir": str(self.run_dir),
            "inputs_dir": str(self.inputs_dir),
            "intermediate_dir": str(self.intermediate_dir),
            "tables_dir": str(self.tables_dir),
            "logs_dir": str(self.logs_dir),
            "manifest_path": str(self.manifest_path),
        }


def generate_run_id(
    target: str,
    purpose: str,
    *,
    when: datetime | None = None,
) -> str:
    """Generate <target>_<purpose>_<YYYYMMDD_HHMMSS>."""

    safe_target = validate_component(target, label="target")
    safe_purpose = validate_component(purpose, label="purpose")
    timestamp = (when or datetime.now()).strftime("%Y%m%d_%H%M%S")
    return f"{safe_target}_{safe_purpose}_{timestamp}"


def create_run_context(
    feature: str | Iterable[str],
    run_id: str,
    *,
    results_root: str | Path | None = None,
    resume: bool = False,
) -> RunContext:
    """Create or explicitly resume one isolated result run."""

    boundary = (
        Path(results_root).expanduser().resolve(strict=False)
        if results_root is not None
        else get_results_root(create=True)
    )
    boundary.mkdir(parents=True, exist_ok=True)
    run_dir = result_run_path(boundary, feature, run_id)
    existed = run_dir.exists()
    if existed and not resume:
        raise FileExistsError(
            f"Run directory already exists; pass resume=True explicitly: {run_dir}"
        )
    run_dir.mkdir(parents=True, exist_ok=resume)
    directories = {
        name: run_dir / name
        for name in ("inputs", "intermediate", "tables", "logs")
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    feature_parts = run_dir.relative_to(boundary).parts[:-1]
    return RunContext(
        feature=tuple(feature_parts),
        run_id=run_id,
        results_root=boundary,
        run_dir=run_dir,
        inputs_dir=directories["inputs"],
        intermediate_dir=directories["intermediate"],
        tables_dir=directories["tables"],
        logs_dir=directories["logs"],
        manifest_path=run_dir / "run_manifest.json",
    )
