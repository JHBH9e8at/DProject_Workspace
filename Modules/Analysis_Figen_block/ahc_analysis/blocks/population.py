"""Pure population merging and cleaning operations migrated from PY-019/PY-020."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


SCRIPT_ID = "PY-102"


def discover_iteration_files(iterations_dir: str | Path) -> tuple[list[Path], Path]:
    """Return merge inputs and the deliberately excluded final score file."""

    files = sorted(Path(iterations_dir).glob("*_scores.csv"))
    if len(files) < 2:
        raise ValueError("At least two *_scores.csv files are required")
    return files[:-1], files[-1]


def merge_iteration_scores(iterations_dir: str | Path) -> tuple[pd.DataFrame, list[Path], Path]:
    """Merge all lexically sorted iteration score files except the final file."""

    included, excluded = discover_iteration_files(iterations_dir)
    merged = pd.concat((pd.read_csv(path) for path in included), ignore_index=True)
    return merged, included, excluded


def _rdkit_components() -> tuple[Any, Any, Any]:
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    except ImportError as exc:
        raise RuntimeError("AHC population cleaning requires RDKit") from exc
    RDLogger.DisableLog("rdApp.*")
    return Chem, FilterCatalog, FilterCatalogParams


def clean_population(
    frame: pd.DataFrame,
    *,
    run_name: str = "PPS",
) -> tuple[pd.DataFrame, list[dict[str, int | str]]]:
    """Apply the legacy validity, uniqueness, canonicalization, and alert filters."""

    required = {"valid", "unique", "smiles", f"{run_name}_r_i_docking_score"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    Chem, FilterCatalog, FilterCatalogParams = _rdkit_components()
    steps: list[dict[str, int | str]] = []

    def record(name: str, before: int, after: int) -> None:
        steps.append({"step": name, "before": before, "after": after, "removed": before - after})

    valid = frame[frame["valid"] == True].copy()  # preserve legacy equality semantics
    record("valid filter", len(frame), len(valid))
    unique = valid[valid["unique"] == True].copy()
    record("unique filter", len(valid), len(unique))

    def canonicalize(smiles: object) -> tuple[str | None, Any | None]:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, None
            canonical = Chem.MolToSmiles(mol, canonical=True)
            return canonical, Chem.MolFromSmiles(canonical)
        except Exception:
            return None, None

    canonicalized = unique.copy()
    converted = canonicalized["smiles"].apply(canonicalize)
    canonicalized["canon_smiles"] = converted.apply(lambda item: item[0])
    canonicalized["_mol"] = converted.apply(lambda item: item[1])
    canonicalized = canonicalized.dropna(subset=["canon_smiles"]).copy()
    record("canonical conversion filter", len(unique), len(canonicalized))

    score_column = f"{run_name}_r_i_docking_score"
    best_indices = canonicalized.groupby("canon_smiles")[score_column].idxmin()
    cleaned = canonicalized.loc[best_indices].copy().reset_index(drop=True)
    record("duplicate removal (best per canon_smiles)", len(canonicalized), len(cleaned))

    pains_params = FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog(pains_params)
    brenk_params = FilterCatalogParams()
    brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    brenk_catalog = FilterCatalog(brenk_params)
    cleaned["PAINS"] = cleaned["_mol"].apply(pains_catalog.HasMatch)
    cleaned["BRENK"] = cleaned["_mol"].apply(brenk_catalog.HasMatch)
    return cleaned.drop(columns=["_mol"]), steps


def format_cleaning_log(
    *, input_label: str, output_label: str, initial_count: int,
    final_count: int, steps: list[dict[str, int | str]],
) -> str:
    lines = [f"input: {input_label}", f"output: {output_label}", f"initial instance count: {initial_count}", ""]
    for step in steps:
        lines.extend([
            f"[{step['step']}]", f"  before: {step['before']}",
            f"  after:  {step['after']}", f"  removed: {step['removed']}", "",
        ])
    lines.extend([f"final instance count: {final_count}", f"total removed: {initial_count - final_count}", ""])
    return "\n".join(lines)
