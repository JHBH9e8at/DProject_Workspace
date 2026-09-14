"""Filter an AHC population and select its best own-state docking results.

Steps:
1. keep valid and unique molecules;
2. canonicalize SMILES with RDKit;
3. keep finite, nonzero docking scores with usable variant IDs;
4. remove canonical duplicates, retaining the most negative score;
5. sort by own-state docking score and select the configured top percentage;
6. optionally override percentage selection with a test-mode top N;
7. annotate PAINS and BRENK matches without removing them.
"""

import argparse
import math
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")


def _true_mask(series):
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def canonicalize_smiles(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None, None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        canonical = Chem.MolToSmiles(mol, canonical=True)
        return canonical, mol
    except Exception:
        return None, None


def annotate_structural_alerts(table, mol_col="_mol"):
    pains_params = FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog(pains_params)

    brenk_params = FilterCatalogParams()
    brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    brenk_catalog = FilterCatalog(brenk_params)

    table = table.copy()
    table["PAINS"] = table[mol_col].apply(pains_catalog.HasMatch)
    table["BRENK"] = table[mol_col].apply(brenk_catalog.HasMatch)
    return table


def filter_population(
    input_csv,
    output_csv,
    run_name,
    best_dscore_top_per=1.0,
    testmode_top_n=None,
):
    """Write a filtered, score-ranked AHC population and return it."""
    run_name = (run_name or "PPS").strip().upper()
    score_col = f"{run_name}_r_i_docking_score"
    variant_col = f"{run_name}_best_variant"
    required = {"valid", "unique", "smiles", score_col, variant_col}

    table = pd.read_csv(input_csv)
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Missing required filtering column(s): {missing}")
    stale_selection_columns = [
        column
        for column in ("crossdock_selection_rank", "canon_smiles", "PAINS", "BRENK")
        if column in table.columns
    ]
    if stale_selection_columns:
        table = table.drop(columns=stale_selection_columns)
    if not 0 < best_dscore_top_per <= 100:
        raise ValueError("best_dscore_top_per must be greater than 0 and at most 100")
    if testmode_top_n is not None and testmode_top_n <= 0:
        raise ValueError("testmode_top_n must be a positive integer")

    steps = []

    def apply_filter(name, working, mask):
        filtered = working.loc[mask].copy()
        steps.append((name, len(working), len(filtered)))
        return filtered

    working = apply_filter("valid molecule", table, _true_mask(table["valid"]))
    working = apply_filter("unique molecule", working, _true_mask(working["unique"]))

    canonicalized = working["smiles"].apply(canonicalize_smiles)
    working["canon_smiles"] = canonicalized.apply(lambda result: result[0])
    working["_mol"] = canonicalized.apply(lambda result: result[1])
    working = apply_filter(
        "canonical SMILES",
        working,
        working["canon_smiles"].notna(),
    )

    working[score_col] = pd.to_numeric(working[score_col], errors="coerce")
    usable_variant = (
        working[variant_col].notna()
        & working[variant_col].astype(str).str.strip().ne("")
        & working[variant_col].astype(str).str.strip().ne("0.0")
    )
    finite_score = working[score_col].map(
        lambda value: pd.notna(value) and math.isfinite(value)
    )
    usable_result = finite_score & working[score_col].ne(0)
    working = apply_filter(
        "usable docking result",
        working,
        usable_result & usable_variant,
    )

    if working.empty:
        raise ValueError("No molecules remain after filtering")

    before_duplicates = len(working)
    best_indices = working.groupby("canon_smiles")[score_col].idxmin()
    working = working.loc[best_indices].copy()
    steps.append(("canonical duplicate removal", before_duplicates, len(working)))

    working = working.sort_values(score_col, ascending=True).reset_index(drop=True)
    working.insert(0, "crossdock_selection_rank", range(1, len(working) + 1))
    before_selection = len(working)
    if testmode_top_n is not None:
        selection_count = min(testmode_top_n, before_selection)
        selection_name = f"test-mode top-{testmode_top_n} override"
        selection_mode = "testmode_top_n"
    else:
        selection_count = max(
            1,
            math.ceil(before_selection * best_dscore_top_per / 100.0),
        )
        selection_name = f"top-{best_dscore_top_per:g}% own-state selection"
        selection_mode = "best_dscore_top_per"
    working = working.head(selection_count).copy()
    steps.append((selection_name, before_selection, len(working)))

    working = annotate_structural_alerts(working)
    working = working.drop(columns=["_mol"])

    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    working.to_csv(output_csv, index=False)

    log_file = output_csv.with_name(f"{output_csv.stem}_filterlog.txt")
    with log_file.open("w", encoding="utf-8") as handle:
        handle.write(f"input: {Path(input_csv)}\n")
        handle.write(f"output: {output_csv}\n")
        handle.write(f"run_name: {run_name}\n")
        handle.write(f"score_column: {score_col}\n")
        handle.write(f"selection_mode: {selection_mode}\n")
        handle.write(f"best_dscore_top_per: {best_dscore_top_per:g}\n")
        handle.write(
            f"testmode_top_n: "
            f"{'' if testmode_top_n is None else testmode_top_n}\n"
        )
        handle.write(f"selected_count: {len(working)}\n")
        handle.write(f"initial_count: {len(table)}\n\n")
        for name, before, after in steps:
            handle.write(f"[{name}]\n")
            handle.write(f"  before: {before}\n")
            handle.write(f"  after: {after}\n")
            handle.write(f"  removed: {before - after}\n\n")
        handle.write(f"PAINS_count: {int(working['PAINS'].sum())}\n")
        handle.write(f"BRENK_count: {int(working['BRENK'].sum())}\n")
        handle.write(f"final_count: {len(working)}\n")

    for name, before, after in steps:
        print(f"{name}: {before} -> {after} (removed {before - after})")
    print(f"PAINS matches: {int(working['PAINS'].sum())}")
    print(f"BRENK matches: {int(working['BRENK'].sum())}")
    print(f"Filtered population: {output_csv}")
    print(f"Filtering log: {log_file}")
    return working


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Raw AHC scores CSV")
    parser.add_argument("--output", required=True, help="Filtered output CSV")
    parser.add_argument("--run-name", required=True, choices=["PPS", "PR"])
    parser.add_argument(
        "--best-dscore-top-per",
        type=float,
        default=1.0,
        help="Top own-state docking-score percentage selected after filtering",
    )
    parser.add_argument(
        "--testmode-top-n",
        type=int,
        help="Testing override: select exactly the top N instead of a percentage",
    )
    args = parser.parse_args(argv)
    filter_population(
        input_csv=args.input,
        output_csv=args.output,
        run_name=args.run_name,
        best_dscore_top_per=args.best_dscore_top_per,
        testmode_top_n=args.testmode_top_n,
    )


if __name__ == "__main__":
    main()
