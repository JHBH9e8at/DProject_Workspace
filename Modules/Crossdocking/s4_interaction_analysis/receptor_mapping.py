"""Map state-specific receptor residue IDs onto a shared PPS/PR numbering."""

import re
from pathlib import Path

import pandas as pd


REQUIRED_MAPPING_COLUMNS = {
    "receptor_state",
    "original_chain",
    "original_residue_number",
    "residue_name",
    "common_residue_id",
}

RESIDUE_PATTERN = re.compile(
    r"^(?P<residue_name>[A-Za-z]+)(?P<residue_number>-?\d+[A-Za-z]?)(?:\.(?P<chain>.*))?$"
)


def load_residue_mapping(mapping_csv):
    """Load and validate a common PPS/PR residue mapping table."""
    mapping_csv = Path(mapping_csv).resolve()
    if not mapping_csv.is_file():
        raise FileNotFoundError(f"Residue mapping CSV not found: {mapping_csv}")

    mapping = pd.read_csv(mapping_csv, dtype=str).fillna("")
    missing = sorted(REQUIRED_MAPPING_COLUMNS - set(mapping.columns))
    if missing:
        raise ValueError(f"Residue mapping missing columns: {missing}")

    mapping["receptor_state"] = mapping["receptor_state"].str.strip().str.upper()
    invalid_states = sorted(set(mapping["receptor_state"]) - {"PPS", "PR"})
    if invalid_states:
        raise ValueError(f"Invalid receptor_state values: {invalid_states}")

    key = [
        "receptor_state",
        "original_chain",
        "original_residue_number",
        "residue_name",
    ]
    duplicated = mapping.duplicated(key, keep=False)
    if duplicated.any():
        examples = mapping.loc[duplicated, key].head(5).to_dict("records")
        raise ValueError(f"Duplicate residue mapping keys found: {examples}")
    return mapping


def parse_prolif_residue(residue_id):
    """Parse common ProLIF residue strings such as ASP129.A."""
    text = str(residue_id).strip()
    match = RESIDUE_PATTERN.match(text)
    if not match:
        return {
            "original_residue_id": text,
            "residue_name": "",
            "original_residue_number": "",
            "original_chain": "",
            "residue_parse_status": "unparsed",
        }
    values = match.groupdict()
    return {
        "original_residue_id": text,
        "residue_name": values["residue_name"].upper(),
        "original_residue_number": values["residue_number"],
        "original_chain": values["chain"] or "",
        "residue_parse_status": "parsed",
    }


def annotate_interactions(interactions, receptor_state, mapping=None):
    """Add parsed residue fields and an optional common residue identifier."""
    annotated = interactions.copy()
    state = str(receptor_state).strip().upper()
    if state not in {"PPS", "PR"}:
        raise ValueError("receptor_state must be PPS or PR")
    annotated["receptor_state"] = state

    parsed_columns = [
        "original_residue_id",
        "residue_name",
        "original_residue_number",
        "original_chain",
        "residue_parse_status",
    ]
    if annotated.empty:
        parsed = pd.DataFrame(index=annotated.index, columns=parsed_columns)
    else:
        parsed = annotated["protein_residue"].map(parse_prolif_residue).apply(pd.Series)
        parsed = parsed.reindex(columns=parsed_columns)
    annotated = pd.concat([annotated.reset_index(drop=True), parsed], axis=1)
    annotated["common_residue_id"] = pd.NA
    annotated["residue_mapping_status"] = "not_requested"

    if mapping is None:
        return annotated

    state_map = mapping[mapping["receptor_state"] == state].copy()
    keys = ["original_chain", "original_residue_number", "residue_name"]
    lookup = state_map[keys + ["common_residue_id"]]
    annotated = annotated.drop(
        columns=["common_residue_id", "residue_mapping_status"]
    ).merge(lookup, on=keys, how="left", validate="many_to_one")
    annotated["residue_mapping_status"] = annotated["common_residue_id"].notna().map(
        {True: "mapped", False: "unmapped"}
    )
    return annotated
