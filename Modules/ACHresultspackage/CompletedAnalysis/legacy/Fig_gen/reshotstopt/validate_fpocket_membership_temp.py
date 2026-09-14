from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(r"Q:\coding_dir\701_Project\Resultsbin\fpocket_res")

POCKETS = {
    "PPS P2": ROOT / "PPS_Protein_out" / "pockets" / "pocket2_atm.pdb",
    "PPS P38": ROOT / "PPS_Protein_out" / "pockets" / "pocket38_atm.pdb",
    "PPS P51": ROOT / "PPS_Protein_out" / "pockets" / "pocket51_atm.pdb",
    "PR P11": ROOT / "PR_potein_out" / "pockets" / "pocket11_atm.pdb",
    "PR P15": ROOT / "PR_potein_out" / "pockets" / "pocket15_atm.pdb",
    "PR P60": ROOT / "PR_potein_out" / "pockets" / "pocket60_atm.pdb",
}

DISPLAY_RESIDUES = [
    "GLU88", "ASP89", "ALA91", "THR94", "SER118", "LYS146", "ARG147",
    "TYR164", "ASP168", "GLU170", "GLU497", "GLU500", "HIS666",
    "PRO710", "ASN711", "ARG712", "ASP717", "TYR722", "LYS762",
    "GLY768", "GLU774",
]


def parse_fixed_width(path: Path) -> tuple[set[str], list[str]]:
    residues: set[str] = set()
    raw_atom_lines: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        raw_atom_lines.append(line)
        name = line[17:20].strip()
        number = line[22:26].strip()
        insertion = line[26:27].strip()
        chain = line[21:22].strip()
        residues.add(f"{name}{number}{insertion}" + (f".{chain}" if chain else ""))
    return residues, raw_atom_lines


def parse_whitespace(path: Path) -> set[str]:
    residues: set[str] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        fields = line.split()
        # Typical fpocket atom row: ATOM serial atom resname chain resnum ...
        if len(fields) < 6:
            continue
        residue_name = fields[3]
        if re.fullmatch(r"[A-Za-z]", fields[4]) and re.fullmatch(r"-?\d+[A-Za-z]?", fields[5]):
            chain = fields[4]
            residue_number = fields[5]
        else:
            chain = ""
            residue_number = fields[4]
        residues.add(f"{residue_name}{residue_number}" + (f".{chain}" if chain else ""))
    return residues


def without_chain(residues: set[str]) -> set[str]:
    return {residue.split(".")[0] for residue in residues}


def main() -> None:
    parsed: dict[str, set[str]] = {}
    print("RAW FPOCKET MEMBERSHIP AUDIT")
    print("=" * 96)
    for label, path in POCKETS.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        fixed, atom_lines = parse_fixed_width(path)
        whitespace = parse_whitespace(path)
        fixed_plain = without_chain(fixed)
        whitespace_plain = without_chain(whitespace)
        parsed[label] = fixed_plain
        print(f"\n{label}: {path}")
        print(f"  ATOM/HETATM lines: {len(atom_lines)}")
        print(f"  Unique residues (fixed-width): {len(fixed_plain)}")
        print(f"  Unique residues (whitespace):  {len(whitespace_plain)}")
        print(f"  Parser agreement: {fixed_plain == whitespace_plain}")
        if fixed_plain != whitespace_plain:
            print(f"  Fixed-only: {sorted(fixed_plain - whitespace_plain)}")
            print(f"  Whitespace-only: {sorted(whitespace_plain - fixed_plain)}")
        print("  Residues: " + ", ".join(sorted(fixed_plain, key=lambda x: int(re.search(r'\d+', x).group()))))

    print("\n" + "=" * 96)
    print("DISPLAYED RESIDUE MEMBERSHIP (directly from raw pocket*_atm.pdb)")
    header = ["Residue", *POCKETS.keys()]
    widths = [10] + [10] * len(POCKETS)
    print(" ".join(name.ljust(width) for name, width in zip(header, widths)))
    print("-" * (sum(widths) + len(widths) - 1))
    for residue in DISPLAY_RESIDUES:
        row = [residue] + [("YES" if residue in parsed[label] else "-") for label in POCKETS]
        print(" ".join(value.ljust(width) for value, width in zip(row, widths)))

    print("\nOVERLAPPING MEMBERSHIPS WITHIN EACH STATE")
    for state in ("PPS", "PR"):
        labels = [label for label in POCKETS if label.startswith(state + " ")]
        for residue in DISPLAY_RESIDUES:
            hits = [label.split()[-1] for label in labels if residue in parsed[label]]
            if len(hits) > 1:
                print(f"  {state} {residue}: {'/'.join(hits)}")


if __name__ == "__main__":
    main()
