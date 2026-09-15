#!/usr/bin/env python3
"""Generate data/unimod.csv from Unimod OBO metadata.

Each original ID has a row. Empty composition means unsupported; reason explains
why. Supported compositions use IsoSpec's formula grammar, including signed
counts. Python and R package sources link to this CSV; releases contain its bytes.
Use --input for an offline refresh from a saved OBO snapshot.
"""
from __future__ import annotations

import csv
import re
import sys
import urllib.request
from pathlib import Path

UNIMOD_OBO_URL = "https://www.unimod.org/obo/unimod.obo"
REPO_ROOT = Path(__file__).resolve().parent.parent
CSV_OUTPUT_PATH = REPO_ROOT / "data" / "unimod.csv"

ID_RE = re.compile(r"^id: UNIMOD:(\d+)$", re.MULTILINE)
NAME_RE = re.compile(r"^name: (.+)$", re.MULTILINE)
MASS_RE = re.compile(r'^xref: delta_mono_mass "([-\d.]+)"$', re.MULTILINE)
COMPOSITION_RE = re.compile(r'^xref: delta_composition "([^"]+)"$', re.MULTILINE)

# One composition token: an element symbol, optionally followed by a signed
# count in parentheses (absent means count 1), e.g. "H(3)", "O(-1)", "N".
TOKEN_RE = re.compile(r"^([A-Za-z][a-z]?)(?:\((-?\d+)\))?$")

# A token whose symbol part starts with a digit is an isotope label (e.g.
# "2H", "13C") -- IsoSpec has no pinned-isotope pseudo-elements for these.
ISOTOPE_TOKEN_RE = re.compile(r"^\d")


def known_element_symbols() -> set[str]:
    """IsoSpec's own supported element symbols, read from the installed
    IsoSpecPy package rather than hardcoded, so this filter can't silently
    drift from element_tables.h."""
    from IsoSpecPy import PeriodicTbl  # noqa: PLC0415 -- deliberately lazy/optional dep of this script only

    return set(PeriodicTbl.symbol_to_masses.keys())


def convert_composition(raw: str, known_symbols: set[str]) -> str:
    """Convert a Unimod delta_composition string into IsoSpec's native
    concatenated formula grammar. Raise ValueError with the exclusion reason
    if a token cannot be represented."""
    out_parts = []
    for raw_token in raw.split():
        if ISOTOPE_TOKEN_RE.match(raw_token):
            raise ValueError(f"isotope-labeled token {raw_token}; pinned-isotope conversion not supported")
        m = TOKEN_RE.match(raw_token)
        if not m:
            raise ValueError(f"non-element token {raw_token}; group-to-formula expansion not supported")
        symbol, count = m.group(1), m.group(2)
        if symbol not in known_symbols:
            raise ValueError(f"token {raw_token} has no supported element symbol or group expansion")
        count = int(count) if count is not None else 1
        out_parts.append(f"{symbol}{count}")
    if not out_parts:
        raise ValueError("empty delta_composition")
    return "".join(out_parts)


def parse_obo(text: str, known_symbols: set[str]) -> list[tuple]:
    records = {}
    for term in text.split("[Term]")[1:]:
        id_match = ID_RE.search(term)
        if not id_match:
            continue
        id = int(id_match.group(1))
        name_match = NAME_RE.search(term)
        mass_match = MASS_RE.search(term)
        composition_match = COMPOSITION_RE.search(term)
        name = name_match.group(1).strip() if name_match else ""
        mass = float(mass_match.group(1)) if mass_match else ""
        composition, reason = "", ""
        if not mass_match:
            reason = "ontology root, not a modification" if id == 0 else "missing delta_mono_mass"
        elif not composition_match:
            reason = "missing delta_composition"
        else:
            try:
                composition = convert_composition(composition_match.group(1), known_symbols)
            except ValueError as error:
                reason = str(error)
        records[id] = (id, name, mass, composition, reason)
    return [records.get(id, (id, "", "", "", "not present in source Unimod snapshot"))
            for id in range(max(records, default=0) + 1)]


def write_csv(rows: list[tuple]) -> None:
    CSV_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_OUTPUT_PATH.open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["id", "name", "mono_mass", "composition", "reason"])
        writer.writerows(rows)
    print(f"wrote {CSV_OUTPUT_PATH}: {len(rows)} IDs, "
          f"{sum(bool(row[3]) for row in rows)} supported", file=sys.stderr)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Local unimod.obo file to use instead of fetching UNIMOD_OBO_URL "
        "(for offline reruns / testing).",
    )
    args = parser.parse_args()

    if args.input is not None:
        print(f"reading {args.input}", file=sys.stderr)
        text = args.input.read_text()
    else:
        print(f"fetching {UNIMOD_OBO_URL}", file=sys.stderr)
        with urllib.request.urlopen(UNIMOD_OBO_URL) as response:
            text = response.read().decode("utf-8")

    known_symbols = known_element_symbols()
    rows = parse_obo(text, known_symbols)
    rows.sort(key=lambda r: r[0])

    write_csv(rows)


if __name__ == "__main__":
    main()
