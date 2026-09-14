#!/usr/bin/env python3
"""Fetch unimod.obo and build data/unimod.csv + src/IsoSpec++/unimod_table_data.h.

Same source/approach as necromerge2's git/sage/scripts/build_unimod_table.py
(https://www.unimod.org/obo/unimod.obo), extended to also capture each term's
xref: delta_composition and convert it into IsoSpec's own native concatenated
formula-string grammar (e.g. "H(3) C(2) N O" -> "H3C2N1O1", "H(-2) O(-1)" ->
"H-2O-1") so loading a table entry's composition is a literal call to
parse_formula (isoSpec++.cpp), not a second bespoke tokenizer. Every element's
count is written out explicitly, including implied count-1 ("N" -> "N1") --
this is load-bearing: parse_formula's scanner greedily consumes all
consecutive letters as one element name bounded by digits on both sides, so a
2-letter symbol (Co, Cl, Ca, ...) sitting next to a 1-letter one in the same
composition is only unambiguous if no element's count is ever omitted.

Only ships entries IsoSpec can represent as a pure element-count delta:
  - skips any term with no delta_composition at all
  - skips any term with an isotope-labeled token (e.g. "2H(8)", "13C(3)") --
    IsoSpec has no pinned-isotope pseudo-elements
  - skips any term with a non-element "brick" token (Unimod's glycan/
    derivatization shorthand: Hex, HexNAc, dHex, NeuAc, Sulf, Ac, ...) --
    each needs its own monosaccharide/group-to-formula expansion table, not
    built here
Verified against a 2026-09 snapshot: 1560 terms have a delta_composition,
979 survive this filter (~63%), including every practically-important common
PTM (Carbamidomethyl, Oxidation, Phospho, Acetyl, Methyl, GG, Deamidated).

Re-run this script to refresh both output files against a newer Unimod
release; they are checked in, not fetched/generated at build time.
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
HEADER_OUTPUT_PATH = REPO_ROOT / "src" / "IsoSpec++" / "unimod_table_data.h"

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


def convert_composition(raw: str, known_symbols: set[str]) -> str | None:
    """Convert a Unimod delta_composition string into IsoSpec's native
    concatenated formula grammar, or return None if any token is an isotope
    label or not a symbol IsoSpec knows (glycan/derivatization "bricks" like
    Hex/HexNAc/dHex/NeuAc/Sulf/Ac, or anything else IsoSpec can't represent)."""
    out_parts = []
    for raw_token in raw.split():
        if ISOTOPE_TOKEN_RE.match(raw_token):
            return None
        m = TOKEN_RE.match(raw_token)
        if not m:
            return None
        symbol, count = m.group(1), m.group(2)
        if symbol not in known_symbols:
            return None
        count = int(count) if count is not None else 1
        out_parts.append(f"{symbol}{count}")
    if not out_parts:
        return None
    return "".join(out_parts)


def parse_obo(text: str, known_symbols: set[str]) -> list[tuple[int, str, float, str]]:
    term_list = text.split("[Term]")
    term_list.pop(0)  # header/version block, not a term

    rows = []
    skipped_no_composition = 0
    skipped_unsupported = 0
    for term in term_list:
        id_match = ID_RE.search(term)
        mass_match = MASS_RE.search(term)
        composition_match = COMPOSITION_RE.search(term)
        if not id_match or not mass_match:
            continue  # e.g. UNIMOD:0 root node has no delta_mono_mass
        if not composition_match:
            skipped_no_composition += 1
            continue
        composition = convert_composition(composition_match.group(1), known_symbols)
        if composition is None:
            skipped_unsupported += 1
            continue
        name_match = NAME_RE.search(term)
        name = name_match.group(1).strip() if name_match else ""
        rows.append((int(id_match.group(1)), name, float(mass_match.group(1)), composition))

    print(
        f"skipped (no delta_composition): {skipped_no_composition}; "
        f"skipped (isotope-labeled or non-element token): {skipped_unsupported}; "
        f"kept: {len(rows)}",
        file=sys.stderr,
    )
    return rows


def write_csv(rows: list[tuple[int, str, float, str]]) -> None:
    CSV_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_OUTPUT_PATH.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "name", "mono_mass", "composition"])
        writer.writerows(rows)
    print(f"wrote {CSV_OUTPUT_PATH}", file=sys.stderr)


def write_header() -> None:
    """Embed data/unimod.csv verbatim as a C++ raw string literal, the same
    checked-in-generated-data pattern element_tables.cpp already uses for the
    periodic table -- no build-time codegen step needed."""
    csv_text = CSV_OUTPUT_PATH.read_text()
    if ")UNIMODCSV\"" in csv_text:
        raise ValueError("data/unimod.csv unexpectedly contains the raw-string delimiter")
    header = f"""// Generated by scripts/build_unimod_table.py -- do not edit by hand.
// Embeds data/unimod.csv verbatim; re-run the script to regenerate both
// together (see that file's docstring and docs/ai/unimod.md).
#pragma once

namespace IsoSpec {{

static const char* const kEmbeddedUnimodCsv = R"UNIMODCSV(
{csv_text})UNIMODCSV";

}}  // namespace IsoSpec
"""
    HEADER_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    HEADER_OUTPUT_PATH.write_text(header)
    print(f"wrote {HEADER_OUTPUT_PATH}", file=sys.stderr)


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
    write_header()


if __name__ == "__main__":
    main()
