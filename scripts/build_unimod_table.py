#!/usr/bin/env python3
"""Fetch unimod.obo and build the CSV, embedded data, and C++ support ledger.

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

Re-run this script to refresh all output files against a newer Unimod
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
SUPPORT_OUTPUT_PATH = REPO_ROOT / "src" / "IsoSpec++" / "unimod_support.h"

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


def parse_obo(text: str, known_symbols: set[str]) -> tuple[list[tuple[int, str, float, str]], dict[int, str]]:
    term_list = text.split("[Term]")
    term_list.pop(0)  # header/version block, not a term

    rows = []
    exclusions = {}
    skipped_no_composition = 0
    skipped_unsupported = 0
    for term in term_list:
        id_match = ID_RE.search(term)
        mass_match = MASS_RE.search(term)
        composition_match = COMPOSITION_RE.search(term)
        if not id_match:
            continue
        id = int(id_match.group(1))
        if not mass_match:
            exclusions[id] = "ontology root, not a modification" if id == 0 else "missing delta_mono_mass"
            continue
        if not composition_match:
            exclusions[id] = "missing delta_composition"
            skipped_no_composition += 1
            continue
        try:
            composition = convert_composition(composition_match.group(1), known_symbols)
        except ValueError as error:
            exclusions[id] = str(error)
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
    return rows, exclusions


def write_csv(rows: list[tuple[int, str, float, str]]) -> None:
    CSV_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_OUTPUT_PATH.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "name", "mono_mass", "composition"])
        writer.writerows(rows)
    print(f"wrote {CSV_OUTPUT_PATH}", file=sys.stderr)


# MSVC's C2026 ("string too big, trailing characters truncated") fires on a
# single large raw string literal well under what its documented 65535-byte
# limit would suggest -- empirically, a single ~35KB embedded-CSV literal
# already tripped it in CI (windows-11-arm, 2026-09-14). GCC/Clang have no
# such limit, so splitting costs them nothing. Chunk well under any
# plausible reading of the limit, split only at line boundaries so each
# chunk stays valid CSV text on its own.
_CHUNK_LIMIT_BYTES = 4000


def _chunk_csv_text(csv_text: str) -> list[str]:
    chunks: list[str] = []
    current: list[str] = []
    current_len = 0
    for line in csv_text.splitlines(keepends=True):
        if current and current_len + len(line) > _CHUNK_LIMIT_BYTES:
            chunks.append("".join(current))
            current = []
            current_len = 0
        current.append(line)
        current_len += len(line)
    if current:
        chunks.append("".join(current))
    return chunks


def write_header() -> None:
    """Embed data/unimod.csv as an array of C++ raw string literal chunks
    (concatenated into one std::string at load time, unimod.cpp's
    embedded_unimod_table()) -- the same checked-in-generated-data pattern
    element_tables.cpp already uses for the periodic table, no build-time
    codegen step needed, just split into pieces small enough for every
    compiler this library targets."""
    csv_text = CSV_OUTPUT_PATH.read_text()
    if ")UNIMODCSV\"" in csv_text:
        raise ValueError("data/unimod.csv unexpectedly contains the raw-string delimiter")

    chunks = _chunk_csv_text(csv_text)
    chunk_literals = ",\n".join(f'    R"UNIMODCSV({chunk})UNIMODCSV"' for chunk in chunks)

    header = f"""// Generated by scripts/build_unimod_table.py -- do not edit by hand.
// Embeds data/unimod.csv, split into chunks well under MSVC's single-string-
// literal size limit (see _CHUNK_LIMIT_BYTES in that script for why this is
// chunked at all -- a single large literal here previously tripped MSVC's
// C2026 "string too big" diagnostic in CI). Concatenated into one
// std::string at load time by unimod.cpp's embedded_unimod_table(). Re-run
// the script to regenerate both data/unimod.csv and this file together (see
// that file's docstring and docs/ai/unimod.md).
#pragma once

#include <cstddef>

namespace IsoSpec {{

static const char* const kEmbeddedUnimodCsvChunks[] = {{
{chunk_literals}
}};
static const std::size_t kEmbeddedUnimodCsvChunkCount =
    sizeof(kEmbeddedUnimodCsvChunks) / sizeof(kEmbeddedUnimodCsvChunks[0]);

}}  // namespace IsoSpec
"""
    HEADER_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    HEADER_OUTPUT_PATH.write_text(header)
    print(f"wrote {HEADER_OUTPUT_PATH} ({len(chunks)} chunks)", file=sys.stderr)


def write_support_header(exclusions: dict[int, str]) -> None:
    """Generate the ledger from the embedded CSV, with reasons from its OBO source."""
    with CSV_OUTPUT_PATH.open(newline="") as f:
        supported_ids = {int(row["id"]) for row in csv.DictReader(f)}
    size = max(supported_ids, default=0) + 1
    lines = []
    for id in range(size):
        if id in supported_ids:
            lines.append(f"    true,   // UNIMOD:{id}")
        else:
            reason = exclusions.get(id, "not present in source Unimod snapshot")
            lines.append(f"    false,  // UNIMOD:{id}: {reason}")
    values_text = "\n".join(lines)
    header = f"""// Generated from data/unimod.csv and its Unimod OBO source by scripts/build_unimod_table.py.
// Index = original Unimod ID; true = supported by the embedded IsoSpec table.
#pragma once

#include <cstddef>
#include <cstdint>

namespace IsoSpec {{

inline constexpr bool unimod_supported[] = {{
{values_text}
}};

//! Tests support in the shipped table without loading or parsing it.
//! Unsupported and unknown IDs return false. Override CSVs do not affect this
//! ledger. A supported modification still needs a compatible base composition.
inline constexpr bool is_unimod_supported(std::uint64_t id) noexcept {{
    return id < sizeof(unimod_supported) / sizeof(unimod_supported[0])
        && unimod_supported[static_cast<std::size_t>(id)];
}}

}}  // namespace IsoSpec
"""
    SUPPORT_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    SUPPORT_OUTPUT_PATH.write_text(header)
    print(f"wrote {SUPPORT_OUTPUT_PATH} ({size} slots)", file=sys.stderr)


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
    rows, exclusions = parse_obo(text, known_symbols)
    rows.sort(key=lambda r: r[0])

    write_csv(rows)
    write_header()
    write_support_header(exclusions)


if __name__ == "__main__":
    main()
