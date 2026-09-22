/*
 *   Copyright (C) 2015-2026 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */

#pragma once

#include <vector>

#include "composition.h"
#include "isoSpec++.h"
#include "unimod.h"

namespace IsoSpec {

//! Parses `sequence` for [UNIMOD:<id>] modification brackets in exactly the
//! notation this monorepo's SAGE fork writes (crates/sage/src/peptide.rs):
//! `[UNIMOD:<id>]-SEQUENCE` for an N-terminal mod, `X[UNIMOD:<id>]` for an
//! internal one (immediately after the modified residue), `SEQUENCE-
//! [UNIMOD:<id>]` for a C-terminal one. Parsing needs no position-awareness
//! to get this right: a `[UNIMOD:<id>]` bracket has the same effect (add
//! that modification's composition delta) wherever it appears, and the
//! surrounding '-' either side of a terminal bracket is just an ordinary
//! ignored character, exactly as it always was in a plain sequence (see
//! parse_fasta's "unrecognized characters are silently ignored" contract,
//! which this function preserves for everything except '[' itself -- a
//! literal '[' or a digit outside a bracket was never valid input before
//! either way).
//!
//! Base residue composition comes from fasta.h's existing
//! aa_symbol_to_elem_counts (CHNOSSe); each recognized mod's resolved
//! composition (from `mods`) is added on top, so the result can include
//! elements no plain amino-acid sequence ever needs (P, Br, Ca, Fe, ...).
//!
//! Throws std::invalid_argument on a malformed bracket (missing "UNIMOD:"
//! prefix, missing/non-numeric id, missing closing ']') or an id not present
//! in `mods` -- covers both a genuinely-unknown id and one of the ids
//! data/unimod.csv deliberately excludes (isotope-labeled, glycan/
//! derivatization "brick" entries) alike, same error, no special-casing.
ElementComposition parse_fasta_with_mods(const char* sequence, const UnimodTable& mods = embedded_unimod_table());

//! As above, plus terminal H2O -- mirrors parse_fasta_full's role for the
//! mods-unaware path.
ElementComposition parse_fasta_with_mods_full(const char* sequence, const UnimodTable& mods = embedded_unimod_table());

//! Same as parse_fasta_with_mods, but fills a caller-owned `out` in place
//! (cleared first) instead of returning a fresh ElementComposition -- the
//! allocation-free form for a hot loop over many sequences (a single `out`
//! reused across the whole loop reallocates at most a handful of times, not
//! once per call). parse_fasta_with_mods/_full are thin convenience wrappers
//! around this pair for one-off callers.
void parse_fasta_with_mods_into(const char* sequence, ElementComposition& out,
                                 const UnimodTable& mods = embedded_unimod_table());

//! As above, plus terminal H2O.
void parse_fasta_with_mods_full_into(const char* sequence, ElementComposition& out,
                                      const UnimodTable& mods = embedded_unimod_table());

}  // namespace IsoSpec
