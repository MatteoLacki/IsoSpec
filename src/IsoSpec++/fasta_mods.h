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

#include "isoSpec++.h"
#include "unimod.h"

namespace IsoSpec {

//! A resolved elemental composition: parallel arrays of
//! (element_table_first_index, count) pairs, one per element actually
//! present (non-zero net count).
//!
//! Counts are *signed*, and parse_fasta_with_mods/_full really can return a
//! negative one: a modification's composition delta may remove more atoms of
//! some element than the bare sequence provides ("G[UNIMOD:11]" nets H-1
//! S-1). That is a legitimate intermediate -- adding a second modification
//! could bring it back up -- so the parser does not reject it. The
//! non-negativity check belongs at the point a composition becomes a
//! molecule, and lives in build_iso_from_composition; anything else
//! consuming an ElementComposition directly (the C ABI's
//! parseFastaWithModsC, R's RParsePeptideSequence) must do its own check
//! before handing the counts to an Iso.
struct ElementComposition {
    std::vector<int> element_first_index;
    std::vector<int> count;

    //! Empties both vectors without releasing their capacity -- a caller
    //! that reuses one ElementComposition across many parse_fasta_with_mods_into
    //! calls (e.g. building a cache over millions of peptides) pays for the
    //! underlying heap allocation at most a handful of times total (until
    //! capacity reaches this feature's real ceiling, ~30 elements), not once
    //! per call the way the value-returning parse_fasta_with_mods does.
    void clear() {
        element_first_index.clear();
        count.clear();
    }
};

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

//! A composition expanded to the flat, per-isotope form the generic
//! Iso(dimNumber, isotopeNumbers, atomCounts, masses, probs) constructor
//! takes: one isotopeNumbers/atomCounts entry per element of the
//! composition, and one isotope_masses/isotope_probabilities entry per
//! isotope, element-major, in the composition's own order. Exactly the
//! layout the C ABI's setupIso expects, since that constructor is what it
//! forwards to.
struct ExpandedComposition {
    std::vector<int> isotopeNumbers;
    std::vector<int> atomCounts;
    std::vector<double> isotope_masses;
    std::vector<double> isotope_probabilities;

    //! Empties all four vectors without releasing their capacity -- same
    //! reuse contract, for the same reason, as ElementComposition::clear
    //! above.
    void clear() {
        isotopeNumbers.clear();
        atomCounts.clear();
        isotope_masses.clear();
        isotope_probabilities.clear();
    }
};

//! Resolves each element of `composition` to its isotopes' masses and
//! probabilities from element_tables.h, filling a caller-owned `out`
//! (cleared first).
//!
//! This is the one implementation of that resolution a binding language is
//! meant to reach: before it existed, IsoSpecPy carried its own copy
//! (IsoParamsFromDict, walking PeriodicTbl's symbol->isotopes dicts) purely
//! because the C ABI offered nothing between setupIso -- which demands that
//! the caller already know every isotope mass -- and the sequence-only
//! isoFromFasta. cwrapper.h's expandCompositionC exposes this, and that
//! copy is gone.
//!
//! Counts are copied through unchanged, negatives included: expanding a
//! composition is not the same as making a molecule out of one, and the
//! non-negativity check belongs at the latter -- see ElementComposition's
//! note above, and build_iso_from_composition, which does check.
void expand_composition_into(const ElementComposition& composition, ExpandedComposition& out,
                             bool use_nominal_masses = false);

//! Resolves a composition (as returned by parse_fasta_with_mods[_full]) into
//! an Iso via the same generic Iso(dimNumber, isotopeNumbers, atomCounts,
//! masses, probs) constructor parse_formula already uses. Throws
//! std::invalid_argument if any accumulated element count is negative (a
//! modification set that removes more atoms of some element than the base
//! sequence plus other modifications provide).
Iso build_iso_from_composition(const ElementComposition& composition, bool use_nominal_masses = false);

}  // namespace IsoSpec
