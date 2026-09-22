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

// The library's vocabulary for "a molecule as element counts", and the one
// implementation of turning that into the isotope masses and probabilities an
// Iso is built from.
//
// This lives in its own header because three unrelated things need it: the
// formula parser (isoSpec++.cpp), the peptide-sequence parsers (fasta_mods.h),
// and the C ABI (cwrapper.cpp). It used to sit inside fasta_mods.h, which
// meant the core formula path could only share the resolution by including a
// peptide-modification header -- backwards, so it didn't, and kept its own
// copy instead. See docs/ai/composition_expansion.md.

#pragma once

#include <vector>

#include "isoSpec++.h"

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
//! The single implementation of that resolution. Everything that turns
//! element counts into an Iso goes through here: parse_formula (the formula
//! string path), build_iso_from_composition (the peptide-sequence path), and
//! cwrapper.h's expandCompositionC -- which exists so that a binding
//! language need not carry its own copy of the periodic table either.
//! IsoSpecPy did carry one (IsoParamsFromDict, walking PeriodicTbl's
//! symbol->isotopes dicts) until that entry point existed.
//!
//! Counts are copied through unchanged, negatives included: expanding a
//! composition is not the same as making a molecule out of one, and the
//! non-negativity check belongs at the latter -- see ElementComposition's
//! note above, and build_iso_from_composition, which does check.
void expand_composition_into(const ElementComposition& composition, ExpandedComposition& out,
                             bool use_nominal_masses = false);

//! Throws std::invalid_argument if any of `composition`'s counts is
//! negative, naming the offending elements. The check every path from a
//! composition to an Iso has to make: a negative atom count reaching Iso's
//! constructor is not merely a wrong answer but undefined behaviour --
//! Marginal::computeModeConf() sizes its configuration buffer from the count
//! and writeInitialConfiguration then walks off the end of it (ASan catches
//! it as a heap-buffer-overflow).
void reject_negative_counts(const ElementComposition& composition);

//! Resolves a composition (as returned by parse_fasta_with_mods[_full]) into
//! an Iso via the same generic Iso(dimNumber, isotopeNumbers, atomCounts,
//! masses, probs) constructor parse_formula uses. Throws
//! std::invalid_argument if any accumulated element count is negative (a
//! modification set that removes more atoms of some element than the base
//! sequence plus other modifications provide).
Iso build_iso_from_composition(const ElementComposition& composition, bool use_nominal_masses = false);

}  // namespace IsoSpec
