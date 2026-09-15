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

#include "fasta_mods.h"

#include <cctype>
#include <cstring>
#include <stdexcept>
#include <string>

#include "element_lookup.h"
#include "element_tables.h"
#include "fasta.h"

namespace IsoSpec {

namespace {

//! Indices into a dense, direct-indexed accumulator: one slot per row of
//! element_tables.h's flat arrays (292), addressed by an element's *first*
//! isotope row -- every element this feature ever touches (the fixed 6 of
//! plain amino acids, plus the 27 Unimod's shipped compositions use) fits
//! this without any map. Sized at compile time, zero-initialized per call.
constexpr int kAccumulatorSize = ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES;

//! First-index lookups for the 6 elements a plain amino-acid sequence can
//! ever contribute (fasta.h's own CHNOSSe order), resolved once via
//! element_lookup.h rather than duplicating fasta.h's separate, narrower
//! aa_elem_masses/aa_isotope_numbers tables -- this path folds the base
//! sequence into the same general per-element-row accumulator a
//! modification's composition delta uses, so both add into one array.
struct AminoAcidElementIndex {
    int first_index[6];
    AminoAcidElementIndex() {
        static const char* const kSymbols[6] = {"C", "H", "N", "O", "S", "Se"};
        for (int i = 0; i < 6; i++)
            first_index[i] = find_element_table_first_index(kSymbols[i], strlen(kSymbols[i]));
    }
};

const AminoAcidElementIndex& amino_acid_element_index() {
    static const AminoAcidElementIndex instance;
    return instance;
}

void add_residue(int accumulator[kAccumulatorSize], char residue) {
    int counts[6] = {0, 0, 0, 0, 0, 0};
    const char one[2] = {residue, '\0'};
    parse_fasta(one, counts);  // zero for anything that isn't a recognized amino-acid letter

    const AminoAcidElementIndex& idx = amino_acid_element_index();
    for (int i = 0; i < 6; i++)
        if (counts[i] != 0)
            accumulator[idx.first_index[i]] += counts[i];
}

//! `p` points at '['. Parses "UNIMOD:<id>]", applies that id's composition
//! delta to `accumulator`, and advances `p` past the closing ']'. Throws
//! std::invalid_argument on anything not exactly this shape, or an id not
//! present in `mods` -- unknown and deliberately-excluded ids are the same
//! error, no special-casing.
void consume_unimod_bracket(const char*& p, const UnimodTable& mods, int accumulator[kAccumulatorSize]) {
    static const char kPrefix[] = "UNIMOD:";
    constexpr size_t kPrefixLen = sizeof(kPrefix) - 1;

    const char* q = p + 1;  // skip '['
    if (strncmp(q, kPrefix, kPrefixLen) != 0)
        throw std::invalid_argument(std::string("Malformed modification bracket (expected '[UNIMOD:<id>]'): ") + p);
    q += kPrefixLen;

    const char* digits_start = q;
    while (isdigit(static_cast<unsigned char>(*q)))
        q++;
    if (q == digits_start)
        throw std::invalid_argument(std::string("Malformed modification bracket (missing numeric id): ") + p);
    if (*q != ']')
        throw std::invalid_argument(std::string("Malformed modification bracket (missing closing ']'): ") + p);

    unsigned long id;
    try {
        id = std::stoul(std::string(digits_start, static_cast<size_t>(q - digits_start)));
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Malformed modification bracket (id out of range): ") + p);
    }

    const UnimodEntry* entry = mods.lookup(id);
    if (entry == nullptr)
        throw std::invalid_argument("Unknown or unsupported UNIMOD id in sequence: UNIMOD:" + std::to_string(id));

    for (size_t i = 0; i < entry->element_first_index.size(); i++)
        accumulator[entry->element_first_index[i]] += entry->element_delta_count[i];

    p = q + 1;  // past ']'
}

void accumulator_to_composition_into(const int accumulator[kAccumulatorSize], ElementComposition& out) {
    out.clear();
    for (int i = 0; i < kAccumulatorSize; i++) {
        if (accumulator[i] != 0) {
            out.element_first_index.push_back(i);
            out.count.push_back(accumulator[i]);
        }
    }
}

}  // namespace

void parse_fasta_with_mods_into(const char* sequence, ElementComposition& out, const UnimodTable& mods) {
    out.clear();

    if (strchr(sequence, '[') == nullptr) {
        // '[' is the only character this parser ever treats specially (it's
        // never a valid FASTA/amino-acid character on its own), so its
        // absence means the sequence cannot contain a [UNIMOD:<id>]
        // reference -- skip the general per-element accumulator below
        // entirely and go through the plain, fixed-CHNOSSe parse_fasta in
        // one pass over the whole string, matching Iso::FromFASTA's cost
        // for this (the common, unmodified) case instead of paying for
        // per-residue re-parsing plus general-element bookkeeping nothing
        // here actually needs.
        int counts[6] = {0, 0, 0, 0, 0, 0};
        parse_fasta(sequence, counts);
        const AminoAcidElementIndex& idx = amino_acid_element_index();
        for (int i = 0; i < 6; i++) {
            if (counts[i] != 0) {
                out.element_first_index.push_back(idx.first_index[i]);
                out.count.push_back(counts[i]);
            }
        }
        return;
    }

    int accumulator[kAccumulatorSize] = {0};

    const char* p = sequence;
    while (*p != '\0') {
        if (*p == '[') {
            // A [UNIMOD:<id>] bracket has the same effect wherever it
            // appears (N-terminal "[id]-SEQUENCE", internal "X[id]",
            // C-terminal "SEQUENCE-[id]") -- the '-' either side of a
            // terminal bracket is just an ordinary character, handled by
            // the branch below like any other. No position-tracking needed.
            consume_unimod_bracket(p, mods, accumulator);
        } else {
            add_residue(accumulator, *p);
            p++;
        }
    }

    accumulator_to_composition_into(accumulator, out);
}

void parse_fasta_with_mods_full_into(const char* sequence, ElementComposition& out, const UnimodTable& mods) {
    parse_fasta_with_mods_into(sequence, out, mods);

    const AminoAcidElementIndex& idx = amino_acid_element_index();
    const int h_index = idx.first_index[1];
    const int o_index = idx.first_index[3];

    bool has_h = false, has_o = false;
    for (size_t i = 0; i < out.element_first_index.size(); i++) {
        if (out.element_first_index[i] == h_index) {
            out.count[i] += 2;
            has_h = true;
        } else if (out.element_first_index[i] == o_index) {
            out.count[i] += 1;
            has_o = true;
        }
    }
    if (!has_h) {
        out.element_first_index.push_back(h_index);
        out.count.push_back(2);
    }
    if (!has_o) {
        out.element_first_index.push_back(o_index);
        out.count.push_back(1);
    }
}

ElementComposition parse_fasta_with_mods(const char* sequence, const UnimodTable& mods) {
    ElementComposition result;
    parse_fasta_with_mods_into(sequence, result, mods);
    return result;
}

ElementComposition parse_fasta_with_mods_full(const char* sequence, const UnimodTable& mods) {
    ElementComposition result;
    parse_fasta_with_mods_full_into(sequence, result, mods);
    return result;
}

Iso build_iso_from_composition(const ElementComposition& composition, bool use_nominal_masses) {
    const size_t dimNumber = composition.element_first_index.size();
    std::vector<int> isotopeNumbers(dimNumber);
    std::vector<int> atomCounts(dimNumber);
    std::vector<double> isotope_masses;
    std::vector<double> isotope_probabilities;

    const double* masses_table = use_nominal_masses ? elem_table_massNo : elem_table_mass;

    for (size_t i = 0; i < dimNumber; i++) {
        int count = composition.count[i];
        if (count < 0)
            throw std::invalid_argument(
                "Modification set removes more atoms of some element than the base sequence "
                "plus other modifications provide");
        int first_idx = composition.element_first_index[i];
        int num_isotopes = element_isotope_count(first_idx);
        isotopeNumbers[i] = num_isotopes;
        atomCounts[i] = count;
        for (int k = 0; k < num_isotopes; k++) {
            isotope_masses.push_back(masses_table[first_idx + k]);
            isotope_probabilities.push_back(elem_table_probability[first_idx + k]);
        }
    }

    return Iso(static_cast<int>(dimNumber), isotopeNumbers.data(), atomCounts.data(),
               isotope_masses.data(), isotope_probabilities.data());
}

Iso Iso::FromFASTAWithMods(const char* sequence, const UnimodTable& mods, bool use_nominal_masses, bool add_water) {
    if (strchr(sequence, '[') == nullptr)
        // Nothing for `mods` to resolve -- go straight through the plain
        // path instead of building a composition via the general
        // dimNumber-vector Iso constructor build_iso_from_composition uses
        // (heap-allocated mass/probability vectors, per-element isotope-
        // count lookups) when Iso::FromFASTA's own fixed static tables
        // already do this cheaper for the unmodified case.
        return Iso::FromFASTA(sequence, use_nominal_masses, add_water);

    ElementComposition composition =
        add_water ? parse_fasta_with_mods_full(sequence, mods) : parse_fasta_with_mods(sequence, mods);
    return build_iso_from_composition(composition, use_nominal_masses);
}

Iso Iso::FromFASTAWithMods(const char* sequence, bool use_nominal_masses, bool add_water, const char* unimod_db_path) {
    if (strchr(sequence, '[') == nullptr)
        // As above, and also skips resolving/loading unimod_db_path's table
        // entirely when there's nothing to look up.
        return Iso::FromFASTA(sequence, use_nominal_masses, add_water);

    const UnimodTable& mods = unimod_table_for_path(unimod_db_path == nullptr ? "" : unimod_db_path);
    return FromFASTAWithMods(sequence, mods, use_nominal_masses, add_water);
}

}  // namespace IsoSpec
