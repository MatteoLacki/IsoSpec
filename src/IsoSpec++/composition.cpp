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

#include "composition.h"

#include <stdexcept>
#include <string>
#include <vector>

#include "element_lookup.h"
#include "element_tables.h"

namespace IsoSpec {

void expand_composition_into(const ElementComposition& composition, ExpandedComposition& out,
                             bool use_nominal_masses) {
    out.clear();

    const size_t dimNumber = composition.element_first_index.size();
    const double* masses_table = use_nominal_masses ? elem_table_massNo : elem_table_mass;

    out.isotopeNumbers.reserve(dimNumber);
    out.atomCounts.reserve(dimNumber);

    for (size_t i = 0; i < dimNumber; i++) {
        const int first_idx = composition.element_first_index[i];
        const int num_isotopes = element_isotope_count(first_idx);
        out.isotopeNumbers.push_back(num_isotopes);
        out.atomCounts.push_back(composition.count[i]);
        for (int k = 0; k < num_isotopes; k++) {
            out.isotope_masses.push_back(masses_table[first_idx + k]);
            out.isotope_probabilities.push_back(elem_table_probability[first_idx + k]);
        }
    }
}

void reject_negative_counts(const ElementComposition& composition) {
    std::string offenders;
    for (size_t i = 0; i < composition.count.size(); i++) {
        if (composition.count[i] < 0) {
            if (!offenders.empty())
                offenders += ", ";
            offenders += elem_table_symbol[composition.element_first_index[i]];
        }
    }

    if (offenders.empty())
        return;

    throw std::invalid_argument(
        "Negative atom count for element(s): " + offenders +
        " -- a molecule can't contain a negative number of atoms. Signed counts are only "
        "meaningful for a modification's composition delta (where a mod may remove more of an "
        "element than the rest of the molecule provides), never for a molecule itself.");
}

Iso build_iso_from_composition(const ElementComposition& composition, bool use_nominal_masses) {
    // Checked here rather than in expand_composition_into: a negative count is
    // a perfectly good intermediate in a composition (a modification's delta
    // can outweigh the bare sequence), and only becomes an error where the
    // composition becomes a molecule -- which is here.
    reject_negative_counts(composition);

    ExpandedComposition expanded;
    expand_composition_into(composition, expanded, use_nominal_masses);

    return Iso(static_cast<int>(expanded.atomCounts.size()), expanded.isotopeNumbers.data(),
               expanded.atomCounts.data(), expanded.isotope_masses.data(),
               expanded.isotope_probabilities.data());
}

}  // namespace IsoSpec
