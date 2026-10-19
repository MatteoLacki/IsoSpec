/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
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

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace IsoSpec{

// We will work with C H N O S Se tuples */
extern const int aa_isotope_numbers[6];

extern const double aa_elem_masses[19];

extern const double aa_elem_nominal_masses[19];

extern const double aa_elem_probabilities[19];

extern const int aa_symbol_to_elem_counts[256*6];

//! aa_symbol_to_elem_counts again, each row packed into one integer: count ii
//! in bits [10*ii, 10*ii+10). Summing packed rows sums all six counts at once,
//! as long as no field passes 1023 -- see parse_fasta().
struct PackedElemCounts { std::uint64_t row[256]; };
extern const PackedElemCounts aa_symbol_to_packed_elem_counts;

constexpr int AA_PACKED_FIELD_BITS = 10;
//! Residues summed in packed form before unpacking; fasta.cpp checks that the
//! largest count in the table times this still fits a field.
constexpr int AA_PACKED_FLUSH_EVERY = 48;

//! Count elemental composition of an unmodificed sequence of amino acids, resulting in CHNOSSe counts.
/*!
    WARNING!!! This function does not add the terminating H and OH groups, resulting in a residue backbone skeleton formula.
    If you don't know what that means, you should probably be using parse_fasta_full function.
*/
inline void parse_fasta(const char* fasta, int atomCounts[6])
{
    memset(atomCounts, 0, sizeof(decltype(atomCounts[0]))*6);

    // One load and one add per residue, using the packed table, instead of
    // six of each: the six-count row add vectorised under clang but not under
    // gcc, which made this loop -- most of the cost of reading a protein --
    // markedly slower there. Packed, both compilers emit the same short loop.
    const unsigned char* p = reinterpret_cast<const unsigned char*>(fasta);
    while(*p != '\0')
    {
        std::uint64_t packed = 0;
        for(int ii = 0; ii < AA_PACKED_FLUSH_EVERY && *p != '\0'; ++ii, ++p)
            packed += aa_symbol_to_packed_elem_counts.row[*p];
        for(int ii = 0; ii < 6; ++ii)
            atomCounts[ii] += static_cast<int>((packed >> (AA_PACKED_FIELD_BITS * ii)) & ((1u << AA_PACKED_FIELD_BITS) - 1));
    }
}

//! Turn an input FASTA aminoacid sequence into atom counts, in CHNOSSe order.
/*!
    Unlike parse_fasta, this function includes the H and OH groups at the N- and C- termini of the skeleton, resulting in a formula of a full (inert) molecule.
*/
inline void parse_fasta_full(const char* fasta, int atomCounts[6])
{
    parse_fasta(fasta, atomCounts);
    // Add terminal water (H2O) for either precursor or fragment.
    const int H_INDEX = 1;  // Indexing: 0=C, 1=H, 2=N, 3=O, 4=S, 5=Se
    const int O_INDEX = 3;

    atomCounts[H_INDEX] += 2;
    atomCounts[O_INDEX] += 1;
}



}  // namespace IsoSpec
