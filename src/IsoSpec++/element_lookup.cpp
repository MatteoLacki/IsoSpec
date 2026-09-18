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

#include "element_lookup.h"

#include <cstring>

#include "element_tables.h"

namespace IsoSpec {

namespace {

//! Direct-indexed [first letter A-Z][second letter a-z, or index 26 for "no
//! second letter"] lookup table, built once from element_tables.cpp's 292
//! rows (verified elsewhere: 87 unique symbols, zero collisions). No map, no
//! per-lookup scan.
struct SymbolIndexTable {
    int table[26][27];

    SymbolIndexTable() {
        for (int i = 0; i < 26; i++)
            for (int j = 0; j < 27; j++)
                table[i][j] = -1;

        for (int row = 0; row < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; row++) {
            const char* symbol = elem_table_symbol[row];
            size_t len = strlen(symbol);
            if (len < 1 || len > 2 || symbol[0] < 'A' || symbol[0] > 'Z')
                continue;  // not a plain [A-Z][a-z]? symbol -- none in this table today
            int first = symbol[0] - 'A';
            int second;
            if (len == 1) {
                second = 26;
            } else {
                if (symbol[1] < 'a' || symbol[1] > 'z')
                    continue;
                second = symbol[1] - 'a';
            }
            // Only the first (lowest) row of each element's contiguous
            // isotope run is a "first index" -- later rows of the same
            // element would just overwrite with the same value if we let
            // them, so no need to special-case that.
            if (table[first][second] == -1)
                table[first][second] = row;
        }
    }
};

const SymbolIndexTable& symbol_index_table() {
    static const SymbolIndexTable instance;
    return instance;
}

}  // namespace

int find_element_table_first_index(const char* symbol, size_t len) {
    if (len < 1 || len > 2)
        return -1;
    if (symbol[0] < 'A' || symbol[0] > 'Z')
        return -1;
    int first = symbol[0] - 'A';
    int second;
    if (len == 1) {
        second = 26;
    } else {
        if (symbol[1] < 'a' || symbol[1] > 'z')
            return -1;
        second = symbol[1] - 'a';
    }
    return symbol_index_table().table[first][second];
}

int element_isotope_count(int first_index) {
    if (first_index < 0)
        return 0;
    int elem_ID = elem_table_ID[first_index];
    int count = 0;
    int idx = first_index;
    while (idx < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES && elem_table_ID[idx] == elem_ID) {
        idx++;
        count++;
    }
    return count;
}

}  // namespace IsoSpec
