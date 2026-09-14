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

#include <cstddef>

namespace IsoSpec {

//! Row index of `symbol`'s first isotope in element_tables.h's flat
//! elem_table_* arrays (grouped by element, contiguous), or -1 if `symbol`
//! isn't a known element. `len` must be 1 or 2 -- IsoSpec's own
//! element_tables.cpp never uses longer symbols. A direct-indexed
//! [26][27] table (first letter x optional second letter), built once,
//! backs this -- no scan, no map.
int find_element_table_first_index(const char* symbol, size_t len);

//! How many contiguous isotope rows starting at `first_index` belong to
//! that element (elem_table_ID stays constant over that run). Returns 0 for
//! a negative `first_index`.
int element_isotope_count(int first_index);

}  // namespace IsoSpec
