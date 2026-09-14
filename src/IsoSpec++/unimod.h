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

#include <string>
#include <vector>

namespace IsoSpec {

//! One Unimod modification's resolved, ready-to-merge atomic composition
//! delta: parallel arrays over each element it touches, `element_first_index`
//! being that element's row in element_tables.h's flat arrays (see
//! element_lookup.h) and `element_delta_count` its signed atom-count change.
struct UnimodEntry {
    bool present = false;
    std::string name;
    double mono_mass = 0.0;
    std::vector<int> element_first_index;
    std::vector<int> element_delta_count;
};

//! Unimod ids -> resolved compositions, dense: a plain `vector<UnimodEntry>`
//! indexed directly by id (sized one past the largest shipped id), not a
//! map -- `lookup(id)` is a bounds check plus a `present` flag, matching
//! Unimod's own id space (ids are assigned once and never reused).
class UnimodTable {
 public:
    size_t size() const { return entries_.size(); }

    //! nullptr for an id that is out of range, or in range but not shipped
    //! (deliberately excluded -- isotope-labeled or glycan/derivatization
    //! "brick" entries, see scripts/build_unimod_table.py) -- both cases are
    //! the same "unknown to this table" outcome, no need to distinguish them.
    const UnimodEntry* lookup(unsigned int id) const {
        if (id >= entries_.size() || !entries_[id].present)
            return nullptr;
        return &entries_[id];
    }

 private:
    std::vector<UnimodEntry> entries_;

    friend UnimodTable parse_unimod_csv(const std::string& csv_text);
};

//! Parses a `data/unimod.csv`-shaped CSV (`id,name,mono_mass,composition`,
//! composition in IsoSpec's own native formula-string grammar -- see
//! scripts/build_unimod_table.py) into a UnimodTable. Throws
//! std::invalid_argument on a malformed row or composition.
UnimodTable parse_unimod_csv(const std::string& csv_text);

//! The default table, embedded at compile time (src/IsoSpec++/
//! unimod_table_data.h) and parsed once, lazily.
const UnimodTable& embedded_unimod_table();

//! Loads and parses `path` (same CSV shape as the embedded default),
//! overriding it. `unimod_db_path == nullptr` or `""` should use
//! `embedded_unimod_table()` instead -- callers check that, not this
//! function. Caches a single (path, table) pair -- not a map: a process
//! realistically uses 0 or 1 distinct override paths, so a full path-keyed
//! map buys nothing over remembering just the last one. IsoSpec is
//! single-threaded by design (see this repo's CLAUDE.md), so this cache is
//! not synchronized; do not call it from more than one thread.
const UnimodTable& unimod_table_for_path(const std::string& path);

}  // namespace IsoSpec
