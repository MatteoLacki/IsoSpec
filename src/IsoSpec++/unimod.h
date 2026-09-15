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

#include <cstdint>
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
    //! One past the largest id in the CSV this was parsed from -- ids are
    //! row indices, gaps and unsupported entries included, so this is not a
    //! count of supported modifications.
    size_t size() const { return entries_.size(); }

    //! Whether this table can resolve `id`. Replaces the generated
    //! compile-time ledger (unimod_support.h, removed): the loaded table is
    //! now the only source of truth, so this cannot drift from it.
    bool supports(std::uint64_t id) const { return lookup(id) != nullptr; }

    //! nullptr for an id that is out of range, or in range but not shipped
    //! (deliberately excluded -- isotope-labeled or glycan/derivatization
    //! "brick" entries, see scripts/build_unimod_table.py) -- both cases are
    //! the same "unknown to this table" outcome, no need to distinguish them.
    const UnimodEntry* lookup(std::uint64_t id) const {
        if (id >= entries_.size() || !entries_[id].present)
            return nullptr;
        return &entries_[id];
    }

 private:
    std::vector<UnimodEntry> entries_;

    friend UnimodTable parse_unimod_csv(const std::string& csv_text);
};

//! Parses a `data/unimod.csv`-shaped CSV into a UnimodTable. Columns are
//! `id,name,mono_mass,composition[,reason]`, composition in IsoSpec's own
//! native formula-string grammar (see scripts/build_unimod_table.py). An
//! empty composition marks an id the shipped table deliberately does not
//! support, with `reason` recording why in human-readable form -- that is
//! how the table doubles as its own support ledger, so `supports()` answers
//! for whatever table was actually loaded rather than for something compiled
//! in. The pre-`reason` four-column shape still loads unchanged, so older
//! hand-written override CSVs keep working. Fields may be double-quoted (with
//! `""` escapes) because names and reasons contain commas; one record must
//! occupy one line. Throws std::invalid_argument on a malformed row or
//! composition.
UnimodTable parse_unimod_csv(const std::string& csv_text);

//! Loads and parses `path`. Caches a single (path, table) pair -- not a map:
//! a process realistically uses 0 or 1 distinct paths, so a full path-keyed
//! map buys nothing over remembering just the last one; loading a different
//! path replaces the cached table, and a caller that needs several tables
//! alive at once owns parse_unimod_csv() results instead. An empty path
//! throws -- there is no compile-time default table to fall back to; callers
//! that want "no modification data" for a sequence with no brackets check
//! that themselves (see Iso::FromFASTAWithMods). IsoSpec is single-threaded
//! by design (see this repo's CLAUDE.md), so this cache is not synchronized;
//! do not call it from more than one thread.
const UnimodTable& unimod_table_for_path(const std::string& path);

}  // namespace IsoSpec
