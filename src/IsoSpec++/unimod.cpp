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

#include "unimod.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "isoSpec++.h"  // parse_formula_tokens -- reused for reading a composition column
#include "unimod_table_data.h"

namespace IsoSpec {

namespace {

// Splits a CSV line into exactly 4 fields (id, name, mono_mass, composition).
// Unimod names never contain a comma in the shipped table, but this stays
// defensive for a hand-edited override CSV: if splitting on ',' yields more
// than 4 pieces, the extra commas are assumed to belong to the name field
// (the only free-text one) and are rejoined into it.
std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ','))
        fields.push_back(field);
    if (fields.size() < 4)
        throw std::invalid_argument("Invalid unimod CSV row (expected 4 fields): " + line);
    if (fields.size() > 4) {
        std::string composition = fields.back();
        fields.pop_back();
        std::string mono_mass = fields.back();
        fields.pop_back();
        std::string id = fields.front();
        std::string name;
        for (size_t i = 1; i < fields.size(); i++) {
            if (i > 1)
                name += ",";
            name += fields[i];
        }
        fields = {id, name, mono_mass, composition};
    }
    return fields;
}

}  // namespace

UnimodTable parse_unimod_csv(const std::string& csv_text) {
    UnimodTable table;

    std::stringstream stream(csv_text);
    std::string line;
    bool first_line = true;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty())
            continue;
        if (first_line) {
            // header row ("id,name,mono_mass,composition")
            first_line = false;
            continue;
        }

        std::vector<std::string> fields = split_csv_line(line);

        unsigned long id;
        try {
            id = std::stoul(fields[0]);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid unimod CSV row (bad id): " + line);
        }

        double mono_mass;
        try {
            mono_mass = std::stod(fields[2]);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid unimod CSV row (bad mono_mass): " + line);
        }

        UnimodEntry entry;
        entry.present = true;
        entry.name = fields[1];
        entry.mono_mass = mono_mass;
        // Reuses the exact same formula-string parser as everything else in
        // IsoSpec (isoSpec++.cpp) -- the composition column is written in
        // that native grammar by scripts/build_unimod_table.py precisely so
        // this call works, rather than a second, bespoke tokenizer for
        // Unimod's own "H(3) C(2) N O" notation.
        parse_formula_tokens(fields[3].c_str(), entry.element_first_index, entry.element_delta_count);

        if (id >= table.entries_.size())
            table.entries_.resize(id + 1);
        table.entries_[id] = std::move(entry);
    }

    return table;
}

const UnimodTable& embedded_unimod_table() {
    static const UnimodTable instance = parse_unimod_csv(kEmbeddedUnimodCsv);
    return instance;
}

const UnimodTable& unimod_table_for_path(const std::string& path) {
    // Single cached (path, table) pair, not a map -- deliberately: a process
    // realistically uses 0 or 1 distinct override paths in its lifetime (set
    // once at startup), so a full path-keyed cache would buy nothing. Not
    // synchronized -- IsoSpec is single-threaded by design (see this repo's
    // CLAUDE.md); do not call this from more than one thread.
    static std::string cached_path;
    static UnimodTable cached_table;
    static bool cached_valid = false;

    if (cached_valid && cached_path == path)
        return cached_table;

    std::ifstream file(path);
    if (!file)
        throw std::invalid_argument("Failed to open unimod db: " + path);
    std::stringstream buffer;
    buffer << file.rdbuf();

    cached_table = parse_unimod_csv(buffer.str());
    cached_path = path;
    cached_valid = true;
    return cached_table;
}

}  // namespace IsoSpec
