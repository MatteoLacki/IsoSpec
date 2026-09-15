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

#include <charconv>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "isoSpec++.h"  // parse_formula_tokens -- reused for reading a composition column

namespace IsoSpec {

namespace {

// CSV metadata can contain quoted commas and escaped quotes. Rows occupy one
// line, as emitted by the generator; composition is the fourth column.
std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
        if (ch == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field += '"';
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (ch == ',' && !quoted) {
            fields.push_back(field);
            field.clear();
        } else {
            field += ch;
        }
    }
    fields.push_back(field);
    if (quoted || (fields.size() != 4 && fields.size() != 5))
        throw std::invalid_argument("Invalid unimod CSV row: " + line);
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
            // header row ("id,name,mono_mass,composition[,reason]")
            first_line = false;
            continue;
        }

        std::vector<std::string> fields = split_csv_line(line);

        unsigned int id;
        const auto parsed = std::from_chars(fields[0].data(), fields[0].data() + fields[0].size(), id);
        if (parsed.ec != std::errc() || parsed.ptr != fields[0].data() + fields[0].size())
            throw std::invalid_argument("Invalid unimod CSV row (bad id): " + line);
        // vector::max_size() is ~2^57 here, so checking against it would let
        // a typo'd or hostile id (say 4000000000, which fits in unsigned int)
        // through to a resize() of ~400GB -- std::bad_alloc, not the clean
        // std::invalid_argument every other malformed row produces. Cap at a
        // bound comfortably above any real Unimod id instead: the 2026-09
        // snapshot's largest is 2147.
        constexpr unsigned int kMaxUnimodId = 1000000;
        if (id > kMaxUnimodId)
            throw std::invalid_argument("Invalid unimod CSV row (id out of range): " + line);
        if (id >= table.entries_.size())
            table.entries_.resize(static_cast<size_t>(id) + 1);
        if (fields[3].empty()) {
            table.entries_[id] = UnimodEntry();
            continue;
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

        table.entries_[id] = std::move(entry);
    }

    return table;
}

const UnimodTable& unimod_table_for_path(const std::string& path) {
    // Single cached (path, table) pair, not a map -- deliberately: a process
    // realistically uses 0 or 1 distinct override paths in its lifetime (set
    // once at startup), so a full path-keyed cache would buy nothing. Not
    // synchronized -- IsoSpec is single-threaded by design (see this repo's
    // CLAUDE.md); do not call this from more than one thread.
    if (path.empty())
        throw std::invalid_argument("A Unimod CSV path is required");
    static std::string cached_path;
    static UnimodTable cached_table;
    static bool cached_valid = false;

    if (cached_valid && cached_path == path)
        return cached_table;

    // The path is interpreted as UTF-8, not the platform's narrow encoding:
    // on Windows, std::filesystem::path(std::string) decodes via the active
    // code page, which mangles a non-ASCII path before the file is ever
    // opened. std::filesystem::u8path would say this more directly but is
    // deprecated in C++20 (a hard error under this repo's -Werror debug
    // build), so construct from char8_t* instead -- same result, no
    // deprecation.
    std::ifstream file(std::filesystem::path(reinterpret_cast<const char8_t*>(path.c_str())));
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
