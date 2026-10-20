// [UNIMOD:<id>] modification-aware peptide sequence parsing (fasta_mods.h,
// unimod.h). See docs/ai/unimod.md for the notation and the shipped table's
// exclusions.

#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

#include "composition.h"
#include "doctest.h"
#include "element_lookup.h"
#include "fasta.h"
#include "fasta_mods.h"
#include "isoSpec++.h"
#include "test_helpers.h"
#include "unimod.h"
#include "unimod_support.h"

using namespace IsoSpec;
using namespace test_helpers;

TEST_CASE("Unimod support ledger preserves IDs and rejects out-of-range IDs") {
    static_assert(is_unimod_supported(4), "Carbamidomethyl must be supported");
    static_assert(!is_unimod_supported(9), "Isotope-labeled ICAT-G is unsupported");
    CHECK_FALSE(is_unimod_supported(0));
    CHECK_FALSE(is_unimod_supported(4294967300ULL));  // must not wrap to ID 4
    CHECK_FALSE(is_unimod_supported(std::numeric_limits<std::uint64_t>::max()));

    const UnimodTable& table = embedded_unimod_table();
    const size_t ledger_size = sizeof(unimod_supported) / sizeof(unimod_supported[0]);
    REQUIRE(ledger_size == table.size());
    for (size_t id = 0; id < ledger_size; ++id) {
        INFO("id=" << id);
        CHECK(is_unimod_supported(id) == (table.lookup(static_cast<unsigned int>(id)) != nullptr));
    }
    CHECK_FALSE(is_unimod_supported(ledger_size));
}

TEST_CASE("FromFASTAWithMods matches FromFASTA for a plain, unmodified sequence") {
    for (const char* seq : {"PEPTIDE", "MKWVTFISLLLLFSSAYSRGV", ""}) {
        INFO("sequence='" << seq << "'");
        Iso with_mods = Iso::FromFASTAWithMods(seq);
        Iso plain = Iso::FromFASTA(seq);
        CHECK(with_mods.getMonoisotopicPeakMass() == doctest::Approx(plain.getMonoisotopicPeakMass()));
        CHECK(with_mods.getTheoreticalAverageMass() == doctest::Approx(plain.getTheoreticalAverageMass()));

        Iso with_mods_dry = Iso::FromFASTAWithMods(seq, false, false);
        Iso plain_dry = Iso::FromFASTA(seq, false, false);
        CHECK(with_mods_dry.getMonoisotopicPeakMass() == doctest::Approx(plain_dry.getMonoisotopicPeakMass()));
    }
}

TEST_CASE("internal UNIMOD mod applies its composition delta") {
    // Carbamidomethyl, UNIMOD:4 -> H3C2N1O1, 57.021464 Da.
    Iso modified = Iso::FromFASTAWithMods("MPEPTC[UNIMOD:4]DEK");
    Iso base = Iso::FromFASTA("MPEPTCDEK");
    double delta = modified.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass();
    CHECK(delta == doctest::Approx(57.021464).epsilon(1e-5));

    // The exact regression this feature exists to fix: before, every letter
    // in "UNIMOD" (U, N, I, M, O, D) was misread as a phantom amino acid --
    // confirmed empirically against the unfixed code this session. The fix
    // must differ from the unmodified baseline by exactly Carbamidomethyl's
    // delta, not by six extra residues (which would be off by ~350+ Da).
    CHECK(delta < 100.0);
}

TEST_CASE("N-terminal UNIMOD mod applies before the first residue") {
    // Acetyl, UNIMOD:1 -> H2C2O1, 42.010565 Da.
    Iso modified = Iso::FromFASTAWithMods("[UNIMOD:1]-PEPTIDE");
    Iso base = Iso::FromFASTA("PEPTIDE");
    CHECK(modified.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass() ==
          doctest::Approx(42.010565).epsilon(1e-5));
}

TEST_CASE("C-terminal UNIMOD mod applies after the last residue") {
    // Oxidation, UNIMOD:35 -> O1, 15.994915 Da.
    Iso modified = Iso::FromFASTAWithMods("PEPTIDE-[UNIMOD:35]");
    Iso base = Iso::FromFASTA("PEPTIDE");
    CHECK(modified.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass() ==
          doctest::Approx(15.994915).epsilon(1e-5));
}

TEST_CASE("a real shipped composition mixing a 2-letter and adjacent 1-letter element parses correctly") {
    // UNIMOD:123 ICAT-H -> H20C15N1O6Cl1 -- exercises the tokenization
    // boundary a 2-letter symbol (Cl) next to a 1-letter one (C, N, O, H)
    // depends on: every element's count must be explicit (never bare "N"),
    // or "...O6Cl1" could misparse. Cross-checked against the hand-built
    // formula for the same composition.
    Iso modified = Iso::FromFASTAWithMods("PEPTC[UNIMOD:123]DEK");
    Iso base = Iso::FromFASTA("PEPTCDEK");
    Iso expected_delta("H20C15N1O6Cl1");
    double observed = modified.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass();
    // Iso(formula) adds no water of its own, so the modification's delta
    // compares directly against expected_delta's raw monoisotopic mass.
    CHECK(observed == doctest::Approx(expected_delta.getMonoisotopicPeakMass()).epsilon(1e-5));
}

TEST_CASE("a net-negative element count is returned by the parser, refused by Iso") {
    // Met->Hsl (UNIMOD:11) is H-4C-1S-1; a lone glycine has neither the
    // hydrogens nor the sulphur to pay for it. The parser's job is
    // arithmetic, so it reports the negative net (a further modification
    // could still bring it back up); the refusal belongs where a composition
    // becomes a molecule. Both halves are contract -- fasta_mods.h says so,
    // and the C ABI and R bindings rely on the first half being true.
    ElementComposition composition = parse_fasta_with_mods("G[UNIMOD:11]");
    bool saw_negative = false;
    for (size_t i = 0; i < composition.count.size(); i++)
        if (composition.count[i] < 0)
            saw_negative = true;
    CHECK(saw_negative);

    CHECK_THROWS_AS(build_iso_from_composition(composition), std::invalid_argument);
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("G[UNIMOD:11]", false, false), std::invalid_argument);
}

TEST_CASE("unknown UNIMOD id throws") {
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:999999999]DEK"), std::invalid_argument);
}

TEST_CASE("deliberately excluded UNIMOD id (isotope-labeled) throws, same as unknown") {
    // UNIMOD:9, ICAT-G:2H(8) -- isotope-labeled, excluded from the shipped
    // table by scripts/build_unimod_table.py (see docs/ai/unimod.md).
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:9]DEK"), std::invalid_argument);
}

TEST_CASE("an id past 32 bits is unknown, not silently truncated") {
    // The id is parsed as a 64-bit value but UnimodTable::lookup takes an
    // unsigned int, so without an explicit range check 2^32 + 4 would wrap to
    // 4 and quietly apply Carbamidomethyl -- a wrong answer with no
    // diagnostic, which is the worst possible outcome for a parser whose
    // whole job is resolving ids.
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:4294967300]DEK"), std::invalid_argument);
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:8589934596]DEK"), std::invalid_argument);

    // And the id it would have wrapped onto still resolves normally.
    Iso modified = Iso::FromFASTAWithMods("PEPTC[UNIMOD:4]DEK");
    Iso base = Iso::FromFASTA("PEPTCDEK");
    CHECK(modified.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass() ==
          doctest::Approx(57.021464).epsilon(1e-5));
}

TEST_CASE("an implausibly large id in a CSV is rejected, not allocated") {
    // entries_ is indexed straight by id, so an unchecked id here is a
    // request for id * sizeof(UnimodEntry) bytes.
    CHECK_THROWS_AS(parse_unimod_csv("id,name,mono_mass,composition\n"
                                     "4000000000,Typo,1.0,H1\n"),
                    std::invalid_argument);
    // The largest id the shipped table actually uses still loads.
    CHECK(embedded_unimod_table().lookup(2147) != nullptr);
}

TEST_CASE("malformed UNIMOD brackets throw") {
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:4DEK"), std::invalid_argument);         // missing ']'
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:]DEK"), std::invalid_argument);          // missing id
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[FOO:4]DEK"), std::invalid_argument);            // wrong prefix
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMODX4]DEK"), std::invalid_argument);         // no colon
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:99999999999999999999999]DEK"),
                    std::invalid_argument);                                                      // id past 64 bits
}

TEST_CASE("non-bracket characters keep the old lenient (silently-ignored) behavior") {
    // The fix targets '[' specifically -- everything parse_fasta already
    // tolerated (spacers, whitespace, indeterminate-formula codes) must keep
    // working exactly as before, including the documented "AE-DA" example.
    Iso base = Iso::FromFASTAWithMods("AEDA");
    for (const char* variant : {"AE-DA", "EAXXDA*", "AE DA", "ae\tda", "A E D A\n", "AEDA?!"}) {
        INFO("variant='" << variant << "'");
        Iso v = Iso::FromFASTAWithMods(variant);
        CHECK(v.getMonoisotopicPeakMass() == doctest::Approx(base.getMonoisotopicPeakMass()));
    }
}

TEST_CASE("unimod_db_path override replaces the embedded table") {
    const char* tmp_path = "build/test_unimod_override.csv";
    {
        std::ofstream f(tmp_path);
        f << "id,name,mono_mass,composition\n";
        f << "4,TestOverride,1.0,H1\n";  // deliberately NOT Carbamidomethyl's real composition
    }

    Iso overridden = Iso::FromFASTAWithMods("PEPTC[UNIMOD:4]DEK", false, true, tmp_path);
    Iso base = Iso::FromFASTA("PEPTCDEK");
    double delta = overridden.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass();
    Iso expected_delta("H1");
    CHECK(delta == doctest::Approx(expected_delta.getMonoisotopicPeakMass()).epsilon(1e-5));

    // And the embedded default is untouched by having loaded an override.
    Iso still_default = Iso::FromFASTAWithMods("PEPTC[UNIMOD:4]DEK");
    double default_delta = still_default.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass();
    CHECK(default_delta == doctest::Approx(57.021464).epsilon(1e-5));

    std::remove(tmp_path);
}

TEST_CASE("unimod_table_for_path caches by path") {
    const char* tmp_path = "build/test_unimod_cache.csv";
    {
        std::ofstream f(tmp_path);
        f << "id,name,mono_mass,composition\n";
        f << "1,X,1.0,H1\n";
    }

    // A named std::string (rather than passing tmp_path inline) avoids a
    // spurious -Wdangling-reference: the returned reference is to a
    // function-static, not to anything owned by the argument, but the
    // compiler's heuristic can't see that through an implicit conversion in
    // the same expression as the reference binding.
    const std::string path(tmp_path);
    const UnimodTable& first = unimod_table_for_path(path);
    const UnimodTable& second = unimod_table_for_path(path);
    CHECK(&first == &second);

    std::remove(tmp_path);
}

TEST_CASE("embedded table agrees with data/unimod.csv on disk") {
    std::ifstream f("../../data/unimod.csv");
    REQUIRE(f.good());
    std::stringstream buffer;
    buffer << f.rdbuf();
    UnimodTable from_disk = parse_unimod_csv(buffer.str());

    for (unsigned int id : {1u, 4u, 7u, 21u, 34u, 35u, 121u, 123u}) {
        INFO("id=" << id);
        const UnimodEntry* embedded_entry = embedded_unimod_table().lookup(id);
        const UnimodEntry* disk_entry = from_disk.lookup(id);
        REQUIRE(embedded_entry != nullptr);
        REQUIRE(disk_entry != nullptr);
        CHECK(embedded_entry->name == disk_entry->name);
        CHECK(embedded_entry->mono_mass == doctest::Approx(disk_entry->mono_mass));
        CHECK(embedded_entry->element_first_index == disk_entry->element_first_index);
        CHECK(embedded_entry->element_delta_count == disk_entry->element_delta_count);
    }
}

namespace {

// Composition keyed by element, independent of the order a parser emits them.
std::map<int, int> by_element(const ElementComposition& c) {
    std::map<int, int> m;
    for (size_t i = 0; i < c.element_first_index.size(); ++i) m[c.element_first_index[i]] += c.count[i];
    return m;
}

int element(const char* symbol) {
    return find_element_table_first_index(symbol, strlen(symbol));
}

}  // namespace

TEST_CASE("override CSV: commas in a name are rejoined into the name") {
    UnimodTable table = parse_unimod_csv("id,name,mono_mass,composition\n"
                                         "7,Oxidation, or something, else,15.994915,O1\n");
    const UnimodEntry* entry = table.lookup(7);
    REQUIRE(entry != nullptr);
    CHECK(entry->name == "Oxidation, or something, else");
    CHECK(entry->mono_mass == doctest::Approx(15.994915));
    REQUIRE(entry->element_first_index.size() == 1);
    CHECK(entry->element_first_index[0] == element("O"));
    CHECK(entry->element_delta_count[0] == 1);
}

TEST_CASE("override CSV: malformed rows are rejected") {
    const std::string header = "id,name,mono_mass,composition\n";
    for (const char* row : {"7,Short,1.0\n",          // three fields
                            "x7,BadId,1.0,H1\n",      // id is not a number
                            "7,BadMass,heavy,H1\n",   // mass is not a number
                            "1000001,PastCap,1.0,H1\n"}) {  // one past the id cap
        INFO("row=" << row);
        CHECK_THROWS_AS(parse_unimod_csv(header + row), std::invalid_argument);
    }
    // The cap itself still loads.
    CHECK(parse_unimod_csv(header + "1000000,AtCap,1.0,H1\n").lookup(1000000) != nullptr);
}

TEST_CASE("unimod_table_for_path: a missing file is an error") {
    const std::string path("build/no_such_unimod_override.csv");
    std::remove(path.c_str());
    CHECK_THROWS_AS(unimod_table_for_path(path), std::invalid_argument);
}

TEST_CASE("the bracket-free fast path agrees with the general per-element parse") {
    // Without a '[' the parser counts residues with parse_fasta; with one it
    // accumulates per element. Adding a C-terminal Carbamidomethyl must change
    // the result by exactly that modification's delta.
    const UnimodEntry* carbamidomethyl = embedded_unimod_table().lookup(4);
    REQUIRE(carbamidomethyl != nullptr);
    ElementComposition reused;
    for (const char* seq : {"PEPTIDE", "MKWVTFISLLLLFSSAYSRGV", "U", "ACDEFGHIKLMNPQRSTVWY"}) {
        INFO("sequence=" << seq);
        std::map<int, int> expected = by_element(parse_fasta_with_mods(seq));
        for (size_t i = 0; i < carbamidomethyl->element_first_index.size(); ++i)
            expected[carbamidomethyl->element_first_index[i]] += carbamidomethyl->element_delta_count[i];

        parse_fasta_with_mods_into((std::string(seq) + "-[UNIMOD:4]").c_str(), reused);
        CHECK(by_element(reused) == expected);

        // And the table overload of FromFASTAWithMods short-circuits to FromFASTA.
        CHECK(Iso::FromFASTAWithMods(seq, embedded_unimod_table()).getMonoisotopicPeakMass() ==
              doctest::Approx(Iso::FromFASTA(seq).getMonoisotopicPeakMass()));
    }
}

TEST_CASE("adding water creates H and O when the modified composition has none") {
    // A modification can cancel every hydrogen (or oxygen) the residues
    // brought; the terminal water must then appear as a new H2 / O1 entry
    // rather than be added to an element that is no longer listed.
    UnimodTable table = parse_unimod_csv("id,name,mono_mass,composition\n"
                                         "1,StripH,-3.0,H-3\n"
                                         "2,Sulfur,31.972071,S1\n");
    // Glycine residue is C2H3N1O1.
    std::map<int, int> g = by_element(parse_fasta_with_mods_full("G[UNIMOD:1]", table));
    CHECK(g == std::map<int, int>{{element("C"), 2}, {element("H"), 2}, {element("N"), 1}, {element("O"), 2}});

    std::map<int, int> s = by_element(parse_fasta_with_mods_full("[UNIMOD:2]", table));
    CHECK(s == std::map<int, int>{{element("S"), 1}, {element("H"), 2}, {element("O"), 1}});
}
