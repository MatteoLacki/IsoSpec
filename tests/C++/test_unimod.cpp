// [UNIMOD:<id>] modification-aware peptide sequence parsing (fasta_mods.h,
// unimod.h). See docs/ai/unimod.md for the notation and the shipped table's
// exclusions.

#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "doctest.h"
#include "element_lookup.h"
#include "fasta.h"
#include "cwrapper.h"
#include "fasta_mods.h"
#include "isoSpec++.h"
#include "test_helpers.h"
#include "unimod.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {
const UnimodTable& shipped_mods() {
    static const UnimodTable table = [] {
        std::ifstream f("../../data/unimod.csv");
        if (!f)
            throw std::runtime_error("Cannot open test Unimod CSV");
        std::stringstream buffer;
        buffer << f.rdbuf();
        return parse_unimod_csv(buffer.str());
    }();
    return table;
}
}

TEST_CASE("Unimod CSV is the support ledger with original IDs") {
    const auto& table = shipped_mods();
    CHECK(table.supports(4));
    CHECK_FALSE(table.supports(9));
    CHECK_FALSE(table.supports(0));
    CHECK_FALSE(table.supports(14));
    CHECK_FALSE(table.supports(4294967300ULL));
    CHECK_FALSE(table.supports(std::numeric_limits<std::uint64_t>::max()));
    CHECK_FALSE(table.supports(table.size()));
    size_t supported = 0;
    for (size_t id = 0; id < table.size(); ++id)
        supported += table.supports(id);
    CHECK(supported == 979);
}

TEST_CASE("CSV overrides accept four columns and unsupported rows with reasons") {
    const auto table = parse_unimod_csv(
        "id,name,mono_mass,composition,reason\n"
        "4,\"Name, with comma\",1.0,H1,\n"
        "9,Unsupported,,,Isotope label\n");
    REQUIRE(table.lookup(4) != nullptr);
    CHECK(table.lookup(4)->name == "Name, with comma");
    CHECK(table.lookup(4)->element_delta_count == std::vector<int>{1});
    CHECK_FALSE(table.supports(9));
    CHECK(parse_unimod_csv("id,name,mono_mass,composition\n4,Override,1.0,H1\n").supports(4));
    CHECK_FALSE(parse_unimod_csv("id,name,mono_mass,composition,reason\n9,,,,Excluded\n").supports(9));
}

TEST_CASE("CSV rejects invalid IDs before indexing the table") {
    for (const char* id : {"-1", "4294967300", "4junk", ""}) {
        CHECK_THROWS_AS(parse_unimod_csv(std::string("id,name,mono_mass,composition\n") +
                        id + ",Bad,1.0,H1\n"), std::invalid_argument);
    }
}

TEST_CASE("FromFASTAWithMods matches FromFASTA for a plain, unmodified sequence") {
    for (const char* seq : {"PEPTIDE", "MKWVTFISLLLLFSSAYSRGV", ""}) {
        INFO("sequence='" << seq << "'");
        Iso with_mods = Iso::FromFASTAWithMods(seq, shipped_mods());
        Iso plain = Iso::FromFASTA(seq);
        CHECK(with_mods.getMonoisotopicPeakMass() == doctest::Approx(plain.getMonoisotopicPeakMass()));
        CHECK(with_mods.getTheoreticalAverageMass() == doctest::Approx(plain.getTheoreticalAverageMass()));

        Iso with_mods_dry = Iso::FromFASTAWithMods(seq, shipped_mods(), false, false);
        Iso plain_dry = Iso::FromFASTA(seq, false, false);
        CHECK(with_mods_dry.getMonoisotopicPeakMass() == doctest::Approx(plain_dry.getMonoisotopicPeakMass()));
    }
}

TEST_CASE("internal UNIMOD mod applies its composition delta") {
    // Carbamidomethyl, UNIMOD:4 -> H3C2N1O1, 57.021464 Da.
    Iso modified = Iso::FromFASTAWithMods("MPEPTC[UNIMOD:4]DEK", shipped_mods());
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
    Iso modified = Iso::FromFASTAWithMods("[UNIMOD:1]-PEPTIDE", shipped_mods());
    Iso base = Iso::FromFASTA("PEPTIDE");
    CHECK(modified.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass() ==
          doctest::Approx(42.010565).epsilon(1e-5));
}

TEST_CASE("C-terminal UNIMOD mod applies after the last residue") {
    // Oxidation, UNIMOD:35 -> O1, 15.994915 Da.
    Iso modified = Iso::FromFASTAWithMods("PEPTIDE-[UNIMOD:35]", shipped_mods());
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
    Iso modified = Iso::FromFASTAWithMods("PEPTC[UNIMOD:123]DEK", shipped_mods());
    Iso base = Iso::FromFASTA("PEPTCDEK");
    Iso expected_delta("H20C15N1O6Cl1");
    double observed = modified.getMonoisotopicPeakMass() - base.getMonoisotopicPeakMass();
    // expected_delta is itself a hydrated molecule (has its own implicit
    // terminal water via the generic formula constructor's defaults -- no,
    // Iso(formula) has no water), so compare against its raw monoisotopic
    // mass directly.
    CHECK(observed == doctest::Approx(expected_delta.getMonoisotopicPeakMass()).epsilon(1e-5));
}

TEST_CASE("unknown UNIMOD id throws") {
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:999999999]DEK", shipped_mods()), std::invalid_argument);
}

TEST_CASE("deliberately excluded UNIMOD id (isotope-labeled) throws, same as unknown") {
    // UNIMOD:9, ICAT-G:2H(8) -- isotope-labeled, excluded from the shipped
    // table by scripts/build_unimod_table.py (see docs/ai/unimod.md).
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:9]DEK", shipped_mods()), std::invalid_argument);
}

TEST_CASE("malformed UNIMOD brackets throw") {
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:4DEK", shipped_mods()), std::invalid_argument);         // missing ']'
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:]DEK", shipped_mods()), std::invalid_argument);          // missing id
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[FOO:4]DEK", shipped_mods()), std::invalid_argument);            // wrong prefix
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMODX4]DEK", shipped_mods()), std::invalid_argument);         // no colon
}

TEST_CASE("non-bracket characters keep the old lenient (silently-ignored) behavior") {
    // The fix targets '[' specifically -- everything parse_fasta already
    // tolerated (spacers, whitespace, indeterminate-formula codes) must keep
    // working exactly as before, including the documented "AE-DA" example.
    Iso base = Iso::FromFASTAWithMods("AEDA", shipped_mods());
    for (const char* variant : {"AE-DA", "EAXXDA*", "AE DA", "ae\tda", "A E D A\n", "AEDA?!"}) {
        INFO("variant='" << variant << "'");
        Iso v = Iso::FromFASTAWithMods(variant, shipped_mods());
        CHECK(v.getMonoisotopicPeakMass() == doctest::Approx(base.getMonoisotopicPeakMass()));
    }
}

TEST_CASE("unimod_db_path override replaces the packaged table") {
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

    // And the shipped table is untouched by having loaded an override.
    Iso still_default = Iso::FromFASTAWithMods("PEPTC[UNIMOD:4]DEK", shipped_mods());
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

TEST_CASE("annotated C++ sequences require an explicit CSV path or table") {
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("A[UNIMOD:4]"), std::invalid_argument);
    CHECK_NOTHROW(Iso::FromFASTAWithMods("PEPTIDE"));
    CHECK_THROWS_AS(Iso::FromFASTAWithMods("AC[UNIMOD:4294967300]", shipped_mods()), std::invalid_argument);
}

TEST_CASE("the C ABI applies the same null-path rule as the C++ overload") {
    // Both entry points must agree on what a null/empty unimod_db_path means:
    // fine for a sequence with nothing to resolve, an error for an annotated
    // one. They previously disagreed -- the C ABI substituted a silently-empty
    // table, so an annotated sequence failed with a misleading "unknown id"
    // for a bracket that could never have resolved in the first place.
    for (const char* path : {static_cast<const char*>(nullptr), ""}) {
        INFO("path=" << (path == nullptr ? "nullptr" : "\"\""));

        void* plain = parseFastaWithModsC("PEPTIDE", path);
        CHECK(plain != nullptr);
        if (plain != nullptr)
            deleteCompositionC(plain);
        CHECK_NOTHROW(Iso::FromFASTAWithMods("PEPTIDE", false, true, path));

        // c_guard turns the escaping std::invalid_argument into a NULL return.
        CHECK(parseFastaWithModsC("PEPTC[UNIMOD:4]DEK", path) == nullptr);
        CHECK_THROWS_AS(Iso::FromFASTAWithMods("PEPTC[UNIMOD:4]DEK", false, true, path),
                        std::invalid_argument);
    }
}
