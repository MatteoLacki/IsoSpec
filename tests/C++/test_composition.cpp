// Composition expansion: turning (element, count) pairs into the flat
// per-isotope arrays Iso's general constructor takes.
//
// expand_composition_into (fasta_mods.h) is the single implementation of that
// resolution, reached from three directions: build_iso_from_composition, the
// C ABI's expandCompositionC, and — through the latter — IsoSpecPy's
// IsoParamsFromDict, which used to carry its own copy walking PeriodicTbl.
// The tests here pin the contract all three depend on.

#include <cstring>
#include <string>
#include <vector>

#include "cwrapper.h"
#include "doctest.h"
#include "element_lookup.h"
#include "element_tables.h"
#include "fasta_mods.h"
#include "isoSpec++.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// Build an ElementComposition from symbols and counts the way a caller of the
// C++ API would — via the same element_lookup.h entry point the C ABI uses.
ElementComposition make_composition(const std::vector<const char*>& symbols,
                                   const std::vector<int>& counts) {
    REQUIRE(symbols.size() == counts.size());
    ElementComposition c;
    for (size_t i = 0; i < symbols.size(); i++) {
        const int idx = find_element_table_first_index(symbols[i], strlen(symbols[i]));
        REQUIRE(idx >= 0);
        c.element_first_index.push_back(idx);
        c.count.push_back(counts[i]);
    }
    return c;
}

// A C-side expanded-composition handle that always gets deleted.
struct ExpandedHandle {
    void* p;
    explicit ExpandedHandle(void* q) : p(q) {}
    ~ExpandedHandle() { deleteExpandedCompositionC(p); }
    ExpandedHandle(const ExpandedHandle&) = delete;
    ExpandedHandle& operator=(const ExpandedHandle&) = delete;
};

}  // namespace

TEST_CASE("expand_composition_into resolves isotopes straight from the element tables") {
    // Water, spelled the way a binding would hand it over.
    ElementComposition c = make_composition({"H", "O"}, {2, 1});

    ExpandedComposition e;
    expand_composition_into(c, e);

    REQUIRE(e.isotopeNumbers.size() == 2);
    REQUIRE(e.atomCounts.size() == 2);
    CHECK(e.isotopeNumbers[0] == 2);  // 1H, 2H
    CHECK(e.isotopeNumbers[1] == 3);  // 16O, 17O, 18O
    CHECK(e.atomCounts[0] == 2);
    CHECK(e.atomCounts[1] == 1);

    // Flat arrays are element-major and exactly sum(isotopeNumbers) long.
    REQUIRE(e.isotope_masses.size() == 5);
    REQUIRE(e.isotope_probabilities.size() == 5);

    const int h = find_element_table_first_index("H", 1);
    const int o = find_element_table_first_index("O", 1);
    for (int k = 0; k < 2; k++) {
        CHECK(e.isotope_masses[k] == elem_table_mass[h + k]);
        CHECK(e.isotope_probabilities[k] == elem_table_probability[h + k]);
    }
    for (int k = 0; k < 3; k++) {
        CHECK(e.isotope_masses[2 + k] == elem_table_mass[o + k]);
        CHECK(e.isotope_probabilities[2 + k] == elem_table_probability[o + k]);
    }
}

TEST_CASE("expand_composition_into preserves the composition's element order") {
    // Same molecule, opposite order: the expansion must follow the input, not
    // sort into any canonical order. Bindings match outputs back up by index.
    ExpandedComposition ho, oh;
    expand_composition_into(make_composition({"H", "O"}, {2, 1}), ho);
    expand_composition_into(make_composition({"O", "H"}, {1, 2}), oh);

    REQUIRE(ho.isotopeNumbers.size() == 2);
    CHECK(ho.isotopeNumbers[0] == 2);
    CHECK(oh.isotopeNumbers[0] == 3);
    CHECK(ho.atomCounts[0] == 2);
    CHECK(oh.atomCounts[0] == 1);
    CHECK(ho.isotope_masses.front() == oh.isotope_masses[3]);
}

TEST_CASE("expand_composition_into honours use_nominal_masses") {
    ElementComposition c = make_composition({"C"}, {1});

    ExpandedComposition real, nominal;
    expand_composition_into(c, real, false);
    expand_composition_into(c, nominal, true);

    const int carbon = find_element_table_first_index("C", 1);
    CHECK(real.isotope_masses[0] == elem_table_mass[carbon]);
    CHECK(nominal.isotope_masses[0] == elem_table_massNo[carbon]);
    CHECK(nominal.isotope_masses[0] == 12.0);
    // Probabilities are unaffected by the mass table choice.
    CHECK(real.isotope_probabilities[0] == nominal.isotope_probabilities[0]);
}

TEST_CASE("expand_composition_into passes negative counts through untouched") {
    // A composition is not a molecule: a Unimod delta can legitimately net a
    // negative count, and the expansion is not where that gets rejected.
    ElementComposition c = make_composition({"H", "S"}, {-1, -1});

    ExpandedComposition e;
    expand_composition_into(c, e);

    REQUIRE(e.atomCounts.size() == 2);
    CHECK(e.atomCounts[0] == -1);
    CHECK(e.atomCounts[1] == -1);
    CHECK(e.isotope_masses.size() == 2u + 4u);  // H has 2 isotopes, S has 4
}

TEST_CASE("expand_composition_into clears a reused out-parameter") {
    // The reuse contract the C ABI's hot path depends on: filling the same
    // ExpandedComposition twice must give the second call's answer, not both
    // appended together.
    ExpandedComposition e;
    expand_composition_into(make_composition({"C", "H", "N"}, {6, 12, 2}), e);
    const size_t first_dim = e.atomCounts.size();
    REQUIRE(first_dim == 3);

    expand_composition_into(make_composition({"O"}, {1}), e);
    CHECK(e.atomCounts.size() == 1);
    CHECK(e.isotopeNumbers.size() == 1);
    CHECK(e.isotope_masses.size() == 3);
    CHECK(e.isotope_probabilities.size() == 3);
    CHECK(e.atomCounts[0] == 1);
}

TEST_CASE("build_iso_from_composition agrees with the equivalent formula") {
    struct Case { std::vector<const char*> symbols; std::vector<int> counts; const char* formula; };
    const Case cases[] = {
        {{"H", "O"}, {2, 1}, "H2O1"},
        {{"C", "H", "N", "O", "S"}, {6, 12, 2, 1, 1}, "C6H12N2O1S1"},
        {{"C"}, {10}, "C10"},
    };

    for (const Case& c : cases) {
        INFO("formula=" << c.formula);
        Iso from_composition = build_iso_from_composition(make_composition(c.symbols, c.counts));
        Iso from_formula(c.formula);

        CHECK(from_composition.getMonoisotopicPeakMass() == from_formula.getMonoisotopicPeakMass());
        CHECK(from_composition.getTheoreticalAverageMass() ==
              doctest::Approx(from_formula.getTheoreticalAverageMass()));
        CHECK(from_composition.getModeLProb() == doctest::Approx(from_formula.getModeLProb()));
    }
}

TEST_CASE("build_iso_from_composition still rejects a negative count") {
    // Moving the expansion out from under it must not have moved this check
    // out with it: a negative atom count reaching Iso's constructor is
    // undefined behaviour, not merely a wrong answer.
    CHECK_THROWS_AS(build_iso_from_composition(make_composition({"H"}, {-1})), std::invalid_argument);
    CHECK_THROWS_AS(build_iso_from_composition(make_composition({"C", "H"}, {6, -2})),
                    std::invalid_argument);
}

TEST_CASE("expandCompositionC feeds setupIso an Iso equal to the formula's") {
    // The round trip a binding actually performs: composition in, flat arrays
    // out, straight into setupIso.
    const char* symbols[] = {"C", "H", "N", "O", "S"};
    const int counts[] = {6, 12, 2, 1, 1};
    const size_t size = 5;

    ExpandedHandle e(expandCompositionC(symbols, counts, size, false));
    REQUIRE(e.p != nullptr);

    const int* isotopeNumbers = expandedIsotopeNumbersC(e.p);
    const int* atomCounts = expandedAtomCountsC(e.p);
    const double* masses = expandedIsotopeMassesC(e.p);
    const double* probs = expandedIsotopeProbabilitiesC(e.p);
    REQUIRE(isotopeNumbers != nullptr);
    REQUIRE(atomCounts != nullptr);
    REQUIRE(masses != nullptr);
    REQUIRE(probs != nullptr);

    for (size_t i = 0; i < size; i++)
        CHECK(atomCounts[i] == counts[i]);

    void* iso = setupIso(static_cast<int>(size), isotopeNumbers, atomCounts, masses, probs);
    REQUIRE(iso != nullptr);
    Iso reference("C6H12N2O1S1");
    CHECK(getMonoisotopicPeakMassIso(iso) == reference.getMonoisotopicPeakMass());
    CHECK(getModeLProbIso(iso) == doctest::Approx(reference.getModeLProb()));
    deleteIso(iso);
}

TEST_CASE("expandCompositionC agrees with the C++ expansion, nominal masses included") {
    const char* symbols[] = {"Se", "P"};
    const int counts[] = {2, 3};

    for (bool nominal : {false, true}) {
        INFO("use_nominal_masses=" << nominal);
        ExpandedHandle e(expandCompositionC(symbols, counts, 2, nominal));
        REQUIRE(e.p != nullptr);

        ExpandedComposition expected;
        expand_composition_into(make_composition({"Se", "P"}, {2, 3}), expected, nominal);

        const int* isotopeNumbers = expandedIsotopeNumbersC(e.p);
        for (size_t i = 0; i < expected.isotopeNumbers.size(); i++)
            CHECK(isotopeNumbers[i] == expected.isotopeNumbers[i]);

        const double* masses = expandedIsotopeMassesC(e.p);
        const double* probs = expandedIsotopeProbabilitiesC(e.p);
        for (size_t k = 0; k < expected.isotope_masses.size(); k++) {
            CHECK(masses[k] == expected.isotope_masses[k]);
            CHECK(probs[k] == expected.isotope_probabilities[k]);
        }
    }
}

TEST_CASE("expandCompositionC returns NULL for an unusable element symbol") {
    const int counts[] = {1};
    // Not an element at all; too long for any symbol; empty; lower-case first
    // letter. Each must come back as NULL rather than unwinding across the ABI.
    for (const char* bad : {"Xx", "Unobtainium", "", "c", "H2"}) {
        INFO("symbol=" << bad);
        const char* symbols[] = {bad};
        CHECK(expandCompositionC(symbols, counts, 1, false) == nullptr);
    }
}

TEST_CASE("expandCompositionC accepts a negative count without complaint") {
    // Same doctrine as the C++ entry point: the composition layer does not
    // adjudicate molecules. A binding that goes on to build an Iso from this
    // is the one that must check.
    const char* symbols[] = {"H", "S"};
    const int counts[] = {-1, -1};

    ExpandedHandle e(expandCompositionC(symbols, counts, 2, false));
    REQUIRE(e.p != nullptr);
    const int* atomCounts = expandedAtomCountsC(e.p);
    CHECK(atomCounts[0] == -1);
    CHECK(atomCounts[1] == -1);
}

TEST_CASE("expandCompositionC handles an empty composition") {
    ExpandedHandle e(expandCompositionC(nullptr, nullptr, 0, false));
    REQUIRE(e.p != nullptr);
    // Zero-length vectors: data() may be null, and the only promise is that
    // nothing is read through it. The getters themselves must not fault.
    CHECK(expandedIsotopeNumbersC(e.p) == expandedIsotopeNumbersC(e.p));
    CHECK(expandedAtomCountsC(e.p) == expandedAtomCountsC(e.p));
}

TEST_CASE("the mods-aware parser's output feeds expandCompositionC unchanged") {
    // The end-to-end path IsoSpecPy takes: parseFastaWithModsC hands back
    // symbols and counts, which go straight into expandCompositionC. The two
    // handle types must agree on the (symbol, count) vocabulary.
    void* composition = parseFastaWithModsC("PEPTC[UNIMOD:4]DEK", nullptr);
    REQUIRE(composition != nullptr);

    const size_t size = compositionSizeC(composition);
    const char* const* symbols = compositionSymbolsC(composition);
    const int* counts = compositionCountsC(composition);
    REQUIRE(size > 0);

    ExpandedHandle e(expandCompositionC(symbols, counts, size, false));
    REQUIRE(e.p != nullptr);

    const int* isotopeNumbers = expandedIsotopeNumbersC(e.p);
    const int* atomCounts = expandedAtomCountsC(e.p);
    for (size_t i = 0; i < size; i++) {
        INFO("element=" << symbols[i]);
        CHECK(atomCounts[i] == counts[i]);
        CHECK(isotopeNumbers[i] ==
              element_isotope_count(find_element_table_first_index(symbols[i], strlen(symbols[i]))));
    }

    void* iso = setupIso(static_cast<int>(size), isotopeNumbers, atomCounts,
                         expandedIsotopeMassesC(e.p), expandedIsotopeProbabilitiesC(e.p));
    REQUIRE(iso != nullptr);
    // Same molecule the C++ side builds directly, modulo the terminal water
    // parseFastaWithModsC deliberately omits.
    Iso reference = Iso::FromFASTAWithMods("PEPTC[UNIMOD:4]DEK", false, false);
    CHECK(getMonoisotopicPeakMassIso(iso) == doctest::Approx(reference.getMonoisotopicPeakMass()));
    deleteIso(iso);

    deleteCompositionC(composition);
}
