// Chemical-formula parsing: what is accepted, what is rejected, and what the
// accepted forms mean.
//
// parse_formula() is documented as "NOT guaranteed to be secure against
// malicious input", so these tests pin the behaviour on *well-formed-ish*
// input and on the malformed shapes the parser explicitly rejects — they do not
// assert anything about hostile input.

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// RAII for the two int arrays parse_formula hands back via out-params.
struct ParseResult {
    unsigned int dimNumber = 0;
    std::vector<double> masses, probs;
    int* isotopeNumbers = nullptr;
    int* atomCounts = nullptr;
    unsigned int confSize = 0;

    ~ParseResult() { delete[] isotopeNumbers; delete[] atomCounts; }
};

void run_parse(ParseResult& r, const char* formula, bool nominal = false) {
    r.dimNumber = parse_formula(formula, r.masses, r.probs,
                                &r.isotopeNumbers, &r.atomCounts, &r.confSize, nominal);
}

}  // namespace

TEST_CASE("parse_formula decomposes a formula into elements") {
    ParseResult r;
    run_parse(r, "C6H12O6");

    REQUIRE(r.dimNumber == 3);
    CHECK(r.confSize == 3 * sizeof(int));
    CHECK(r.atomCounts[0] == 6);
    CHECK(r.atomCounts[1] == 12);
    CHECK(r.atomCounts[2] == 6);
    CHECK(r.isotopeNumbers[0] == 2);  // C
    CHECK(r.isotopeNumbers[1] == 2);  // H
    CHECK(r.isotopeNumbers[2] == 3);  // O
    CHECK(r.masses.size() == 7);
    CHECK(r.probs.size() == 7);

    // The flattened tables must match the raw element tables, element by element.
    const std::vector<ElementSpec> es = elements_of("C6H12O6");
    std::size_t idx = 0;
    for (const ElementSpec& e : es)
        for (const Isotope& i : e.isotopes) {
            CHECK(r.masses[idx] == i.mass);
            CHECK(r.probs[idx] == i.prob);
            ++idx;
        }
    CHECK(idx == r.masses.size());
}

TEST_CASE("parse_formula: nominal masses are the mass numbers") {
    ParseResult r;
    run_parse(r, "C1H1", true);
    REQUIRE(r.masses.size() == 4);
    CHECK(r.masses[0] == doctest::Approx(12.0));
    CHECK(r.masses[1] == doctest::Approx(13.0));
    CHECK(r.masses[2] == doctest::Approx(1.0));
    CHECK(r.masses[3] == doctest::Approx(2.0));
}

TEST_CASE("parse_formula: two-letter symbols and multi-digit counts") {
    ParseResult r;
    run_parse(r, "Se2Sn10Cl100");
    REQUIRE(r.dimNumber == 3);
    CHECK(r.atomCounts[0] == 2);
    CHECK(r.atomCounts[1] == 10);
    CHECK(r.atomCounts[2] == 100);
    CHECK(r.isotopeNumbers[0] == static_cast<int>(element_isotopes("Se").size()));
    CHECK(r.isotopeNumbers[1] == static_cast<int>(element_isotopes("Sn").size()));
    CHECK(r.isotopeNumbers[2] == static_cast<int>(element_isotopes("Cl").size()));
}

TEST_CASE("a repeated element becomes two independent marginals") {
    // C1C1 is not normalized into C2 — it stays two dimensions. Both describe
    // the same molecule, so the distributions must agree.
    Iso split("C1C1");
    Iso joined("C2");
    CHECK(split.getDimNumber() == 2);
    CHECK(joined.getDimNumber() == 1);
    CHECK(split.getTheoreticalAverageMass() == doctest::Approx(joined.getTheoreticalAverageMass()));
    CHECK(split.variance() == doctest::Approx(joined.variance()));

    CHECK(peaks_close(
        merge_equal_masses(enumerate_threshold_full("C1C1")),
        merge_equal_masses(enumerate_threshold_full("C2")), 1e-12));
}

TEST_CASE("zero-count elements are allowed and contribute nothing") {
    Iso zero("C0");
    CHECK(zero.getDimNumber() == 1);
    CHECK(zero.getMonoisotopicPeakMass() == doctest::Approx(0.0));
    CHECK(zero.getTheoreticalAverageMass() == doctest::Approx(0.0));
    CHECK(zero.variance() == doctest::Approx(0.0));

    FixedEnvelope env = FixedEnvelope::FromThreshold(std::move(zero), 0.0, true, false);
    CHECK(env.confs_no() == 1);
    CHECK(env.probs()[0] == doctest::Approx(1.0));
    CHECK(env.masses()[0] == doctest::Approx(0.0));

    // A zero-count element mixed with real ones must not perturb the result.
    CHECK(peaks_close(enumerate_threshold_full("H2O1N0"), enumerate_threshold_full("H2O1")));
}

TEST_CASE("malformed formulas are rejected with invalid_argument") {
    const char* const bad[] = {
        "",           // empty
        "H2O",        // trailing element without a count
        "H",          // ditto, single element
        "2H2",        // leading digit: element name parses as empty
        "H2 O1",      // space is neither digit nor alpha
        "H2-O1",      // ditto
        "H2O1!",      // ditto
        "Xx2",        // unknown element
        "Q1",         // unknown element, single letter
        "h2o1",       // symbols are case-sensitive
    };
    for (const char* f : bad) {
        INFO("formula='" << f << "'");
        CHECK_THROWS_AS(Iso{f}, std::invalid_argument);
    }
}

TEST_CASE("a failed parse does not leave a half-built Iso behind") {
    // Constructing from a bad formula must throw cleanly (leaks and double-frees
    // here are what the ASan configuration is for).
    for (int i = 0; i < 100; ++i)
        CHECK_THROWS_AS(Iso("NotAnElement1"), std::invalid_argument);
}

TEST_CASE("element tables are internally consistent") {
    // Not a formula test as such, but this is the data every formula resolves
    // against: if it is inconsistent, every mass in the library is suspect.
    CHECK(isospec_number_of_isotopic_entries == ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES);

    int i = 0;
    int elements_seen = 0;
    while (i < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES) {
        const int id = elem_table_ID[i];
        const int atomicNo = elem_table_atomicNo[i];
        // IDs >= 1000 are the pseudo-elements (electron, missing electron,
        // protonation): real chemistry does not apply to them.
        const bool real_element = id < 1000;
        double psum = 0.0;
        int j = i;
        for (; j < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES && elem_table_ID[j] == id; ++j) {
            INFO("entry " << j << " symbol=" << std::string(elem_table_symbol[j]));
            // All entries of one element share symbol, name and atomic number.
            CHECK(std::string(elem_table_symbol[j]) == std::string(elem_table_symbol[i]));
            CHECK(std::string(elem_table_element[j]) == std::string(elem_table_element[i]));
            CHECK(elem_table_atomicNo[j] == atomicNo);

            CHECK(elem_table_probability[j] > 0.0);
            CHECK(elem_table_probability[j] <= 1.0);
            // log_probability is the log of probability.
            CHECK(elem_table_log_probability[j] ==
                  doctest::Approx(std::log(elem_table_probability[j])).epsilon(1e-12));
            if (real_element) {
                CHECK(elem_table_mass[j] > 0.0);
                // The real mass is within a nucleon of the mass number.
                CHECK(std::fabs(elem_table_mass[j] - elem_table_massNo[j]) < 1.0);
            }
            // extraNeutrons is offset by a per-element constant (uranium counts
            // from U-233, which is not tabulated), so only the *differences*
            // are pinned down: one extra neutron per unit of mass number.
            CHECK(elem_table_extraNeutrons[j] - elem_table_extraNeutrons[i] ==
                  static_cast<int>(elem_table_massNo[j] - elem_table_massNo[i]));
            psum += elem_table_probability[j];
        }
        INFO("element " << std::string(elem_table_symbol[i]));
        CHECK(psum == doctest::Approx(1.0).epsilon(1e-9));
        ++elements_seen;
        i = j;
    }
    CHECK(elements_seen > 80);  // sanity: the table covers most of the periodic table
}

TEST_CASE("the pseudo-elements behave as documented") {
    // "E" (electron), "Me" (missing electron) and "Pn" (protonation) exist so
    // that charge states can be written into a formula.
    Iso electron("E1");
    CHECK(electron.getMonoisotopicPeakMass() == doctest::Approx(0.000548579909).epsilon(1e-6));
    Iso missing("Me1");
    CHECK(missing.getMonoisotopicPeakMass() ==
          doctest::Approx(-electron.getMonoisotopicPeakMass()));

    // A protonated molecule is heavier than the neutral one by ~1 Da.
    Iso neutral("C6H12O6");
    Iso protonated("C6H12O6Pn1");
    const double d = protonated.getMonoisotopicPeakMass() - neutral.getMonoisotopicPeakMass();
    CHECK(d > 1.0);
    CHECK(d < 1.01);

    // Adding an electron shifts the mass down by an electron mass, and does not
    // change the number of peaks (E is monoisotopic).
    CHECK(enumerate_threshold_full("C2E1").size() == enumerate_threshold_full("C2").size());
}

TEST_CASE("every tabulated element parses and yields a sane distribution") {
    // Walk the whole table: each element must be usable in a formula.
    int i = 0;
    while (i < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES) {
        const int id = elem_table_ID[i];
        const std::string sym = elem_table_symbol[i];
        int j = i;
        while (j < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES && elem_table_ID[j] == id) ++j;
        const int n_isotopes = j - i;

        INFO("element=" << sym);
        const std::string formula = sym + "1";
        Iso iso(formula);
        CHECK(iso.getDimNumber() == 1);
        CHECK(iso.getAllDim() == n_isotopes);

        // One atom: the full support is exactly the isotope list.
        FixedEnvelope env = FixedEnvelope::FromThreshold(std::move(iso), 0.0, true, false);
        CHECK(env.confs_no() == static_cast<std::size_t>(n_isotopes));
        CHECK(env.get_total_prob() == doctest::Approx(1.0).epsilon(1e-9));

        i = j;
    }
}
