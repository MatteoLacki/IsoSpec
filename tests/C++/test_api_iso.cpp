// The Iso class: every public accessor, every constructor, and the invariants
// that tie them together.
//
// Where possible the expected value is computed straight from the isotope
// tables (see test_helpers::element_isotopes) rather than from another IsoSpec
// call, so these are oracles and not tautologies.

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// Sum over elements of count * f(isotopes), the shape every Iso aggregate takes.
double expected_lightest_mass(const std::vector<ElementSpec>& es) {
    double m = 0.0;
    for (const ElementSpec& e : es) {
        double lightest = std::numeric_limits<double>::infinity();
        for (const Isotope& i : e.isotopes) lightest = (std::min)(lightest, i.mass);
        m += lightest * e.count;
    }
    return m;
}

double expected_heaviest_mass(const std::vector<ElementSpec>& es) {
    double m = 0.0;
    for (const ElementSpec& e : es) {
        double heaviest = 0.0;
        for (const Isotope& i : e.isotopes) heaviest = (std::max)(heaviest, i.mass);
        m += heaviest * e.count;
    }
    return m;
}

double expected_monoisotopic_mass(const std::vector<ElementSpec>& es) {
    double m = 0.0;
    for (const ElementSpec& e : es) {
        const Isotope* best = &e.isotopes[0];
        for (const Isotope& i : e.isotopes) if (i.prob > best->prob) best = &i;
        m += best->mass * e.count;
    }
    return m;
}

double expected_average_mass(const std::vector<ElementSpec>& es) {
    double m = 0.0;
    for (const ElementSpec& e : es) {
        double avg = 0.0;
        for (const Isotope& i : e.isotopes) avg += i.mass * i.prob;
        m += avg * e.count;
    }
    return m;
}

double expected_variance(const std::vector<ElementSpec>& es) {
    double v = 0.0;
    for (const ElementSpec& e : es) {
        double avg = 0.0;
        for (const Isotope& i : e.isotopes) avg += i.mass * i.prob;
        double var = 0.0;
        for (const Isotope& i : e.isotopes) var += i.prob * (i.mass - avg) * (i.mass - avg);
        v += var * e.count;
    }
    return v;
}

const char* const kFormulas[] = {"H2O1", "C6H12O6", "C1", "S1", "Sn1", "C10H20N2O3S1"};

}  // namespace

TEST_CASE("Iso dimensions match the formula") {
    for (const char* f : kFormulas) {
        INFO("formula=" << f);
        const std::vector<ElementSpec> es = elements_of(f);
        Iso iso(f);
        CHECK(iso.getDimNumber() == static_cast<int>(es.size()));
        int allDim = 0;
        for (const ElementSpec& e : es) allDim += static_cast<int>(e.isotopes.size());
        CHECK(iso.getAllDim() == allDim);
    }
}

TEST_CASE("Iso peak masses agree with the isotope tables") {
    for (const char* f : kFormulas) {
        INFO("formula=" << f);
        const std::vector<ElementSpec> es = elements_of(f);
        Iso iso(f);

        CHECK(iso.getLightestPeakMass() == doctest::Approx(expected_lightest_mass(es)).epsilon(1e-12));
        CHECK(iso.getHeaviestPeakMass() == doctest::Approx(expected_heaviest_mass(es)).epsilon(1e-12));
        CHECK(iso.getMonoisotopicPeakMass() == doctest::Approx(expected_monoisotopic_mass(es)).epsilon(1e-12));
        CHECK(iso.getTheoreticalAverageMass() == doctest::Approx(expected_average_mass(es)).epsilon(1e-12));
        CHECK(iso.variance() == doctest::Approx(expected_variance(es)).epsilon(1e-9));
        CHECK(iso.stddev() == doctest::Approx(std::sqrt(iso.variance())));

        // Ordering invariants that must hold for any molecule.
        CHECK(iso.getLightestPeakMass() <= iso.getMonoisotopicPeakMass());
        CHECK(iso.getMonoisotopicPeakMass() <= iso.getHeaviestPeakMass());
        CHECK(iso.getLightestPeakMass() <= iso.getTheoreticalAverageMass());
        CHECK(iso.getTheoreticalAverageMass() <= iso.getHeaviestPeakMass());
        CHECK(iso.variance() >= 0.0);
    }
}

TEST_CASE("Iso peak log-probabilities are ordered and bounded") {
    for (const char* f : kFormulas) {
        INFO("formula=" << f);
        Iso iso(f);
        const double mode = iso.getModeLProb();

        // The mode is the most probable configuration: nothing may beat it.
        CHECK(iso.getLightestPeakLProb() <= mode + 1e-12);
        CHECK(iso.getHeaviestPeakLProb() <= mode + 1e-12);
        CHECK(iso.getMonoisotopicPeakLProb() <= mode + 1e-12);
        CHECK(iso.getUnlikeliestPeakLProb() <= mode + 1e-12);
        // ... and the least likely peak is no likelier than any named one.
        CHECK(iso.getUnlikeliestPeakLProb() <= iso.getHeaviestPeakLProb() + 1e-12);
        CHECK(mode <= 0.0);

        CHECK(std::isfinite(iso.getModeMass()));
        CHECK(iso.getModeMass() >= iso.getLightestPeakMass());
        CHECK(iso.getModeMass() <= iso.getHeaviestPeakMass());
    }
}

TEST_CASE("monoisotopic peak of a single-element molecule is the mode") {
    // For one element there is exactly one marginal, so the mode of the whole
    // molecule is the mode of that multinomial; for an element whose commonest
    // isotope dominates (C, H, N, O) that is the all-commonest-isotope config,
    // i.e. the monoisotopic peak.
    for (const char* f : {"C10", "H20", "N5", "O7"}) {
        INFO("formula=" << f);
        Iso iso(f);
        CHECK(iso.getMonoisotopicPeakLProb() == doctest::Approx(iso.getModeLProb()));
        CHECK(iso.getMonoisotopicPeakMass() == doctest::Approx(iso.getModeMass()));
    }
}

TEST_CASE("peak signatures are valid configurations") {
    for (const char* f : kFormulas) {
        INFO("formula=" << f);
        const std::vector<ElementSpec> es = elements_of(f);
        Iso iso(f);
        const int allDim = iso.getAllDim();

        std::vector<int> sig(allDim, -1);

        auto check_signature = [&](const char* which, double expected_mass) {
            INFO("signature=" << which);
            std::size_t off = 0;
            double mass = 0.0;
            for (const ElementSpec& e : es) {
                int total = 0;
                for (std::size_t j = 0; j < e.isotopes.size(); ++j) {
                    CHECK(sig[off + j] >= 0);
                    total += sig[off + j];
                    mass += sig[off + j] * e.isotopes[j].mass;
                }
                // Every atom of the element must be assigned to some isotope.
                CHECK(total == e.count);
                off += e.isotopes.size();
            }
            CHECK(off == static_cast<std::size_t>(allDim));
            CHECK(mass == doctest::Approx(expected_mass).epsilon(1e-12));
        };

        iso.getLightestPeakSignature(sig.data());
        check_signature("lightest", iso.getLightestPeakMass());

        iso.getHeaviestPeakSignature(sig.data());
        check_signature("heaviest", iso.getHeaviestPeakMass());

        iso.getMonoisotopicPeakSignature(sig.data());
        check_signature("monoisotopic", iso.getMonoisotopicPeakMass());
    }
}

TEST_CASE("Iso from std::string equals Iso from const char*") {
    const std::string f = "C6H12O6";
    Iso a(f);
    Iso b(f.c_str());
    CHECK(a.getDimNumber() == b.getDimNumber());
    CHECK(a.getAllDim() == b.getAllDim());
    CHECK(a.getMonoisotopicPeakMass() == b.getMonoisotopicPeakMass());
    CHECK(a.getModeLProb() == b.getModeLProb());
}

TEST_CASE("nominal masses give integer-valued peaks") {
    Iso real("C6H12O6", false);
    Iso nominal("C6H12O6", true);

    const double nm = nominal.getMonoisotopicPeakMass();
    CHECK(nm == doctest::Approx(std::round(nm)));
    CHECK(nominal.getMonoisotopicPeakMass() == doctest::Approx(180.0));  // C6H12O6 nominal
    // The real mass is close to, but not equal to, the nominal one.
    CHECK(std::fabs(real.getMonoisotopicPeakMass() - nm) < 1.0);
    CHECK(real.getMonoisotopicPeakMass() != nm);
    // Probabilities are unaffected by the mass convention.
    CHECK(real.getModeLProb() == doctest::Approx(nominal.getModeLProb()));
}

TEST_CASE("raw-array Iso reproduces the formula-built one") {
    // Build C6H12O6 by hand out of the raw tables and check it matches Iso("C6H12O6").
    const std::vector<ElementSpec> es = elements_of("C6H12O6");

    std::vector<int> isotopeNumbers, atomCounts;
    std::vector<double> masses, probs;
    for (const ElementSpec& e : es) {
        isotopeNumbers.push_back(static_cast<int>(e.isotopes.size()));
        atomCounts.push_back(e.count);
        for (const Isotope& i : e.isotopes) { masses.push_back(i.mass); probs.push_back(i.prob); }
    }

    Iso manual(static_cast<int>(es.size()), isotopeNumbers.data(), atomCounts.data(),
               masses.data(), probs.data());
    Iso parsed("C6H12O6");

    CHECK(manual.getDimNumber() == parsed.getDimNumber());
    CHECK(manual.getAllDim() == parsed.getAllDim());
    CHECK(manual.getMonoisotopicPeakMass() == doctest::Approx(parsed.getMonoisotopicPeakMass()));
    CHECK(manual.getTheoreticalAverageMass() == doctest::Approx(parsed.getTheoreticalAverageMass()));
    CHECK(manual.getModeLProb() == doctest::Approx(parsed.getModeLProb()));

    // ... and so does the array-of-arrays overload.
    std::vector<const double*> mass_ptrs, prob_ptrs;
    std::size_t off = 0;
    for (const ElementSpec& e : es) {
        mass_ptrs.push_back(masses.data() + off);
        prob_ptrs.push_back(probs.data() + off);
        off += e.isotopes.size();
    }
    Iso nested(static_cast<int>(es.size()), isotopeNumbers.data(), atomCounts.data(),
               mass_ptrs.data(), prob_ptrs.data());
    CHECK(nested.getAllDim() == parsed.getAllDim());
    CHECK(nested.getMonoisotopicPeakMass() == doctest::Approx(parsed.getMonoisotopicPeakMass()));
    CHECK(nested.getModeLProb() == doctest::Approx(parsed.getModeLProb()));
}

TEST_CASE("addElement extends the molecule") {
    Iso iso("C2");
    const int dim0 = iso.getDimNumber();
    const int all0 = iso.getAllDim();
    const double mass0 = iso.getMonoisotopicPeakMass();

    const std::vector<Isotope> hydrogen = element_isotopes("H");
    REQUIRE(hydrogen.size() == 2);
    const double h_masses[2] = {hydrogen[0].mass, hydrogen[1].mass};
    const double h_probs[2] = {hydrogen[0].prob, hydrogen[1].prob};

    iso.addElement(6, 2, h_masses, h_probs);

    CHECK(iso.getDimNumber() == dim0 + 1);
    CHECK(iso.getAllDim() == all0 + 2);
    CHECK(iso.getMonoisotopicPeakMass() ==
          doctest::Approx(mass0 + 6 * hydrogen[0].mass));

    // The result must be indistinguishable from parsing C2H6 directly.
    Iso direct("C2H6");
    CHECK(iso.getMonoisotopicPeakMass() == doctest::Approx(direct.getMonoisotopicPeakMass()));
    CHECK(iso.getTheoreticalAverageMass() == doctest::Approx(direct.getTheoreticalAverageMass()));
    CHECK(iso.getModeLProb() == doctest::Approx(direct.getModeLProb()));

    // ... and it enumerates the same distribution.
    CHECK(peaks_close(envelope_peaks(FixedEnvelope::FromThreshold(std::move(iso), 0.0, true, false)),
                      enumerate_threshold_full("C2H6")));
}

TEST_CASE("addElement onto a default-constructed Iso") {
    // dimNumber == 0 is a legal starting point for incremental construction.
    Iso iso;
    CHECK(iso.getDimNumber() == 0);
    CHECK(iso.getAllDim() == 0);

    const std::vector<Isotope> carbon = element_isotopes("C");
    const double c_masses[2] = {carbon[0].mass, carbon[1].mass};
    const double c_probs[2] = {carbon[0].prob, carbon[1].prob};
    iso.addElement(4, 2, c_masses, c_probs);

    CHECK(iso.getDimNumber() == 1);
    CHECK(iso.getAllDim() == 2);
    CHECK(iso.getMonoisotopicPeakMass() == doctest::Approx(4 * carbon[0].mass));

    CHECK(peaks_close(envelope_peaks(FixedEnvelope::FromThreshold(std::move(iso), 0.0, true, false)),
                      enumerate_threshold_full("C4")));
}

TEST_CASE("copy constructor with fullcopy=true is independent of the source") {
    Iso original("C6H12O6");
    const double mass = original.getMonoisotopicPeakMass();

    Iso copy(original, true);
    CHECK(copy.getMonoisotopicPeakMass() == mass);
    CHECK(copy.getDimNumber() == original.getDimNumber());

    // Consume the copy; the original must survive and still be usable.
    FixedEnvelope from_copy = FixedEnvelope::FromThreshold(std::move(copy), 1e-6, true, false);
    CHECK(from_copy.confs_no() > 0);
    CHECK(original.getMonoisotopicPeakMass() == mass);

    FixedEnvelope from_original = FixedEnvelope::FromThreshold(original, 1e-6, true, false);
    CHECK(from_original.confs_no() == from_copy.confs_no());
}

TEST_CASE("shallow copy (fullcopy=false) shares the source's marginals") {
    Iso original("C10H20");
    {
        Iso shallow(original, false);
        // A shallow copy is disowned: destroying it must not free the marginals.
        CHECK(shallow.getMonoisotopicPeakMass() == original.getMonoisotopicPeakMass());
        CHECK(shallow.getAllDim() == original.getAllDim());
    }
    // Original still intact after the shallow copy died.
    CHECK(original.getModeLProb() == doctest::Approx(Iso("C10H20").getModeLProb()));
}

TEST_CASE("move construction transfers ownership") {
    Iso source("C10H20O2");
    const double mass = source.getMonoisotopicPeakMass();
    Iso moved(std::move(source));
    CHECK(moved.getMonoisotopicPeakMass() == mass);
    CHECK(moved.getDimNumber() == 3);
    // `source` is now disowned; it must still be destructible (checked by ASan
    // at scope exit) — we deliberately do not read from it.
}

TEST_CASE("saveMarginalLogSizeEstimates ranks marginals sensibly") {
    Iso iso("C100H2O1");
    std::vector<double> pri(iso.getDimNumber(), 0.0);

    for (double target : {0.5, 0.9, 0.99, 0.999}) {
        CAPTURE(target);
        iso.saveMarginalLogSizeEstimates(pri.data(), target);
        for (double p : pri) CHECK_FALSE(std::isnan(p));
        // C100 dominates: a hundred atoms of a two-isotope element has far more
        // reachable subisotopologues than H2 or O1.
        CHECK(pri[0] > pri[1]);
        CHECK(pri[0] > pri[2]);
    }

    // A larger target coverage cannot shrink the estimate.
    std::vector<double> small(iso.getDimNumber()), large(iso.getDimNumber());
    iso.saveMarginalLogSizeEstimates(small.data(), 0.5);
    iso.saveMarginalLogSizeEstimates(large.data(), 0.999);
    for (int i = 0; i < iso.getDimNumber(); ++i)
        CHECK(large[i] >= small[i] - 1e-12);
}

TEST_CASE("single-isotope element has a degenerate log-size estimate") {
    // Getting a -inf estimate for a marginal that cannot vary is the documented
    // behaviour of Marginal::getLogSizeEstimate.
    Iso iso("F10");  // fluorine is monoisotopic
    double pri = 0.0;
    iso.saveMarginalLogSizeEstimates(&pri, 0.99);
    CHECK(pri == -std::numeric_limits<double>::infinity());

    // ... and its distribution is a single peak.
    FixedEnvelope env = FixedEnvelope::FromThreshold(std::move(iso), 0.0, true, false);
    CHECK(env.confs_no() == 1);
    CHECK(env.probs()[0] == doctest::Approx(1.0));
}
