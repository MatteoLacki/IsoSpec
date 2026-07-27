// The Marginal family (marginalTrek++.h): the per-element subisotopologue
// distributions every generator is assembled from.
//
// A marginal of one element with n atoms is exactly a multinomial, so its
// contents can be checked against the closed form.

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "doctest.h"
#include "marginalTrek++.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// Chlorine: two isotopes, comfortably unequal probabilities.
const double kClMasses[2] = {34.96885268, 36.96590259};
const double kClProbs[2] = {0.7576, 0.2424};

// Oxygen-like: three isotopes.
const double kOMasses[3] = {15.99491461956, 16.99913170, 17.9991610};
const double kOProbs[3] = {0.9975716, 0.0003804, 0.002048};

// operators.h forward-declares a PrecalculatedMarginal in the *global*
// namespace as well, so the unqualified name is ambiguous here.
using PrecalcMarginal = IsoSpec::PrecalculatedMarginal;

Marginal make_marginal(const double* masses, const double* probs, int isotopeNo, int atoms) {
    return Marginal(masses, probs, isotopeNo, atoms);
}

// MarginalTrek, PrecalculatedMarginal and LayeredMarginal all read mode_conf in
// their constructors, so the Marginal they consume must have had its mode
// computed first — IsoGenerator's constructor does exactly this for every
// marginal before specializing it.  Hand-built marginals must do it themselves.
Marginal make_specializable(const double* masses, const double* probs, int isotopeNo, int atoms) {
    Marginal m(masses, probs, isotopeNo, atoms);
    m.ensureModeConf();
    return m;
}

// log P(conf) for a multinomial with the given per-atom probabilities.
double multinomial_lprob(const std::vector<int>& conf, const double* probs, int n) {
    double lp = std::lgamma(n + 1.0);
    for (std::size_t i = 0; i < conf.size(); ++i)
        lp += conf[i] * std::log(probs[i]) - std::lgamma(conf[i] + 1.0);
    return lp;
}

}  // namespace

TEST_CASE("Marginal reports the extremes of a subisotopologue") {
    const int atoms = 7;
    Marginal m = make_marginal(kClMasses, kClProbs, 2, atoms);

    CHECK(m.get_isotopeNo() == 2);
    CHECK(m.getLightestConfMass() == doctest::Approx(atoms * kClMasses[0]));
    CHECK(m.getHeaviestConfMass() == doctest::Approx(atoms * kClMasses[1]));
    CHECK(m.getMonoisotopicConfMass() == doctest::Approx(atoms * kClMasses[0]));
    CHECK(m.getLightestAtomIndex() == 0);
    CHECK(m.getHeaviestAtomIndex() == 1);
    CHECK(m.getMonoisotopicAtomIndex() == 0);

    // All-of-one-isotope configurations: probability p^atoms.
    CHECK(m.getLightestConfLProb() == doctest::Approx(atoms * std::log(kClProbs[0])));
    CHECK(m.getHeaviestConfLProb() == doctest::Approx(atoms * std::log(kClProbs[1])));
    CHECK(m.getMonoisotopicConfLProb() == doctest::Approx(m.getLightestConfLProb()));
    // The least likely subisotopologue is all of the rarest isotope.
    CHECK(m.getSmallestLProb() == doctest::Approx(atoms * std::log(kClProbs[1])));

    // Mean and variance of a single atom, scaled by the atom count.
    const double atom_mean = kClMasses[0] * kClProbs[0] + kClMasses[1] * kClProbs[1];
    CHECK(m.getAtomAverageMass() == doctest::Approx(atom_mean));
    CHECK(m.getTheoreticalAverageMass() == doctest::Approx(atoms * atom_mean));
    double atom_var = 0.0;
    for (int i = 0; i < 2; ++i)
        atom_var += kClProbs[i] * (kClMasses[i] - atom_mean) * (kClMasses[i] - atom_mean);
    CHECK(m.variance() == doctest::Approx(atoms * atom_var));

    // get_lProbs exposes the per-atom log-probabilities.
    CHECK(m.get_lProbs()[0] == doctest::Approx(std::log(kClProbs[0])));
    CHECK(m.get_lProbs()[1] == doctest::Approx(std::log(kClProbs[1])));
}

TEST_CASE("Marginal mode is the most probable configuration") {
    for (int atoms : {1, 2, 5, 20, 100}) {
        CAPTURE(atoms);
        Marginal m = make_marginal(kClMasses, kClProbs, 2, atoms);
        m.ensureModeConf();
        const double mode_lprob = m.getModeLProb();
        CHECK(m.fastGetModeLProb() == mode_lprob);

        // Brute force over every configuration of a two-isotope element.
        double best = -std::numeric_limits<double>::infinity();
        double best_mass = 0.0;
        for (int k = 0; k <= atoms; ++k) {
            const std::vector<int> conf = {atoms - k, k};
            const double lp = multinomial_lprob(conf, kClProbs, atoms);
            if (lp > best) {
                best = lp;
                best_mass = (atoms - k) * kClMasses[0] + k * kClMasses[1];
            }
        }
        CHECK(mode_lprob == doctest::Approx(best).epsilon(1e-12));
        CHECK(m.getModeMass() == doctest::Approx(best_mass));

        // computeModeConf hands back a caller-owned configuration.
        Conf conf = m.computeModeConf();
        REQUIRE(conf != nullptr);
        CHECK(conf[0] + conf[1] == atoms);
        delete[] conf;
    }
}

TEST_CASE("Marginal copy and move constructors") {
    Marginal original = make_marginal(kOMasses, kOProbs, 3, 10);
    original.ensureModeConf();
    const double mode = original.getModeLProb();
    const double avg = original.getTheoreticalAverageMass();

    Marginal copy(original);
    CHECK(copy.getModeLProb() == mode);
    CHECK(copy.getTheoreticalAverageMass() == doctest::Approx(avg));
    CHECK(copy.get_lProbs() != original.get_lProbs());  // deep copy

    // Copying before the mode is computed is also legal.
    Marginal fresh = make_marginal(kOMasses, kOProbs, 3, 10);
    Marginal fresh_copy(fresh);
    CHECK(fresh_copy.getModeLProb() == doctest::Approx(mode));

    Marginal moved(std::move(copy));
    CHECK(moved.getModeLProb() == mode);
    CHECK(moved.getTheoreticalAverageMass() == doctest::Approx(avg));
}

TEST_CASE("Marginal rejects a zero-probability isotope") {
    const double probs[2] = {0.0, 1.0};
    CHECK_THROWS_AS(make_marginal(kClMasses, probs, 2, 5), std::invalid_argument);

    const double negative[2] = {-0.5, 1.5};
    CHECK_THROWS_AS(make_marginal(kClMasses, negative, 2, 5), std::invalid_argument);
}

TEST_CASE("Marginal rejects an impossibly large atom count") {
    // The memoized log-factorial table bounds the number of atoms of a single
    // element; exceeding it must be a clean length_error.
    CHECK_THROWS_AS(make_marginal(kClMasses, kClProbs, 2, ISOSPEC_G_FACT_TABLE_SIZE),
                    std::length_error);
}

TEST_CASE("Marginal log-size estimate is degenerate for one isotope") {
    const double mass[1] = {18.998403};
    const double prob[1] = {1.0};
    Marginal m = make_marginal(mass, prob, 1, 10);
    CHECK(m.getLogSizeEstimate(1.0) == -std::numeric_limits<double>::infinity());
    CHECK(m.getLightestConfMass() == doctest::Approx(10 * mass[0]));
    CHECK(m.getHeaviestConfMass() == doctest::Approx(10 * mass[0]));
    CHECK(m.variance() == doctest::Approx(0.0));

    // A polyisotopic marginal gives a finite, monotone estimate.
    Marginal poly = make_marginal(kClMasses, kClProbs, 2, 100);
    const double small = poly.getLogSizeEstimate(std::log(1.0));
    const double large = poly.getLogSizeEstimate(std::log(10.0));
    CHECK(std::isfinite(small));
    CHECK(large > small);
}

TEST_CASE("PrecalculatedMarginal holds every configuration above the cutoff") {
    const int atoms = 12;
    const double cutoff = std::log(1e-6);
    PrecalcMarginal pm(make_specializable(kClMasses, kClProbs, 2, atoms), cutoff, true);

    CHECK(pm.get_no_confs() > 0);
    CHECK(pm.inRange(0));
    CHECK(pm.inRange(pm.get_no_confs() - 1));
    CHECK_FALSE(pm.inRange(pm.get_no_confs()));

    // Sorted descending by log-probability, and consistent across accessors.
    for (unsigned int i = 0; i < pm.get_no_confs(); ++i) {
        REQUIRE(pm.get_lProb(i) >= cutoff);
        REQUIRE(pm.get_prob(i) == doctest::Approx(std::exp(pm.get_lProb(i))));
        const Conf c = pm.get_conf(i);
        REQUIRE(c[0] + c[1] == atoms);
        REQUIRE(pm.get_mass(i) ==
                doctest::Approx(c[0] * kClMasses[0] + c[1] * kClMasses[1]).epsilon(1e-12));
        REQUIRE(pm.get_lProb(i) ==
                doctest::Approx(multinomial_lprob({c[0], c[1]}, kClProbs, atoms)).epsilon(1e-9));
        if (i > 0) REQUIRE(pm.get_lProb(i) <= pm.get_lProb(i - 1));
        // The pointer accessors must agree with the indexed ones.
        REQUIRE(pm.get_lProbs_ptr()[i] == pm.get_lProb(i));
        REQUIRE(pm.get_masses_ptr()[i] == pm.get_mass(i));
    }

    // Exactly the configurations above the cutoff, no more.
    int expected = 0;
    for (int k = 0; k <= atoms; ++k)
        if (multinomial_lprob({atoms - k, k}, kClProbs, atoms) >= cutoff) ++expected;
    CHECK(pm.get_no_confs() == static_cast<unsigned int>(expected));

    // The mode is the first entry when sorted.
    CHECK(pm.get_lProb(0) == doctest::Approx(pm.getModeLProb()));
}

TEST_CASE("PrecalculatedMarginal without sorting holds the same set") {
    const double cutoff = std::log(1e-8);
    PrecalcMarginal sorted(make_specializable(kOMasses, kOProbs, 3, 8), cutoff, true);
    PrecalcMarginal unsorted(make_specializable(kOMasses, kOProbs, 3, 8), cutoff, false);

    REQUIRE(sorted.get_no_confs() == unsorted.get_no_confs());

    std::vector<double> a, b;
    for (unsigned int i = 0; i < sorted.get_no_confs(); ++i) {
        a.push_back(sorted.get_lProb(i));
        b.push_back(unsorted.get_lProb(i));
    }
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(a[i] == doctest::Approx(b[i]));
}

TEST_CASE("PrecalculatedMarginal with an unreachable cutoff keeps only the mode") {
    // A cutoff at the mode's own probability admits just the mode (and any tie).
    Marginal m = make_marginal(kClMasses, kClProbs, 2, 30);
    m.ensureModeConf();
    const double mode = m.getModeLProb();
    PrecalcMarginal pm(std::move(m), mode, true);
    CHECK(pm.get_no_confs() >= 1);
    for (unsigned int i = 0; i < pm.get_no_confs(); ++i)
        CHECK(pm.get_lProb(i) == doctest::Approx(mode));
}

TEST_CASE("LayeredMarginal extends monotonically") {
    LayeredMarginal lm(make_specializable(kClMasses, kClProbs, 2, 40));
    CHECK(lm.get_no_confs() == 0);

    // The guardian at index -1 must compare as more probable than anything.
    CHECK(lm.get_lProb(-1) == std::numeric_limits<double>::infinity());

    double threshold = lm.getModeLProb() - 1.0;
    unsigned int previous = 0;
    for (int step = 0; step < 6; ++step) {
        CAPTURE(step);
        const bool extended = lm.extend(threshold);
        CHECK(extended);
        CHECK(lm.get_no_confs() >= previous);
        previous = lm.get_no_confs();

        for (unsigned int i = 0; i < lm.get_no_confs(); ++i) {
            REQUIRE(lm.get_lProb(i) >= threshold);
            REQUIRE(lm.get_prob(i) == doctest::Approx(std::exp(lm.get_lProb(i))));
            const Conf c = lm.get_conf(i);
            REQUIRE(c[0] + c[1] == 40);
            REQUIRE(lm.get_mass(i) ==
                    doctest::Approx(c[0] * kClMasses[0] + c[1] * kClMasses[1]).epsilon(1e-12));
        }
        CHECK(lm.get_min_mass() <= lm.get_max_mass());
        CHECK(lm.get_min_mass() >= 40 * kClMasses[0] - 1e-9);
        CHECK(lm.get_max_mass() <= 40 * kClMasses[1] + 1e-9);
        threshold -= 5.0;
    }

    // Eventually the whole support is in, and extending further reports that
    // there is nothing left on the fringe.
    while (lm.extend(threshold)) threshold -= 20.0;
    CHECK(lm.get_no_confs() == 41);  // 40 atoms, two isotopes -> 41 configurations
}

TEST_CASE("LayeredMarginal masses and probabilities stay paired") {
    LayeredMarginal lm(make_specializable(kOMasses, kOProbs, 3, 6));
    lm.extend(lm.getModeLProb() - 30.0);
    REQUIRE(lm.get_no_confs() > 3);

    for (unsigned int i = 0; i < lm.get_no_confs(); ++i) {
        const Conf c = lm.get_conf(i);
        int total = 0;
        double mass = 0.0, lprob_terms = 0.0;
        for (int j = 0; j < 3; ++j) {
            total += c[j];
            mass += c[j] * kOMasses[j];
            lprob_terms += c[j] * std::log(kOProbs[j]) - std::lgamma(c[j] + 1.0);
        }
        REQUIRE(total == 6);
        REQUIRE(lm.get_mass(i) == doctest::Approx(mass).epsilon(1e-12));
        REQUIRE(lm.get_lProb(i) ==
                doctest::Approx(std::lgamma(7.0) + lprob_terms).epsilon(1e-9));
    }

    // conf_lprobs()/conf_masses() are the same data in bulk form.
    CHECK(lm.conf_masses().size() == lm.get_no_confs());
}

TEST_CASE("MarginalTrek walks configurations in descending probability") {
    MarginalTrek trek(make_specializable(kClMasses, kClProbs, 2, 25));

    int idx = 0;
    double last = std::numeric_limits<double>::infinity();
    while (trek.probeConfigurationIdx(idx)) {
        const double lp = trek.conf_lprobs()[idx];
        REQUIRE(lp <= last + 1e-12);
        last = lp;
        const Conf c = trek.confs()[idx];
        REQUIRE(c[0] + c[1] == 25);
        REQUIRE(trek.conf_masses()[idx] ==
                doctest::Approx(c[0] * kClMasses[0] + c[1] * kClMasses[1]).epsilon(1e-12));
        REQUIRE(lp == doctest::Approx(multinomial_lprob({c[0], c[1]}, kClProbs, 25)).epsilon(1e-9));
        ++idx;
    }

    // Two isotopes, 25 atoms: 26 configurations, no more.
    CHECK(idx == 26);
    CHECK_FALSE(trek.probeConfigurationIdx(26));
    CHECK_FALSE(trek.probeConfigurationIdx(1000));
    CHECK(trek.conf_lprobs()[0] == doctest::Approx(trek.getModeLProb()));
}

TEST_CASE("MarginalTrek on a single-isotope element yields one configuration") {
    const double mass[1] = {18.998403};
    const double prob[1] = {1.0};
    MarginalTrek trek(make_specializable(mass, prob, 1, 10));
    CHECK(trek.probeConfigurationIdx(0));
    CHECK_FALSE(trek.probeConfigurationIdx(1));
    CHECK(trek.conf_lprobs()[0] == doctest::Approx(0.0));
    CHECK(trek.conf_masses()[0] == doctest::Approx(10 * mass[0]));
}

TEST_CASE("marginals with zero atoms are a single empty configuration") {
    Marginal m = make_marginal(kClMasses, kClProbs, 2, 0);
    CHECK(m.getLightestConfMass() == doctest::Approx(0.0));
    CHECK(m.getHeaviestConfMass() == doctest::Approx(0.0));
    CHECK(m.getTheoreticalAverageMass() == doctest::Approx(0.0));
    CHECK(m.variance() == doctest::Approx(0.0));
    CHECK(m.getModeLProb() == doctest::Approx(0.0));

    PrecalcMarginal pm(make_specializable(kClMasses, kClProbs, 2, 0), -1000.0, true);
    CHECK(pm.get_no_confs() == 1);
    CHECK(pm.get_prob(0) == doctest::Approx(1.0));
    CHECK(pm.get_mass(0) == doctest::Approx(0.0));
}
