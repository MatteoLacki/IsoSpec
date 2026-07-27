// Property-based checks against an independent oracle.
//
// The oracle is brute-force convolution (test_helpers::brute_force_distribution):
// start from a single zero-mass peak and convolve in one atom at a time, using
// only the raw isotope tables.  It shares no algorithm with IsoSpec, so an
// agreement between the two is real evidence.  It is exponential in the atom
// count, hence the tiny molecules.
//
// The randomized cases use a fixed seed: a failure is reproducible, and the
// suite does not become flaky.

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// Formulas small enough for brute-force convolution (product of isotope counts
// over all atoms stays manageable).
const char* const kTiny[] = {
    "H2O1", "C2", "C3", "H4", "N2O1", "C2H4", "S1", "S2", "Cl2", "Br2",
    "C1H2N1O1", "Se1", "Sn1", "C4", "O3", "K2", "Mg2", "Si2",
};

// Number of configurations in the full support: the product over elements of
// multiset coefficients C(atoms + isotopes - 1, isotopes - 1).
double support_size(const std::vector<ElementSpec>& es) {
    double total = 1.0;
    for (const ElementSpec& e : es) {
        const int k = static_cast<int>(e.isotopes.size()) - 1;
        double c = 1.0;
        for (int i = 1; i <= k; ++i) c = c * (e.count + i) / i;
        total *= c;
    }
    return total;
}

// Random formulas over a handful of elements, with a fixed seed.
std::vector<std::string> random_formulas(unsigned seed, int count) {
    const char* const elements[] = {"C", "H", "N", "O", "S", "Cl", "Fe", "Si", "B", "Li"};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> n_elements(1, 4);
    std::uniform_int_distribution<int> which(0, 9);
    std::uniform_int_distribution<int> atoms(0, 12);

    std::vector<std::string> out;
    for (int i = 0; i < count; ++i) {
        std::string f;
        const int k = n_elements(rng);
        for (int j = 0; j < k; ++j)
            f += std::string(elements[which(rng)]) + std::to_string(atoms(rng));
        out.push_back(f);
    }
    return out;
}

}  // namespace

TEST_CASE("full enumeration equals brute-force convolution") {
    for (const char* f : kTiny) {
        INFO("formula=" << f);
        std::vector<Peak> oracle = merge_equal_masses(brute_force_distribution(elements_of(f)));
        std::vector<Peak> got = merge_equal_masses(enumerate_threshold_full(f));

        REQUIRE(got.size() == oracle.size());
        bool ok = true;
        for (std::size_t i = 0; i < got.size(); ++i)
            ok = ok && std::fabs(got[i].mass - oracle[i].mass) < 1e-9 &&
                 std::fabs(got[i].prob - oracle[i].prob) < 1e-12;
        CHECK(ok);
        CHECK(total_prob(oracle) == doctest::Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("theoretical moments equal the brute-force moments") {
    for (const char* f : kTiny) {
        INFO("formula=" << f);
        std::vector<Peak> oracle = brute_force_distribution(elements_of(f));
        Iso iso(f);

        const double mean = mean_mass(oracle);
        double var = 0.0;
        for (const Peak& p : oracle) var += p.prob * (p.mass - mean) * (p.mass - mean);

        CHECK(iso.getTheoreticalAverageMass() == doctest::Approx(mean).epsilon(1e-9));
        CHECK(iso.variance() == doctest::Approx(var).epsilon(1e-6));

        // ... and so do the mode, lightest and heaviest peaks.
        double best_prob = -1.0, lightest = 1e300, heaviest = -1e300;
        for (const Peak& p : merge_equal_masses(oracle)) {
            lightest = (std::min)(lightest, p.mass);
            heaviest = (std::max)(heaviest, p.mass);
        }
        for (const Peak& p : enumerate_threshold_full(f)) best_prob = (std::max)(best_prob, p.prob);
        CHECK(iso.getLightestPeakMass() == doctest::Approx(lightest).epsilon(1e-9));
        CHECK(iso.getHeaviestPeakMass() == doctest::Approx(heaviest).epsilon(1e-9));
        CHECK(std::exp(iso.getModeLProb()) == doctest::Approx(best_prob).epsilon(1e-9));
    }
}

TEST_CASE("randomized formulas: all three generators agree") {
    int checked = 0;
    for (const std::string& f : random_formulas(20260727, 60)) {
        INFO("formula=" << f);
        if (support_size(elements_of(f)) > 30000) continue;  // keep the suite quick
        ++checked;

        std::vector<Peak> thr = enumerate_threshold_full(f.c_str());
        std::vector<Peak> ord = enumerate_ordered_full(f.c_str());
        std::vector<Peak> lay = enumerate_layered_full(f.c_str());

        REQUIRE(thr.size() > 0);
        CHECK(thr.size() == ord.size());
        CHECK(thr.size() == lay.size());

        // Compared in the mass domain, merging near-degenerate peaks: the
        // generators sum the same terms in different orders, so masses can
        // differ in the last bits — enough to reorder two isotopologues whose
        // masses coincide to ~1e-13, which a peak-by-peak comparison would
        // then mismatch.
        CHECK(peaks_close(merge_equal_masses(thr, 1e-9), merge_equal_masses(ord, 1e-9), 1e-9));
        CHECK(peaks_close(merge_equal_masses(thr, 1e-9), merge_equal_masses(lay, 1e-9), 1e-9));
        CHECK(total_prob(thr) == doctest::Approx(1.0).epsilon(1e-9));
    }
    CHECK(checked > 20);
}

TEST_CASE("randomized formulas: enumeration matches brute-force convolution") {
    for (const std::string& f : random_formulas(31337, 25)) {
        INFO("formula=" << f);
        const std::vector<ElementSpec> es = elements_of(f);

        // Skip anything whose brute-force expansion would be enormous.
        double work = 1.0;
        for (const ElementSpec& e : es) work *= std::pow(e.isotopes.size(), e.count);
        if (work > 2e6) continue;

        std::vector<Peak> oracle = merge_equal_masses(brute_force_distribution(es));
        std::vector<Peak> got = merge_equal_masses(enumerate_threshold_full(f.c_str()));

        REQUIRE(got.size() == oracle.size());
        bool ok = true;
        for (std::size_t i = 0; i < got.size(); ++i)
            ok = ok && std::fabs(got[i].mass - oracle[i].mass) < 1e-9 &&
                 std::fabs(got[i].prob - oracle[i].prob) < 1e-12;
        CHECK(ok);
    }
}

TEST_CASE("lowering the threshold only adds peaks") {
    for (const char* f : {"C6H12O6", "C20H40N4O6S2", "Sn2Cl4"}) {
        INFO("formula=" << f);
        std::vector<Peak> previous;
        for (double thr : {1e-2, 1e-4, 1e-6, 1e-8}) {
            CAPTURE(thr);
            FixedEnvelope env = FixedEnvelope::FromThreshold(Iso(f), thr, true, false);
            std::vector<Peak> current = envelope_peaks(env);
            sort_peaks(current);

            CHECK(current.size() >= previous.size());
            CHECK(env.get_total_prob() >= total_prob(previous) - 1e-12);

            // Every peak of the coarser envelope survives into the finer one.
            // Matched with a tolerance: reordering the marginals (which the
            // generator does based on the cutoff) changes the summation order,
            // so the same isotopologue can come out differing in the last bits.
            bool superset = true;
            for (const Peak& p : previous) {
                bool found = false;
                for (const Peak& q : current)
                    if (std::fabs(p.mass - q.mass) < 1e-9 &&
                        std::fabs(p.prob - q.prob) <= 1e-12 * (1.0 + p.prob)) {
                        found = true;
                        break;
                    }
                superset = superset && found;
            }
            CHECK(superset);
            previous.swap(current);
        }
    }
}

TEST_CASE("FromTotalProb returns the most probable peaks") {
    // The trimmed envelope must be an optimal set: nothing left out may be more
    // probable than something kept.
    for (const char* f : {"C6H12O6", "C10H20O2", "S4Cl2"}) {
        INFO("formula=" << f);
        std::vector<Peak> all = enumerate_threshold_full(f);
        std::sort(all.begin(), all.end(),
                  [](const Peak& a, const Peak& b) { return a.prob > b.prob; });

        for (double target : {0.5, 0.9, 0.99}) {
            CAPTURE(target);
            FixedEnvelope env = FixedEnvelope::FromTotalProb(Iso(f), target, true, false);
            REQUIRE(env.confs_no() > 0);
            REQUIRE(env.confs_no() <= all.size());

            double smallest_kept = 1.0;
            for (std::size_t i = 0; i < env.confs_no(); ++i)
                smallest_kept = (std::min)(smallest_kept, env.prob(i));

            // The n-th most probable peak overall, where n is the envelope size:
            // everything outside must be no more probable than what is inside.
            const double largest_dropped =
                env.confs_no() < all.size() ? all[env.confs_no()].prob : 0.0;
            CHECK(largest_dropped <= smallest_kept * (1.0 + 1e-9));
            CHECK(env.get_total_prob() >= target - 1e-9);
        }
    }
}

TEST_CASE("convolving sub-molecules reproduces the whole") {
    // IsoSpec's distribution for A+B must equal the convolution of its
    // distributions for A and for B.
    struct Split { const char* a; const char* b; const char* whole; };
    const Split splits[] = {
        {"C2", "H6", "C2H6"},
        {"C3H4", "O2", "C3H4O2"},
        {"S1", "Cl2", "S1Cl2"},
        {"N2", "O1", "N2O1"},
    };
    for (const Split& s : splits) {
        INFO("whole=" << s.whole);
        FixedEnvelope a = FixedEnvelope::FromThreshold(Iso(s.a), 0.0, true, false);
        FixedEnvelope b = FixedEnvelope::FromThreshold(Iso(s.b), 0.0, true, false);
        FixedEnvelope whole = FixedEnvelope::FromThreshold(Iso(s.whole), 0.0, true, false);

        FixedEnvelope product = a * b;
        CHECK(peaks_close(merge_equal_masses(envelope_peaks(product)),
                          merge_equal_masses(envelope_peaks(whole)), 1e-12));
        CHECK(product.get_total_prob() == doctest::Approx(whole.get_total_prob()).epsilon(1e-9));
    }
}

TEST_CASE("adding an element multiplies the number of configurations") {
    // A monoisotopic element adds mass but no configurations; a polyisotopic one
    // multiplies them (the marginals are independent).
    const std::size_t base = enumerate_threshold_full("C3").size();
    CHECK(enumerate_threshold_full("C3F5").size() == base);
    CHECK(enumerate_threshold_full("C3H2").size() == base * enumerate_threshold_full("H2").size());
    CHECK(enumerate_threshold_full("C3S1").size() == base * enumerate_threshold_full("S1").size());
}

TEST_CASE("mass and probability scale correctly with molecule size") {
    // Means add, variances add: n copies of an element have n times the mean and
    // n times the variance of one.
    for (const char* element : {"C", "S", "Cl", "Sn"}) {
        INFO("element=" << element);
        Iso one(std::string(element) + "1");
        for (int n : {2, 5, 17}) {
            CAPTURE(n);
            Iso many(std::string(element) + std::to_string(n));
            CHECK(many.getTheoreticalAverageMass() ==
                  doctest::Approx(n * one.getTheoreticalAverageMass()).epsilon(1e-12));
            CHECK(many.variance() == doctest::Approx(n * one.variance()).epsilon(1e-9));
            CHECK(many.getLightestPeakMass() ==
                  doctest::Approx(n * one.getLightestPeakMass()).epsilon(1e-12));
            CHECK(many.getHeaviestPeakMass() ==
                  doctest::Approx(n * one.getHeaviestPeakMass()).epsilon(1e-12));
        }
    }
}

TEST_CASE("binning is a partition of the spectrum") {
    // Every peak lands in exactly one bin, and the bin it lands in is the one
    // containing its mass.
    for (const char* f : {"C6H12O6", "C20H40N4O6S2"}) {
        INFO("formula=" << f);
        for (double width : {0.05, 0.5, 1.0}) {
            CAPTURE(width);
            FixedEnvelope env = FixedEnvelope::FromThreshold(Iso(f), 1e-8, true, false);
            std::vector<Peak> peaks = envelope_peaks(env);
            FixedEnvelope binned = env.bin(width, 0.0);

            CHECK(binned.get_total_prob() == doctest::Approx(total_prob(peaks)).epsilon(1e-12));

            // Bin middles are distinct and increasing.
            bool increasing = true;
            for (std::size_t i = 1; i < binned.confs_no(); ++i)
                increasing = increasing && binned.mass(i - 1) < binned.mass(i);
            CHECK(increasing);

            // Every peak is within half a bin width of some bin middle.
            bool covered = true;
            for (const Peak& p : peaks) {
                bool found = false;
                for (std::size_t i = 0; i < binned.confs_no() && !found; ++i)
                    found = std::fabs(p.mass - binned.mass(i)) <= width * 0.5 + 1e-9;
                covered = covered && found;
            }
            CHECK(covered);
        }
    }
}
