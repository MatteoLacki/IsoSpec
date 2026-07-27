// FixedEnvelope: the materialized-spectrum type.
//
// Construction (all four named constructors), value semantics, sorting,
// statistics, arithmetic, binning and resampling.  Small hand-built envelopes
// are used wherever the expected answer can be written down exactly.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// FixedEnvelope takes ownership of malloc'd arrays; build one from literals.
FixedEnvelope make_envelope(const std::vector<double>& masses,
                            const std::vector<double>& probs) {
    REQUIRE(masses.size() == probs.size());
    const std::size_t n = masses.size();
    double* m = reinterpret_cast<double*>(malloc(sizeof(double) * (n == 0 ? 1 : n)));
    double* p = reinterpret_cast<double*>(malloc(sizeof(double) * (n == 0 ? 1 : n)));
    REQUIRE(m != nullptr);
    REQUIRE(p != nullptr);
    if (n > 0) {
        memcpy(m, masses.data(), sizeof(double) * n);
        memcpy(p, probs.data(), sizeof(double) * n);
    }
    return FixedEnvelope(m, p, n);
}

bool is_sorted_by(const double* arr, std::size_t n) {
    for (std::size_t i = 1; i < n; ++i)
        if (arr[i - 1] > arr[i]) return false;
    return true;
}

}  // namespace

TEST_CASE("default-constructed envelope is empty") {
    FixedEnvelope env;
    CHECK(env.confs_no() == 0);
    CHECK(env.getAllDim() == 0);
    CHECK(env.get_total_prob() == doctest::Approx(0.0));
    // Sorting and binning an empty envelope must not crash.
    CHECK_NOTHROW(env.sort_by_mass());
    CHECK_NOTHROW(env.sort_by_prob());
    CHECK_NOTHROW(env.shift_mass(1.0));
    CHECK_NOTHROW(env.scale(2.0));
    FixedEnvelope binned = env.bin(1.0);
    CHECK(binned.confs_no() == 0);
    CHECK_THROWS_AS(env.resample(10), std::logic_error);
}

TEST_CASE("envelope built from arrays reports what it was given") {
    FixedEnvelope env = make_envelope({3.0, 1.0, 2.0}, {0.5, 0.25, 0.25});
    CHECK(env.confs_no() == 3);
    CHECK(env.mass(0) == 3.0);
    CHECK(env.prob(0) == 0.5);
    CHECK(env.masses()[1] == 1.0);
    CHECK(env.probs()[1] == 0.25);
    CHECK(env.confs() == nullptr);
    CHECK(env.get_total_prob() == doctest::Approx(1.0));
}

TEST_CASE("copy and move semantics") {
    FixedEnvelope original = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-6, true, true);
    REQUIRE(original.confs_no() > 1);
    const std::size_t n = original.confs_no();
    const double m0 = original.mass(0);
    const int allDim = original.getAllDim();

    SUBCASE("copy constructor deep-copies") {
        FixedEnvelope copy(original);
        CHECK(copy.confs_no() == n);
        CHECK(copy.getAllDim() == allDim);
        CHECK(copy.masses() != original.masses());
        CHECK(copy.probs() != original.probs());
        CHECK(copy.confs() != original.confs());
        for (std::size_t i = 0; i < n; ++i) {
            CHECK(copy.mass(i) == original.mass(i));
            CHECK(copy.prob(i) == original.prob(i));
        }
        CHECK(memcmp(copy.confs(), original.confs(), sizeof(int) * n * allDim) == 0);

        // Mutating the copy leaves the original untouched.
        copy.shift_mass(100.0);
        CHECK(original.mass(0) == m0);
    }

    SUBCASE("copy assignment deep-copies") {
        FixedEnvelope copy;
        copy = original;
        CHECK(copy.confs_no() == n);
        CHECK(copy.masses() != original.masses());
        for (std::size_t i = 0; i < n; ++i) CHECK(copy.mass(i) == original.mass(i));
    }

    SUBCASE("move constructor steals the buffers") {
        const double* buf = original.masses();
        FixedEnvelope moved(std::move(original));
        CHECK(moved.confs_no() == n);
        CHECK(moved.masses() == buf);
    }

    SUBCASE("move assignment steals the buffers") {
        const double* buf = original.masses();
        FixedEnvelope moved;
        moved = std::move(original);
        CHECK(moved.confs_no() == n);
        CHECK(moved.masses() == buf);
    }

    SUBCASE("self-assignment is harmless") {
        FixedEnvelope& alias = original;
        original = alias;
        CHECK(original.confs_no() == n);
        CHECK(original.mass(0) == m0);
    }
}

TEST_CASE("release_* hands ownership to the caller") {
    FixedEnvelope env = FixedEnvelope::FromThreshold(Iso("H2O1"), 0.0, true, true);
    const std::size_t n = env.confs_no();
    REQUIRE(n > 0);

    double* masses = env.release_masses();
    double* probs = env.release_probs();
    int* confs = env.release_confs();
    REQUIRE(masses != nullptr);
    REQUIRE(probs != nullptr);
    REQUIRE(confs != nullptr);
    CHECK(env.masses() == nullptr);
    CHECK(env.probs() == nullptr);
    CHECK(env.confs() == nullptr);
    // The data is still valid after the envelope stops owning it.
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) sum += probs[i];
    CHECK(sum == doctest::Approx(1.0));
    free(masses);
    free(probs);
    free(confs);

    FixedEnvelope env2 = FixedEnvelope::FromThreshold(Iso("H2O1"), 0.0, true, false);
    double* m2 = env2.masses() ? const_cast<double*>(env2.masses()) : nullptr;
    double* p2 = const_cast<double*>(env2.probs());
    env2.release_everything();
    CHECK(env2.masses() == nullptr);
    CHECK(env2.probs() == nullptr);
    free(m2);
    free(p2);
}

TEST_CASE("sort_by_mass and sort_by_prob order the whole envelope") {
    FixedEnvelope env = FixedEnvelope::FromThreshold(Iso("C10H20O2"), 1e-8, true, true);
    REQUIRE(env.confs_no() > 10);
    const std::size_t n = env.confs_no();
    const int allDim = env.getAllDim();

    // Peaks as a set are invariant under sorting.
    std::vector<Peak> before = envelope_peaks(env);
    const double total = env.get_total_prob();

    env.sort_by_mass();
    CHECK(is_sorted_by(env.masses(), n));
    CHECK(peaks_close(before, envelope_peaks(env), 0.0));
    CHECK(env.get_total_prob() == doctest::Approx(total));
    // Sorting twice is a no-op (the cached flag must not corrupt anything).
    env.sort_by_mass();
    CHECK(is_sorted_by(env.masses(), n));

    // Configurations travel with their peaks: the mass implied by each
    // signature must still match the stored mass.
    const std::vector<ElementSpec> es = elements_of("C10H20O2");
    for (std::size_t i = 0; i < n; ++i) {
        const int* c = env.conf(i);
        double mass = 0.0;
        int off = 0;
        for (const ElementSpec& e : es)
            for (std::size_t j = 0; j < e.isotopes.size(); ++j)
                mass += c[off++] * e.isotopes[j].mass;
        REQUIRE(off == allDim);
        REQUIRE(env.mass(i) == doctest::Approx(mass).epsilon(1e-12));
    }

    env.sort_by_prob();
    CHECK(is_sorted_by(env.probs(), n));
    CHECK(peaks_close(before, envelope_peaks(env), 0.0));

    // Back to mass order: the flags must not short-circuit this.
    env.sort_by_mass();
    CHECK(is_sorted_by(env.masses(), n));
}

TEST_CASE("get_total_prob, scale and normalize") {
    FixedEnvelope env = make_envelope({1.0, 2.0, 3.0}, {0.2, 0.3, 0.5});
    CHECK(env.get_total_prob() == doctest::Approx(1.0));

    env.scale(2.0);
    CHECK(env.get_total_prob() == doctest::Approx(2.0));
    CHECK(env.prob(0) == doctest::Approx(0.4));

    env.normalize();
    CHECK(env.get_total_prob() == doctest::Approx(1.0));
    CHECK(env.prob(0) == doctest::Approx(0.2));
    CHECK(env.prob(2) == doctest::Approx(0.5));

    // Normalizing an already-normalized envelope changes nothing.
    const double p0 = env.prob(0);
    env.normalize();
    CHECK(env.prob(0) == p0);

    // A truncated envelope normalizes to 1 without changing peak ratios.
    FixedEnvelope truncated = FixedEnvelope::FromThreshold(Iso("C100"), 1e-3, true, false);
    const double before = truncated.get_total_prob();
    CHECK(before < 1.0);
    const double ratio = truncated.prob(0) / truncated.prob(1);
    truncated.normalize();
    CHECK(truncated.get_total_prob() == doctest::Approx(1.0));
    CHECK(truncated.prob(0) / truncated.prob(1) == doctest::Approx(ratio));
}

TEST_CASE("shift_mass translates the spectrum") {
    FixedEnvelope env = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-6, true, false);
    const double avg = env.empiric_average_mass();
    const double var = env.empiric_variance();
    const double total = env.get_total_prob();

    env.shift_mass(-1000.0);
    CHECK(env.empiric_average_mass() == doctest::Approx(avg - 1000.0));
    CHECK(env.empiric_variance() == doctest::Approx(var).epsilon(1e-6));
    CHECK(env.get_total_prob() == doctest::Approx(total));

    env.shift_mass(1000.0);
    CHECK(env.empiric_average_mass() == doctest::Approx(avg));
}

TEST_CASE("empiric statistics match a hand computation") {
    // Deliberately unnormalized: the statistics must weight by probability.
    FixedEnvelope env = make_envelope({10.0, 20.0, 30.0}, {1.0, 2.0, 1.0});
    const double mean = (10.0 * 1 + 20.0 * 2 + 30.0 * 1) / 4.0;
    CHECK(env.empiric_average_mass() == doctest::Approx(mean));
    const double var = (1 * 100.0 + 2 * 0.0 + 1 * 100.0) / 4.0;
    CHECK(env.empiric_variance() == doctest::Approx(var));
    CHECK(env.empiric_stddev() == doctest::Approx(std::sqrt(var)));
}

TEST_CASE("empiric statistics of a full envelope match the theoretical ones") {
    for (const char* f : {"C10H20O2", "C6H12O6", "S4", "Sn2"}) {
        INFO("formula=" << f);
        Iso iso(f);
        FixedEnvelope env = FixedEnvelope::FromThreshold(iso, 0.0, true, false);
        CHECK(env.get_total_prob() == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(env.empiric_average_mass() ==
              doctest::Approx(iso.getTheoreticalAverageMass()).epsilon(1e-9));
        CHECK(env.empiric_variance() == doctest::Approx(iso.variance()).epsilon(1e-6));
        CHECK(env.empiric_stddev() == doctest::Approx(iso.stddev()).epsilon(1e-6));
    }
}

TEST_CASE("operator+ concatenates two spectra") {
    FixedEnvelope a = make_envelope({1.0, 2.0}, {0.3, 0.7});
    FixedEnvelope b = make_envelope({5.0}, {1.0});
    FixedEnvelope sum = a + b;

    CHECK(sum.confs_no() == 3);
    CHECK(sum.get_total_prob() == doctest::Approx(2.0));
    std::vector<Peak> expect = {{1.0, 0.3}, {2.0, 0.7}, {5.0, 1.0}};
    CHECK(peaks_close(envelope_peaks(sum), expect, 1e-12));

    // Adding an empty envelope is the identity.
    FixedEnvelope empty;
    FixedEnvelope same = a + empty;
    CHECK(same.confs_no() == a.confs_no());
    CHECK(peaks_close(envelope_peaks(same), envelope_peaks(a), 0.0));
}

TEST_CASE("operator* convolves two spectra") {
    FixedEnvelope a = make_envelope({1.0, 2.0}, {0.5, 0.5});
    FixedEnvelope b = make_envelope({10.0, 20.0}, {0.25, 0.75});
    FixedEnvelope prod = a * b;

    REQUIRE(prod.confs_no() == 4);
    std::vector<Peak> expect = {{11.0, 0.125}, {21.0, 0.375}, {12.0, 0.125}, {22.0, 0.375}};
    CHECK(peaks_close(envelope_peaks(prod), expect, 1e-12));
    CHECK(prod.get_total_prob() == doctest::Approx(1.0));
    // Mean masses add under convolution.
    CHECK(prod.empiric_average_mass() ==
          doctest::Approx(a.empiric_average_mass() + b.empiric_average_mass()));
}

TEST_CASE("convolving element spectra reproduces the molecule's spectrum") {
    // C2 * H6 must equal C2H6 — a real, non-trivial check of operator*.
    FixedEnvelope carbon = FixedEnvelope::FromThreshold(Iso("C2"), 0.0, true, false);
    FixedEnvelope hydrogen = FixedEnvelope::FromThreshold(Iso("H6"), 0.0, true, false);
    FixedEnvelope product = carbon * hydrogen;
    FixedEnvelope direct = FixedEnvelope::FromThreshold(Iso("C2H6"), 0.0, true, false);

    CHECK(product.confs_no() == direct.confs_no());
    CHECK(product.get_total_prob() == doctest::Approx(direct.get_total_prob()));
    CHECK(peaks_close(merge_equal_masses(envelope_peaks(product)),
                      merge_equal_masses(envelope_peaks(direct)), 1e-12));
}

TEST_CASE("LinearCombination mixes spectra with given intensities") {
    FixedEnvelope a = make_envelope({1.0, 2.0}, {0.5, 0.5});
    FixedEnvelope b = make_envelope({3.0}, {1.0});

    const std::vector<const FixedEnvelope*> spectra = {&a, &b};
    const std::vector<double> intensities = {2.0, 10.0};
    FixedEnvelope mix = FixedEnvelope::LinearCombination(spectra, intensities);

    REQUIRE(mix.confs_no() == 3);
    std::vector<Peak> expect = {{1.0, 1.0}, {2.0, 1.0}, {3.0, 10.0}};
    CHECK(peaks_close(envelope_peaks(mix), expect, 1e-12));
    CHECK(mix.get_total_prob() == doctest::Approx(12.0));

    // The pointer/size overload must agree.
    const FixedEnvelope* raw[2] = {&a, &b};
    FixedEnvelope mix2 = FixedEnvelope::LinearCombination(raw, intensities.data(), 2);
    CHECK(peaks_close(envelope_peaks(mix2), envelope_peaks(mix), 0.0));

    // Zero intensity contributes zero probability but still contributes peaks.
    const std::vector<double> zero = {0.0, 1.0};
    FixedEnvelope mix3 = FixedEnvelope::LinearCombination(spectra, zero);
    CHECK(mix3.confs_no() == 3);
    CHECK(mix3.get_total_prob() == doctest::Approx(1.0));
}

TEST_CASE("bin() with zero width merges exactly-equal masses") {
    FixedEnvelope env = make_envelope({2.0, 1.0, 2.0, 3.0, 1.0}, {0.1, 0.2, 0.3, 0.15, 0.25});
    FixedEnvelope binned = env.bin(0.0);

    REQUIRE(binned.confs_no() == 3);
    std::vector<Peak> expect = {{1.0, 0.45}, {2.0, 0.4}, {3.0, 0.15}};
    CHECK(peaks_close(envelope_peaks(binned), expect, 1e-12));
    CHECK(binned.get_total_prob() == doctest::Approx(env.get_total_prob()));
}

TEST_CASE("bin() groups peaks into fixed-width bins") {
    // Bins of width 1 centred on integers: [x-0.5, x+0.5].
    FixedEnvelope env = make_envelope({0.9, 1.1, 1.4, 2.6, 7.0}, {0.1, 0.2, 0.3, 0.2, 0.2});
    FixedEnvelope binned = env.bin(1.0, 0.0);

    // Probability is conserved by binning.
    CHECK(binned.get_total_prob() == doctest::Approx(1.0));
    // Every bin middle is an integer multiple of the width, offset by `middle`.
    for (std::size_t i = 0; i < binned.confs_no(); ++i)
        CHECK(binned.mass(i) == doctest::Approx(std::round(binned.mass(i))));
    // 0.9, 1.1, 1.4 -> bin 1; 2.6 -> bin 3; 7.0 -> bin 7.
    std::vector<Peak> expect = {{1.0, 0.6}, {3.0, 0.2}, {7.0, 0.2}};
    CHECK(peaks_close(envelope_peaks(binned), expect, 1e-12));

    // A non-zero `middle` shifts the bin grid.
    FixedEnvelope shifted = env.bin(1.0, 0.5);
    CHECK(shifted.get_total_prob() == doctest::Approx(1.0));
    for (std::size_t i = 0; i < shifted.confs_no(); ++i)
        CHECK(shifted.mass(i) == doctest::Approx(std::round(shifted.mass(i) - 0.5) + 0.5));
}

TEST_CASE("bin() of a real spectrum conserves probability and reduces peaks") {
    for (double width : {0.01, 0.1, 1.0, 10.0}) {
        CAPTURE(width);
        FixedEnvelope env = FixedEnvelope::FromThreshold(Iso("C50H100N20O20S5"), 1e-10, true, false);
        const double total = env.get_total_prob();
        const double avg = env.empiric_average_mass();

        FixedEnvelope binned = env.bin(width);
        CHECK(binned.confs_no() <= env.confs_no());
        CHECK(binned.get_total_prob() == doctest::Approx(total).epsilon(1e-12));
        // Binning perturbs the mean by at most half a bin.
        CHECK(std::fabs(binned.empiric_average_mass() - avg) <= width * 0.5 + 1e-9);
        CHECK(is_sorted_by(binned.masses(), binned.confs_no()));
    }
}

TEST_CASE("Binned() named constructor") {
    Iso iso("C6H12O6");
    FixedEnvelope binned = FixedEnvelope::Binned(iso, 0.99, 1.0, 0.0);
    CHECK(binned.confs_no() > 0);
    CHECK(binned.get_total_prob() >= 0.99);
    for (std::size_t i = 0; i < binned.confs_no(); ++i)
        CHECK(binned.mass(i) == doctest::Approx(std::round(binned.mass(i))));

    // Equivalent to computing the envelope and binning it by hand.
    FixedEnvelope manual = FixedEnvelope::FromTotalProb(iso, 0.99, true, false).bin(1.0, 0.0);
    CHECK(binned.get_total_prob() == doctest::Approx(manual.get_total_prob()).epsilon(1e-6));

    // The move overload consumes the Iso.
    FixedEnvelope from_move = FixedEnvelope::Binned(Iso("H2O1"), 0.9, 0.5, 0.0);
    CHECK(from_move.confs_no() > 0);
    CHECK(from_move.get_total_prob() >= 0.9);
}

TEST_CASE("FromThreshold: absolute vs relative thresholds") {
    Iso iso("C50H100N20O20S5");
    const double mode_prob = std::exp(iso.getModeLProb());

    FixedEnvelope absolute = FixedEnvelope::FromThreshold(iso, 1e-4, true, false);
    FixedEnvelope relative = FixedEnvelope::FromThreshold(iso, 1e-4, false, false);

    // Relative cuts at 1e-4 * mode, which is a lower bar than absolute 1e-4
    // whenever the mode is below 1 — so it must keep at least as many peaks.
    CHECK(mode_prob < 1.0);
    CHECK(relative.confs_no() >= absolute.confs_no());

    for (std::size_t i = 0; i < absolute.confs_no(); ++i) CHECK(absolute.prob(i) >= 1e-4);
    for (std::size_t i = 0; i < relative.confs_no(); ++i)
        CHECK(relative.prob(i) >= 1e-4 * mode_prob * (1.0 - 1e-9));

    // A relative threshold of 1 keeps only the mode-probability peaks.
    FixedEnvelope only_mode = FixedEnvelope::FromThreshold(iso, 1.0, false, false);
    CHECK(only_mode.confs_no() >= 1);
    for (std::size_t i = 0; i < only_mode.confs_no(); ++i)
        CHECK(only_mode.prob(i) == doctest::Approx(mode_prob));

    // A threshold above 1 (absolute) keeps nothing.
    FixedEnvelope nothing = FixedEnvelope::FromThreshold(iso, 2.0, true, false);
    CHECK(nothing.confs_no() == 0);
    CHECK(nothing.get_total_prob() == doctest::Approx(0.0));
}

TEST_CASE("FromThreshold with negative or zero threshold enumerates everything") {
    for (double thr : {0.0, -1.0}) {
        CAPTURE(thr);
        FixedEnvelope env = FixedEnvelope::FromThreshold(Iso("C4H8O2"), thr, true, false);
        CHECK(env.get_total_prob() == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(env.confs_no() == enumerate_threshold_full("C4H8O2").size());
    }
}

TEST_CASE("FromTotalProb covers the requested probability") {
    Iso iso("C100H200N40O40");
    for (double target : {0.0, 0.1, 0.5, 0.9, 0.99, 0.999, 0.99999}) {
        CAPTURE(target);
        FixedEnvelope trimmed = FixedEnvelope::FromTotalProb(iso, target, true, false);
        FixedEnvelope untrimmed = FixedEnvelope::FromTotalProb(iso, target, false, false);

        CHECK(trimmed.get_total_prob() >= target - 1e-9);
        CHECK(untrimmed.get_total_prob() >= target - 1e-9);
        // Trimming drops the peaks that overshoot the target, so it can only be
        // smaller.
        CHECK(trimmed.confs_no() <= untrimmed.confs_no());
        CHECK(trimmed.get_total_prob() <= untrimmed.get_total_prob() + 1e-12);
    }

    // Larger targets need at least as many peaks.
    std::size_t previous = 0;
    for (double target : {0.5, 0.9, 0.99, 0.999}) {
        FixedEnvelope env = FixedEnvelope::FromTotalProb(iso, target, true, false);
        CHECK(env.confs_no() >= previous);
        previous = env.confs_no();
    }
}

TEST_CASE("FromTotalProb with target >= 1 gives the full distribution") {
    FixedEnvelope env = FixedEnvelope::FromTotalProb(Iso("C4H8O2"), 1.0, true, false);
    CHECK(env.get_total_prob() == doctest::Approx(1.0).epsilon(1e-9));

    FixedEnvelope over = FixedEnvelope::FromTotalProb(Iso("C4H8O2"), 2.0, true, false);
    CHECK(over.get_total_prob() == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("envelopes with configurations expose valid signatures") {
    Iso iso("C6H12O6");
    FixedEnvelope env = FixedEnvelope::FromTotalProb(iso, 0.99, true, true);
    REQUIRE(env.confs_no() > 0);
    CHECK(env.getAllDim() == iso.getAllDim());

    const std::vector<ElementSpec> es = elements_of("C6H12O6");
    for (std::size_t i = 0; i < env.confs_no(); ++i) {
        const int* c = env.conf(i);
        int off = 0;
        double mass = 0.0;
        for (const ElementSpec& e : es) {
            int total = 0;
            for (std::size_t j = 0; j < e.isotopes.size(); ++j) {
                REQUIRE(c[off] >= 0);
                total += c[off];
                mass += c[off] * e.isotopes[j].mass;
                ++off;
            }
            REQUIRE(total == e.count);
        }
        REQUIRE(env.mass(i) == doctest::Approx(mass).epsilon(1e-12));
    }
}

TEST_CASE("FromStochastic samples the requested number of molecules") {
    for (std::size_t n : {std::size_t(1), std::size_t(100), std::size_t(100000)}) {
        CAPTURE(n);
        FixedEnvelope env = FixedEnvelope::FromStochastic(Iso("C100H200N40O40"), n);
        CHECK(env.confs_no() > 0);
        // Every "probability" is an integer count, and they sum to n.
        double total = 0.0;
        for (std::size_t i = 0; i < env.confs_no(); ++i) {
            CHECK(env.prob(i) > 0.0);
            CHECK(env.prob(i) == doctest::Approx(std::round(env.prob(i))));
            total += env.prob(i);
        }
        CHECK(total == doctest::Approx(static_cast<double>(n)));
    }
}

TEST_CASE("FromStochastic peaks lie in the molecule's mass range") {
    Iso iso("C100H200N40O40");
    const double lo = iso.getLightestPeakMass();
    const double hi = iso.getHeaviestPeakMass();
    FixedEnvelope env = FixedEnvelope::FromStochastic(iso, 10000, 0.9999, 5.0, true);
    REQUIRE(env.confs_no() > 0);
    CHECK(env.getAllDim() == iso.getAllDim());
    for (std::size_t i = 0; i < env.confs_no(); ++i) {
        CHECK(env.mass(i) >= lo - 1e-9);
        CHECK(env.mass(i) <= hi + 1e-9);
    }
    // The sample mean must be close to the theoretical mean (10k draws).
    CHECK(env.empiric_average_mass() ==
          doctest::Approx(iso.getTheoreticalAverageMass()).epsilon(1e-4));
}

TEST_CASE("resample redistributes a fixed ionic current") {
    FixedEnvelope env = FixedEnvelope::FromTotalProb(Iso("C100H200N40O40"), 0.999, true, false);
    env.normalize();
    const std::size_t n = env.confs_no();
    REQUIRE(n > 10);

    env.resample(100000);

    // resample() replaces probabilities with counts, so the cached total_prob
    // from the normalize() above must have been invalidated.
    CHECK(env.get_total_prob() == doctest::Approx(100000.0));

    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(env.prob(i) >= 0.0);
        CHECK(env.prob(i) == doctest::Approx(std::round(env.prob(i))));
        total += env.prob(i);
    }
    CHECK(total == doctest::Approx(100000.0));

    // Resampling to a single molecule puts all the current in one peak.
    FixedEnvelope one = FixedEnvelope::FromTotalProb(Iso("C10H20"), 0.999, true, false);
    one.normalize();
    one.resample(1);
    double sum = 0.0;
    for (std::size_t i = 0; i < one.confs_no(); ++i) sum += one.prob(i);
    CHECK(sum == doctest::Approx(1.0));
}

TEST_CASE("FromThreshold const-Iso overloads do not consume the Iso") {
    Iso iso("C6H12O6");
    FixedEnvelope a = FixedEnvelope::FromThreshold(iso, 1e-6, true, false);
    FixedEnvelope b = FixedEnvelope::FromThreshold(iso, 1e-6, true, false);
    FixedEnvelope c = FixedEnvelope::FromTotalProb(iso, 0.99, true, false);
    FixedEnvelope d = FixedEnvelope::FromStochastic(iso, 1000);
    CHECK(a.confs_no() == b.confs_no());
    CHECK(c.confs_no() > 0);
    CHECK(d.confs_no() > 0);
    // The Iso is still usable afterwards.
    CHECK(iso.getMonoisotopicPeakMass() == doctest::Approx(Iso("C6H12O6").getMonoisotopicPeakMass()));
}
