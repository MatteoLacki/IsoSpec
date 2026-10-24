// FixedEnvelope: the materialized-spectrum type.
//
// Construction (all four named constructors), value semantics, sorting,
// statistics, arithmetic, binning and resampling.  Small hand-built envelopes
// are used wherever the expected answer can be written down exactly.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

// --- ABI guard ---------------------------------------------------------------
//
// FixedEnvelope is returned by value from FromThreshold()/FromTotalProb()/…, so
// callers stack-allocate it and its layout is part of the shared library's ABI.
// Growing it (2.4.0 took it from 96 to 216 bytes) breaks every program still
// linked against the old library, and nothing else in the build notices.
//
// The soname carries MAJOR.MINOR (see CMakeLists.txt), so a minor bump already
// covers a break like that one.  What this assertion is for is the case that
// scheme does NOT cover: a layout change slipped into a *patch* release, where
// the soname stays put.  If it fires, either revert the layout change or make
// the release a minor one.
//
// Pinned only on LP64 Unix, the reference platform: the size legitimately
// differs elsewhere (aligned_unique_ptr's VM bookkeeping carries one more word
// on Windows, and pointers are narrower on 32-bit).
#if defined(__unix__) && defined(__LP64__)
static_assert(sizeof(FixedEnvelope) == 216,
              "FixedEnvelope's layout changed: this is an ABI break. See the note above.");
#endif

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

TEST_CASE("resample with few molecules jumps over configurations") {
    // While samples * (the current configuration's share of what is left) is
    // below beta_bias, resample() places molecules one at a time by beta
    // jumps, and a jump can carry it past configurations that receive
    // nothing. A large resample is almost all binomial steps; a handful of
    // molecules over hundreds of small peaks takes the jumping path on nearly
    // every draw, and walks onto the end-of-array sentinel when the last
    // molecule lands in the tail.
    FixedEnvelope source = FixedEnvelope::FromTotalProb(Iso("C100H200N40O40S2"), 0.9999, true, false);
    source.normalize();
    const std::size_t n = source.confs_no();
    REQUIRE(n > 100);
    const double true_mean = source.empiric_average_mass();

    double mass_sum = 0.0;
    double molecules = 0.0;
    for (std::size_t samples : {1, 2, 5, 20}) {
        for (int rep = 0; rep < 100; ++rep) {
            CAPTURE(samples);
            CAPTURE(rep);
            FixedEnvelope env = source;
            env.resample(samples);
            REQUIRE(env.confs_no() == n);
            double total = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                // Includes the sentinel slot, which must come back as a count.
                REQUIRE(env.prob(i) >= 0.0);
                REQUIRE(env.prob(i) <= static_cast<double>(samples));
                REQUIRE(env.prob(i) == std::round(env.prob(i)));
                total += env.prob(i);
                mass_sum += env.prob(i) * env.mass(i);
            }
            REQUIRE(total == static_cast<double>(samples));
            molecules += total;
        }
    }
    // 2800 molecules with a mass spread of 1.39 Da: the standard error of the
    // mean is 0.026 Da, and the tolerance below is ~6 of them.
    CHECK(mass_sum / molecules == doctest::Approx(true_mean).epsilon(0.15 / true_mean));
}

TEST_CASE("FromTotalProb's trim keeps exactly the smallest covering set") {
    // quicktrim narrows the last layer down by partitioning; a segment
    // narrowed to a single configuration is kept or dropped on the prefix sum
    // alone. Sweeping the target through many values puts that configuration
    // on both sides of the boundary.
    for (const char* formula : {"C10H20O5", "C20H30N5O5S1"}) {
        CAPTURE(formula);
        std::vector<double> probs;
        for (const Peak& p : enumerate_threshold_full(formula)) probs.push_back(p.prob);
        std::sort(probs.begin(), probs.end(), std::greater<double>());
        std::vector<double> prefix(probs.size());
        double sum = 0.0;
        for (std::size_t i = 0; i < probs.size(); ++i) prefix[i] = (sum += probs[i]);

        for (int step = 1; step < 400; ++step) {
            const double target = step / 400.0;
            CAPTURE(target);
            // Summation order differs from the library's; skip targets that
            // a rounding difference could put on the other side of a prefix.
            const auto it = std::lower_bound(prefix.begin(), prefix.end(), target);
            REQUIRE(it != prefix.end());
            if (std::fabs(*it - target) < 1e-12 || (it != prefix.begin() && std::fabs(*(it - 1) - target) < 1e-12))
                continue;
            const std::size_t expected = static_cast<std::size_t>(it - prefix.begin()) + 1;
            CHECK(FixedEnvelope::FromTotalProb(Iso(formula), target, true, false).confs_no() == expected);
        }
    }
}

TEST_CASE("Binned: a non-positive target gives an empty envelope") {
    for (double target : {0.0, -1.0}) {
        CAPTURE(target);
        CHECK(FixedEnvelope::Binned(Iso("C2"), target, 1.0, 0.0).confs_no() == 0);
    }
}

TEST_CASE("Binned skips configurations whose probability underflows to zero") {
    // {1, 1e-200} sums to exactly 1 in double. A target above 1 is never
    // reached, so the layered fill runs through every configuration,
    // including the two whose product underflows to 0.0. Those must not
    // leave an (empty) bin behind.
    int isotope_numbers[1] = {2};
    int atom_counts[1] = {3};
    double masses[2] = {1.0, 2.0};
    double probs[2] = {1.0, 1e-200};
    FixedEnvelope binned = FixedEnvelope::Binned(Iso(1, isotope_numbers, atom_counts, masses, probs), 2.0, 1.0, 0.0);
    REQUIRE(binned.confs_no() == 2);
    binned.sort_by_mass();
    CHECK(binned.mass(0) == doctest::Approx(3.0));
    CHECK(binned.prob(0) == doctest::Approx(1.0));
    CHECK(binned.mass(1) == doctest::Approx(4.0));
    CHECK(binned.prob(1) == doctest::Approx(3e-200));
}

TEST_CASE("FromTotalProb with an unreachable target keeps every configuration") {
    // A target above the molecule's total of 1: the layers run out before it
    // is met, and whatever was stored is the answer, trimmed or not.
    int isotope_numbers[1] = {2};
    int atom_counts[1] = {3};
    double masses[2] = {1.0, 2.0};
    double probs[2] = {1.0, 1e-200};
    for (bool optimize : {true, false})
        for (bool get_confs : {true, false}) {
            CAPTURE(optimize);
            CAPTURE(get_confs);
            FixedEnvelope env = FixedEnvelope::FromTotalProb(Iso(1, isotope_numbers, atom_counts, masses, probs),
                                                             2.0, optimize, get_confs);
            CHECK(env.confs_no() == 4);
            CHECK(env.get_total_prob() == doctest::Approx(1.0));
        }
}

TEST_CASE("Binned rejects a grid whose last bin index would overflow") {
    // The first bin index fits in ptrdiff_t, but the ~10^5 bins above it do
    // not: the second range check, not the first, has to catch this.
    Iso iso("C100");
    const double width = 1e-3;
    const double middle = iso.getLightestPeakMass() + 0.5 * width - (std::ldexp(1.0, 63) - 50000.0) * width;
    CHECK_THROWS_AS(FixedEnvelope::Binned(Iso("C100"), 0.99, width, middle), std::length_error);
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

TEST_CASE("LinearCombination rejects unequal vector lengths") {
    auto env = make_envelope({1.0}, {1.0});
    const std::vector<const FixedEnvelope*> spectra{&env, &env};
    for(const std::vector<double>& intensities :
            {std::vector<double>{}, {1.0}, {1.0, 2.0, 3.0}})
        CHECK_THROWS_AS(FixedEnvelope::LinearCombination(spectra, intensities), std::invalid_argument);
}

TEST_CASE("binning validates grid parameters before processing peaks") {
    auto env = make_envelope({1.0, 1.0, 2.0}, {0.2, 0.3, 0.5});
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for(double width : {-1.0, nan, inf, -inf}) {
        CHECK_THROWS_AS(env.bin(width), std::invalid_argument);
        CHECK_THROWS_AS(FixedEnvelope::Binned(Iso("C2"), 0.99, width), std::invalid_argument);
        CHECK_THROWS_AS(FixedEnvelope::Binned(Iso("C2"), 0.0, width), std::invalid_argument);
    }
    for(double middle : {nan, inf, -inf}) {
        CHECK_THROWS_AS(env.bin(1.0, middle), std::invalid_argument);
        CHECK_THROWS_AS(FixedEnvelope::Binned(Iso("C2"), 0.99, 1.0, middle), std::invalid_argument);
    }
    CHECK_THROWS_AS(FixedEnvelope::Binned(Iso("C2"), 0.99, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(FixedEnvelope::Binned(Iso("C2"), 0.99, 1e-300), std::length_error);
    CHECK_THROWS_AS(FixedEnvelope::Binned(Iso("C2"), 0.99, 1.0, 1e300), std::length_error);
    auto exact = env.bin(0.0);
    REQUIRE(exact.confs_no() == 2);
    CHECK(exact.mass(0) == 1.0);
    CHECK(exact.prob(0) == doctest::Approx(0.5));
    CHECK(exact.mass(1) == 2.0);
    CHECK(exact.prob(1) == doctest::Approx(0.5));
}

TEST_CASE("Scalar parameters that would loop forever or reach undefined behaviour are rejected") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const char* f = "C10H20O5S1";
    std::mt19937 rng = seeded_rng(1);

    SUBCASE("threshold: NaN rejected, <= 0 still means everything") {
        CHECK_THROWS_AS(IsoThresholdGenerator(Iso(f), nan), std::invalid_argument);
        CHECK_THROWS_AS(FixedEnvelope::FromThreshold(Iso(f), nan, false), std::invalid_argument);
        CHECK(FixedEnvelope::FromThreshold(Iso(f), -1.0, true).confs_no() ==
              FixedEnvelope::FromThreshold(Iso(f), 0.0, true).confs_no());
    }
    SUBCASE("target total probability: NaN rejected, out-of-range values still clamp") {
        CHECK_THROWS_AS(FixedEnvelope::FromTotalProb(Iso(f), nan, true), std::invalid_argument);
        CHECK_THROWS_AS(FixedEnvelope::Binned(Iso(f), nan, 1.0), std::invalid_argument);
        CHECK(FixedEnvelope::FromTotalProb(Iso(f), -1.0, true).confs_no() == 0);
        CHECK(FixedEnvelope::FromTotalProb(Iso(f), 2.0, true).get_total_prob() == doctest::Approx(1.0));
    }
    SUBCASE("stochastic precision in (0, 1], beta_bias nonnegative") {
        for (double p : {nan, -inf, -1.0, 0.0, 1.0 + 1e-15, 2.0, inf})
            CHECK_THROWS_AS(FixedEnvelope::FromStochastic(Iso(f), 1000, p, 5.0, false, rng), std::invalid_argument);
        for (double b : {nan, -1.0, -inf})
            CHECK_THROWS_AS(IsoStochasticGenerator(Iso(f), 1000, 0.9999, b, rng), std::invalid_argument);
        CHECK(FixedEnvelope::FromStochastic(Iso(f), 1000, 1.0, 5.0, false, rng).get_total_prob() == 1000.0);
        CHECK(FixedEnvelope::FromStochastic(Iso(f), 1000, 0.9999, inf, false, rng).get_total_prob() == 1000.0);
    }
    SUBCASE("envelope operations") {
        FixedEnvelope e = FixedEnvelope::FromTotalProb(Iso(f), 0.999, true);
        FixedEnvelope other = FixedEnvelope::FromTotalProb(Iso("C10H22O5S1"), 0.999, true);
        for (double x : {nan, inf, -inf, -1.0})
        {
            CHECK_THROWS_AS(e.scale(x), std::invalid_argument);
            const FixedEnvelope* s[1] = {&e};
            CHECK_THROWS_AS(FixedEnvelope::LinearCombination(s, &x, 1), std::invalid_argument);
            CHECK_THROWS_AS(e.AbyssalWassersteinDistance(other, 1.0, x), std::invalid_argument);
            CHECK_THROWS_AS(e.WassersteinMatch(other, 0.1, x), std::invalid_argument);
        }
        for (double x : {nan, inf, -inf})
            CHECK_THROWS_AS(e.shift_mass(x), std::invalid_argument);
        for (double x : {nan, -inf, -1.0})
        {
            CHECK_THROWS_AS(e.AbyssalWassersteinDistance(other, x), std::invalid_argument);
            CHECK_THROWS_AS(e.WassersteinMatch(other, x), std::invalid_argument);
            CHECK_THROWS_AS(e.resample(1000, x, rng), std::invalid_argument);
        }
        // Infinite depth and flow distance stay meaningful, as do a zero factor and beta_bias.
        CHECK(e.AbyssalWassersteinDistance(other, inf) > 0.0);
        CHECK(std::get<1>(e.WassersteinMatch(other, inf)) >= 0.0);
        e.normalize();
        e.resample(1000, 0.0, rng);
        CHECK(e.get_total_prob() == 1000.0);
        e.scale(0.0);
        CHECK(e.get_total_prob() == 0.0);
    }
}
