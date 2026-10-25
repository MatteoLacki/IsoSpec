// Bit-identity between the scalar and SIMD fill paths.
//
// FixedEnvelope::FromThreshold(..., tgetConfs=true) runs the scalar per-config
// fill; tgetConfs=false runs the SIMD batch fill (simd_massprobs).  Both emit
// configurations in identical order, so they must agree index-by-index and
// bit-for-bit — the SIMD path performs the same FP ops (partialProbs[1]*p,
// partialMasses[1]+m) in the same order.
//
// This is the regression net that guards the SIMD generator work and MUST stay
// green before any change to the shared marginal/guardian machinery.  Binned()'s
// full-enumeration filler is batched over the same generator and is covered at
// the bottom of this file.  The layered generator's fills are batched too
// (FromTotalProb without configurations, and Binned()'s probability-bounded
// filler, both via IsoLayeredGenerator::batch_fill_run()), and are held to the
// same bit-identity below.

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using test_helpers::big_formulas;
using test_helpers::small_formulas;

namespace {

// Compare the two fill paths for one (formula, threshold).  Returns the config
// count so callers can also sanity-check it is nonzero.
std::size_t check_threshold_paths(const char* formula, double threshold) {
    FixedEnvelope scalar = FixedEnvelope::FromThreshold(Iso(formula), threshold, true, true);
    FixedEnvelope simd   = FixedEnvelope::FromThreshold(Iso(formula), threshold, true, false);

    REQUIRE(scalar.confs_no() == simd.confs_no());

    const std::size_t n = scalar.confs_no();
    const double* sm = scalar.masses();
    const double* sp = scalar.probs();
    const double* im = simd.masses();
    const double* ip = simd.probs();

    for (std::size_t i = 0; i < n; ++i) {
        // Exact equality: same order, same arithmetic.
        REQUIRE(sm[i] == im[i]);
        REQUIRE(sp[i] == ip[i]);
    }
    // total_prob is a downstream rescan; it too must match exactly.
    REQUIRE(scalar.get_total_prob() == simd.get_total_prob());
    return n;
}

}  // namespace

TEST_CASE("scalar vs SIMD FromThreshold: small molecules") {
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);
        for (double thr : {1e-2, 1e-4, 1e-6, 1e-8, 1e-12}) {
            CAPTURE(thr);
            check_threshold_paths(f, thr);
        }
    }
}

TEST_CASE("scalar vs SIMD FromThreshold: large multi-run molecules") {
    // These have many marginal-0 runs (many carries), exercising the post-carry
    // convention bridge and the unaligned load/store tail in the SIMD path.
    for (const char* f : big_formulas()) {
        INFO("formula=" << f);
        for (double thr : {1e-4, 1e-8, 1e-10}) {
            CAPTURE(thr);
            std::size_t n = check_threshold_paths(f, thr);
            CHECK(n > 0);
        }
    }
}

TEST_CASE("scalar vs SIMD: single leftover config (< W tail only)") {
    // A threshold high enough to keep just the monoisotopic peak: exercises the
    // path where the SIMD batch never fires and only the scalar tail runs.
    for (const char* f : {"C1", "H2O1", "C10"}) {
        INFO("formula=" << f);
        check_threshold_paths(f, 0.999999);
    }
}

// ---------------------------------------------------------------------------
// FromTotalProb: with configurations the fill is scalar, without them it is
// batched (FixedEnvelope::store_run()), including stopping partway through a
// marginal-0 run once the target is met -- the batch overshoots and is trimmed
// back.  Both must walk the layers identically, stop at the same configuration
// and agree bit-for-bit.
// ---------------------------------------------------------------------------

namespace {

std::size_t check_total_prob_paths(const char* formula, double target) {
    FixedEnvelope scalar = FixedEnvelope::FromTotalProb(Iso(formula), target, false, true);
    FixedEnvelope simd   = FixedEnvelope::FromTotalProb(Iso(formula), target, false, false);

    REQUIRE(scalar.confs_no() == simd.confs_no());

    const std::size_t n = scalar.confs_no();
    for (std::size_t i = 0; i < n; ++i) {
        REQUIRE(scalar.masses()[i] == simd.masses()[i]);
        REQUIRE(scalar.probs()[i] == simd.probs()[i]);
    }
    REQUIRE(scalar.get_total_prob() == simd.get_total_prob());
    return n;
}

// With optimize=true the fill is followed by quicktrim, whose pivots come from
// a global Mersenne twister, so the two runs end up with the kept
// configurations in different orders.  The set of them, however, must match.
void check_total_prob_paths_optimized(const char* formula, double target) {
    FixedEnvelope scalar = FixedEnvelope::FromTotalProb(Iso(formula), target, true, true);
    FixedEnvelope simd   = FixedEnvelope::FromTotalProb(Iso(formula), target, true, false);

    REQUIRE(scalar.confs_no() == simd.confs_no());

    std::vector<test_helpers::Peak> s, v;
    for (std::size_t i = 0; i < scalar.confs_no(); ++i)
        s.push_back({scalar.masses()[i], scalar.probs()[i]});
    for (std::size_t i = 0; i < simd.confs_no(); ++i)
        v.push_back({simd.masses()[i], simd.probs()[i]});
    test_helpers::sort_peaks(s);
    test_helpers::sort_peaks(v);

    for (std::size_t i = 0; i < s.size(); ++i) {
        REQUIRE(s[i].mass == v[i].mass);
        REQUIRE(s[i].prob == v[i].prob);
    }
}

}  // namespace

TEST_CASE("FromTotalProb, both fills agree: small molecules") {
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);
        for (double p : {0.5, 0.9, 0.99, 0.9999}) {
            CAPTURE(p);
            check_total_prob_paths(f, p);
            check_total_prob_paths_optimized(f, p);
        }
    }
}

TEST_CASE("FromTotalProb, both fills agree: large multi-layer molecules") {
    for (const char* f : big_formulas()) {
        INFO("formula=" << f);
        for (double p : {0.9, 0.999, 0.99999}) {
            CAPTURE(p);
            std::size_t n = check_total_prob_paths(f, p);
            CHECK(n > 0);
            check_total_prob_paths_optimized(f, p);
        }
    }
}

TEST_CASE("FromTotalProb, both fills agree: targets met almost immediately") {
    // Tiny targets are met within the first few configurations, so the fill
    // stops in the middle of the very first layer.
    for (const char* f : {"C1", "H2O1", "C10", "C6H12O6", "C1000"}) {
        INFO("formula=" << f);
        for (double p : {1e-9, 0.01, 0.1, 0.3}) {
            CAPTURE(p);
            check_total_prob_paths(f, p);
        }
    }
}

TEST_CASE("FromTotalProb, both fills agree: degenerate targets") {
    // >= 1.0 delegates to threshold_init; <= 0.0 produces nothing at all.
    for (const char* f : {"H2O1", "C6H12O6"}) {
        INFO("formula=" << f);
        check_total_prob_paths(f, 1.0);
        check_total_prob_paths(f, 0.0);
    }
}

// ---------------------------------------------------------------------------
// Binned: its target>=1.0 filler runs over the threshold generator and *is*
// batched (measured: 89% of configurations, worth up to 2x on a large full
// enumeration into wide bins); its target<1.0 one runs over the layered
// generator and is not.  There is no scalar twin to compare either against (no
// tgetConfs variant), so the reference is built here, from the very generators
// Binned() drives, with the scatter done one configuration at a time in the same
// order.  That pins down the batched fill, the vectorised bin-index arithmetic,
// and the layer stepping alike.
// ---------------------------------------------------------------------------

namespace {

std::ptrdiff_t ref_bin_index(double mass, double hwmm, double inv_bin_width) {
    return static_cast<std::ptrdiff_t>(std::floor((mass + hwmm) * inv_bin_width));
}

// Unlike threshold_init()'s fill, the binned ones do not merely store what the
// generator yields, they accumulate it: `acc[bin] += partialProb * marginalProb`.
// That is a multiply-add, which a compiler with -ffp-contract=fast (the default,
// and what this project builds with) fuses in the scalar rendering and cannot
// fuse in the batched one, where the multiply has already happened in a vector
// register.  So the two agree to within a rounding of the last bit per
// accumulation rather than exactly -- verified by rebuilding this whole file,
// library included, with -ffp-contract=off, under which every comparison below
// holds bit-for-bit.  Bin count and bin masses are pure integer-grid arithmetic
// with no such freedom and are still required to match exactly; any mistake in
// the batching or the stopping rule shows up there, or as a probability
// difference orders of magnitude beyond this tolerance.
constexpr double kBinProbRelTolerance = 1e-12;

// Bins accumulated the slow way: an ordered list of (bin index, probability),
// built by a scalar scatter, then compacted the way Binned() compacts its
// accumulator.
std::vector<test_helpers::Peak> reference_binned(const char* formula, double target_total_prob,
                                                 double bin_width, double bin_middle) {
    std::vector<test_helpers::Peak> out;
    if (target_total_prob <= 0.0)
        return out;

    // Binned()'s grid: the middle reduced mod the width, and the origin moved
    // to two bins below the lightest mass, so that the index arithmetic yields
    // accumulator slots directly.
    Iso iso(formula);
    bin_middle = std::fmod(bin_middle, bin_width);
    const double rounding = iso.getMassRoundingBound();
    const double lo = iso.getLightestPeakMass() - rounding;
    const double hi = iso.getHeaviestPeakMass() + rounding;
    const double inv_bin_width = 1.0 / bin_width;
    const std::ptrdiff_t idx_min = ref_bin_index(lo, 0.5 * bin_width - bin_middle, inv_bin_width) - 2;
    const double hwmm = (0.5 * bin_width - bin_middle) - static_cast<double>(idx_min) * bin_width;
    const std::size_t no_bins = static_cast<std::size_t>(ref_bin_index(hi, hwmm, inv_bin_width)) + 2;
    const std::ptrdiff_t idx_max = static_cast<std::ptrdiff_t>(no_bins);

    std::vector<double> acc(no_bins, 0.0);
    auto at = [&](std::ptrdiff_t i) -> double& { return acc.at(static_cast<std::size_t>(i)); };

    std::ptrdiff_t nonzero_idx = 0;
    bool seeded = false;

    if (target_total_prob >= 1.0) {
        IsoThresholdGenerator g(Iso(formula), 0.0, true);
        while (g.advanceToNextConfiguration()) {
            const double prob = g.prob();
            const std::ptrdiff_t bin = ref_bin_index(g.mass(), hwmm, inv_bin_width);
            if (!seeded) {
                if (prob == 0.0) continue;  // the seed must land on a populated bin
                nonzero_idx = bin;
                seeded = true;
                at(bin) = prob;
                continue;
            }
            at(bin) += prob;
        }
    } else {
        IsoLayeredGenerator g(Iso(formula), 1000, 1000, true,
                              std::min<double>(target_total_prob, 0.9999));
        double prob_so_far = 0.0;
        const double sum_above = std::log1p(-target_total_prob) - 2.3025850929940455;
        bool done = false;
        do {
            while (g.advanceToNextConfigurationWithinLayer()) {
                const double prob = g.prob();
                if (prob == 0.0) continue;
                const std::ptrdiff_t bin = ref_bin_index(g.mass(), hwmm, inv_bin_width);
                at(bin) += prob;
                if (!seeded) { nonzero_idx = bin; seeded = true; }
                prob_so_far += prob;
                if (prob_so_far >= target_total_prob) { done = true; break; }
            }
            if (done) break;
            double layer_delta = sum_above - std::log1p(-prob_so_far);
            layer_delta = std::max(std::min(layer_delta, -0.1),
                                   static_cast<double>(ISOSPEC_BINNED_LAYER_MAXSTEP));
            if (!g.nextLayer(layer_delta)) break;
        } while (true);
    }

    if (!seeded)
        return out;

    const std::size_t distance_10da = static_cast<std::size_t>(10.0 / bin_width) + 1;
    std::size_t empty_steps = 0;
    for (std::ptrdiff_t ii = nonzero_idx; empty_steps < distance_10da; ) {
        if (at(ii) > 0.0) {
            empty_steps = 0;
            out.push_back({static_cast<double>(ii + idx_min) * bin_width + bin_middle, at(ii)});
        } else {
            empty_steps++;
        }
        if (ii == 0) break;
        ii--;
    }
    empty_steps = 0;
    for (std::ptrdiff_t ii = nonzero_idx + 1; ii < idx_max && empty_steps < distance_10da; ii++) {
        if (at(ii) > 0.0) {
            empty_steps = 0;
            out.push_back({static_cast<double>(ii + idx_min) * bin_width + bin_middle, at(ii)});
        } else {
            empty_steps++;
        }
    }
    return out;
}

void check_binned(const char* formula, double target, double bin_width, double bin_middle = 0.0) {
    FixedEnvelope binned = FixedEnvelope::Binned(Iso(formula), target, bin_width, bin_middle);
    std::vector<test_helpers::Peak> ref = reference_binned(formula, target, bin_width, bin_middle);

    REQUIRE(binned.confs_no() == ref.size());
    for (std::size_t i = 0; i < ref.size(); ++i) {
        REQUIRE(binned.masses()[i] == ref[i].mass);
        REQUIRE(std::abs(binned.probs()[i] - ref[i].prob)
                <= kBinProbRelTolerance * std::abs(ref[i].prob));
    }
}

}  // namespace

TEST_CASE("Binned vs reference scatter: full enumeration (threshold-driven, batched)") {
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);
        for (double w : {1.0, 0.1, 0.01}) {
            CAPTURE(w);
            check_binned(f, 1.0, w);
            check_binned(f, 1.0, w, 0.5 * w);
        }
    }
}

TEST_CASE("Binned vs reference scatter: probability-bounded (layer-driven)") {
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);
        for (double p : {0.5, 0.99, 0.9999}) {
            CAPTURE(p);
            check_binned(f, p, 1.0);
            check_binned(f, p, 0.05);
        }
    }
    for (const char* f : big_formulas()) {
        INFO("formula=" << f);
        for (double p : {0.9, 0.999}) {
            CAPTURE(p);
            check_binned(f, p, 1.0);
        }
    }
}
