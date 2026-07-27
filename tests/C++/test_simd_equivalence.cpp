// Bit-identity between the scalar and SIMD fill paths.
//
// FixedEnvelope::FromThreshold(..., tgetConfs=true) runs the scalar per-config
// fill; tgetConfs=false runs the SIMD batch fill (simd_massprobs).  Both emit
// configurations in identical order, so they must agree index-by-index and
// bit-for-bit — the SIMD path performs the same FP ops (partialProbs[1]*p,
// partialMasses[1]+m) in the same order.
//
// This is the regression net that guards the SIMD generator work and MUST stay
// green before any change to the shared marginal/guardian machinery (e.g. the
// planned LayeredMarginal SIMD extension).

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
