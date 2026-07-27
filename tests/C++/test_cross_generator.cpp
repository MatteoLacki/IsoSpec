// Cross-generator invariants.
//
// The threshold, ordered, and layered generators traverse the configuration
// space in different orders and with different stopping rules, but over the FULL
// support of a molecule made of real elements they must all enumerate the same
// multiset of (mass, prob) peaks.  This is the strongest correctness check the
// suite has that does not depend on a single generator being its own oracle.

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

TEST_CASE("full support agrees across threshold / ordered / layered") {
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);

        std::vector<Peak> thr = enumerate_threshold_full(f);
        std::vector<Peak> ord = enumerate_ordered_full(f);
        std::vector<Peak> lay = enumerate_layered_full(f);

        REQUIRE(thr.size() > 0);
        CHECK(thr.size() == ord.size());
        CHECK(thr.size() == lay.size());

        CHECK(peaks_close(thr, ord));
        CHECK(peaks_close(thr, lay));

        // Full support sums to 1.
        CHECK(std::fabs(total_prob(thr) - 1.0) < 1e-9);
        CHECK(std::fabs(total_prob(ord) - 1.0) < 1e-9);
        CHECK(std::fabs(total_prob(lay) - 1.0) < 1e-9);
    }
}

TEST_CASE("ordered generator yields strictly non-increasing probability") {
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);
        IsoOrderedGenerator g(f);
        double last = 2.0;  // > any probability
        std::size_t n = 0;
        while (g.advanceToNextConfiguration()) {
            double p = g.prob();
            CHECK(p <= last);
            last = p;
            ++n;
        }
        CHECK(n > 0);
    }
}

TEST_CASE("FixedEnvelope::FromThreshold materialization matches the generator") {
    // The materialized envelope must contain exactly the generator's output.
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);
        std::vector<Peak> gen = enumerate_threshold_full(f);

        FixedEnvelope env = FixedEnvelope::FromThreshold(Iso(f), 0.0, true, false);
        std::vector<Peak> mat;
        mat.reserve(env.confs_no());
        for (std::size_t i = 0; i < env.confs_no(); ++i)
            mat.push_back({env.masses()[i], env.probs()[i]});

        CHECK(gen.size() == mat.size());
        CHECK(peaks_close(gen, mat));
    }
}

TEST_CASE("FromTotalProb captures at least the requested probability") {
    for (const char* f : big_formulas()) {
        INFO("formula=" << f);
        for (double target : {0.5, 0.9, 0.99, 0.999}) {
            CAPTURE(target);
            FixedEnvelope env = FixedEnvelope::FromTotalProb(Iso(f), target, true, false);
            CHECK(env.get_total_prob() >= target - 1e-9);
        }
    }
}
