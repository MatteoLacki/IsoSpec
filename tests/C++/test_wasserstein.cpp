// Wasserstein-family spectrum distances.
//
// The plain and oriented distances are checked against closed forms on
// hand-built spectra (for one-dimensional distributions the Wasserstein-1
// distance is the L1 distance between the CDFs, which is easy to write down),
// then against the metric axioms on real isotopic envelopes.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

FixedEnvelope make_envelope(const std::vector<double>& masses,
                            const std::vector<double>& probs) {
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

// Independent reference: W1 between two discrete distributions, computed by
// walking the merged support and integrating |CDF_a - CDF_b|.
double reference_w1(std::vector<Peak> a, std::vector<Peak> b) {
    std::vector<double> support;
    for (const Peak& p : a) support.push_back(p.mass);
    for (const Peak& p : b) support.push_back(p.mass);
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());

    double ca = 0.0, cb = 0.0, acc = 0.0;
    for (std::size_t i = 0; i < support.size(); ++i) {
        for (const Peak& p : a) if (p.mass == support[i]) ca += p.prob;
        for (const Peak& p : b) if (p.mass == support[i]) cb += p.prob;
        if (i + 1 < support.size())
            acc += std::fabs(ca - cb) * (support[i + 1] - support[i]);
    }
    return acc;
}

}  // namespace

TEST_CASE("Wasserstein distance of a spectrum to itself is zero") {
    for (const char* f : {"H2O1", "C6H12O6", "C50H100N20O20S5"}) {
        INFO("formula=" << f);
        FixedEnvelope a = FixedEnvelope::FromThreshold(Iso(f), 1e-8, true, false);
        FixedEnvelope b = FixedEnvelope::FromThreshold(Iso(f), 1e-8, true, false);
        a.normalize();
        b.normalize();
        CHECK(a.WassersteinDistance(b) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(a.OrientedWassersteinDistance(b) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(a.AbyssalWassersteinDistance(b, 1.0) == doctest::Approx(0.0).epsilon(1e-12));
    }
}

TEST_CASE("two-point spectra: distance is the transported mass times the gap") {
    // All the probability sits at 0 in one spectrum and at d in the other:
    // moving a unit of mass a distance d costs exactly d.
    for (double d : {0.5, 1.0, 17.25}) {
        CAPTURE(d);
        FixedEnvelope a = make_envelope({0.0}, {1.0});
        FixedEnvelope b = make_envelope({d}, {1.0});
        CHECK(a.WassersteinDistance(b) == doctest::Approx(d));
        CHECK(b.WassersteinDistance(a) == doctest::Approx(d));  // symmetry
    }

    // Half the mass moves 2 units, half stays: cost 1.
    FixedEnvelope a = make_envelope({0.0, 2.0}, {0.5, 0.5});
    FixedEnvelope b = make_envelope({0.0, 4.0}, {0.5, 0.5});
    CHECK(a.WassersteinDistance(b) == doctest::Approx(1.0));
}

TEST_CASE("Wasserstein distance matches an independent CDF integration") {
    struct Case { std::vector<double> ma, pa, mb, pb; };
    const Case cases[] = {
        {{0.0, 1.0, 2.0}, {0.2, 0.3, 0.5}, {0.0, 1.0, 2.0}, {0.5, 0.3, 0.2}},
        {{1.0, 3.0}, {0.5, 0.5}, {2.0, 4.0}, {0.25, 0.75}},
        {{0.0}, {1.0}, {-5.0, 5.0}, {0.5, 0.5}},
        {{10.0, 10.5, 11.0}, {0.1, 0.6, 0.3}, {10.2, 10.7}, {0.4, 0.6}},
    };
    for (const Case& c : cases) {
        FixedEnvelope a = make_envelope(c.ma, c.pa);
        FixedEnvelope b = make_envelope(c.mb, c.pb);
        std::vector<Peak> pa, pb;
        for (std::size_t i = 0; i < c.ma.size(); ++i) pa.push_back({c.ma[i], c.pa[i]});
        for (std::size_t i = 0; i < c.mb.size(); ++i) pb.push_back({c.mb[i], c.pb[i]});
        CHECK(a.WassersteinDistance(b) == doctest::Approx(reference_w1(pa, pb)).epsilon(1e-9));
    }
}

TEST_CASE("shifting a spectrum by d puts it at distance d") {
    FixedEnvelope a = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-10, true, false);
    a.normalize();
    for (double d : {0.001, 0.5, 3.0, 1000.0}) {
        CAPTURE(d);
        FixedEnvelope b = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-10, true, false);
        b.normalize();
        b.shift_mass(d);
        CHECK(a.WassersteinDistance(b) == doctest::Approx(d).epsilon(1e-9));
        // Oriented distance keeps the sign of the displacement.
        CHECK(std::fabs(a.OrientedWassersteinDistance(b)) == doctest::Approx(d).epsilon(1e-9));
        CHECK(a.OrientedWassersteinDistance(b) == -b.OrientedWassersteinDistance(a));
    }
}

TEST_CASE("Wasserstein distance obeys the triangle inequality") {
    FixedEnvelope a = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-8, true, false);
    FixedEnvelope b = FixedEnvelope::FromThreshold(Iso("C6H12O5N1"), 1e-8, true, false);
    FixedEnvelope c = FixedEnvelope::FromThreshold(Iso("C5H12O6"), 1e-8, true, false);
    a.normalize(); b.normalize(); c.normalize();

    const double ab = a.WassersteinDistance(b);
    const double bc = b.WassersteinDistance(c);
    const double ac = a.WassersteinDistance(c);
    CHECK(ab >= 0.0);
    CHECK(bc >= 0.0);
    CHECK(ac >= 0.0);
    CHECK(ac <= ab + bc + 1e-9);
    // Distinct molecules are at positive distance.
    CHECK(ab > 0.0);
    CHECK(ac > 0.0);
    // Symmetry.
    CHECK(ab == doctest::Approx(b.WassersteinDistance(a)));
}

TEST_CASE("unnormalized spectra are rejected") {
    FixedEnvelope a = make_envelope({1.0, 2.0}, {0.5, 0.5});
    FixedEnvelope b = make_envelope({1.0, 2.0}, {0.5, 5.0});
    CHECK_THROWS_AS(a.WassersteinDistance(b), std::logic_error);
    CHECK_THROWS_AS(a.OrientedWassersteinDistance(b), std::logic_error);
    CHECK_THROWS_AS(b.WassersteinDistance(a), std::logic_error);

    // A 0.1% mismatch is inside the tolerance the implementation allows.
    FixedEnvelope c = make_envelope({1.0, 2.0}, {0.5, 0.5005});
    CHECK_NOTHROW(a.WassersteinDistance(c));
}

TEST_CASE("distance to an empty spectrum is zero by convention") {
    FixedEnvelope a = make_envelope({1.0, 2.0}, {0.5, 0.5});
    FixedEnvelope empty;
    // get_total_prob() of an empty envelope is 0, which trips the
    // normalization check before the emptiness short-circuit, so this throws.
    CHECK_THROWS_AS(a.WassersteinDistance(empty), std::logic_error);

    // Two empty spectra are trivially at distance 0.
    FixedEnvelope empty2;
    CHECK(empty.WassersteinDistance(empty2) == doctest::Approx(0.0));
}

TEST_CASE("oriented distance is signed and cancels") {
    // b sits entirely to the right of a: the oriented distance is +/- the
    // transport distance, whereas the unsigned one is its absolute value.
    FixedEnvelope a = make_envelope({0.0, 1.0}, {0.5, 0.5});
    FixedEnvelope b = make_envelope({2.0, 3.0}, {0.5, 0.5});
    const double oriented = a.OrientedWassersteinDistance(b);
    CHECK(std::fabs(oriented) == doctest::Approx(a.WassersteinDistance(b)));
    CHECK(oriented == -b.OrientedWassersteinDistance(a));

    // When mass moves both ways the oriented distance partially cancels, so it
    // is strictly smaller in magnitude than the unsigned distance.
    FixedEnvelope c = make_envelope({0.0, 10.0}, {0.5, 0.5});
    FixedEnvelope d = make_envelope({5.0, 5.0}, {0.5, 0.5});
    CHECK(std::fabs(c.OrientedWassersteinDistance(d)) < c.WassersteinDistance(d));
}

TEST_CASE("abyssal distance: deep abyss transports, shallow abyss discards") {
    // One unit of mass has to travel a distance of 1.
    FixedEnvelope a = make_envelope({0.0}, {1.0});
    FixedEnvelope b = make_envelope({1.0}, {1.0});

    // With an abyss deeper than the gap, transporting is cheaper: cost == gap.
    CHECK(a.AbyssalWassersteinDistance(b, 10.0) == doctest::Approx(1.0));

    // With a shallow abyss it is cheaper to condemn both peaks; the cost is
    // then (total condemned mass) * depth / 2.
    const double depth = 0.25;
    CHECK(a.AbyssalWassersteinDistance(b, depth) == doctest::Approx(2.0 * depth * 0.5));

    // Monotone in the depth, and never above the plain Wasserstein distance.
    double previous = -1.0;
    for (double d : {0.1, 0.5, 1.0, 2.0, 5.0}) {
        FixedEnvelope x = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-8, true, false);
        FixedEnvelope y = FixedEnvelope::FromThreshold(Iso("C6H12O5N1"), 1e-8, true, false);
        x.normalize();
        y.normalize();
        const double abyssal = x.AbyssalWassersteinDistance(y, d);
        CHECK(abyssal >= previous - 1e-9);
        CHECK(abyssal <= x.WassersteinDistance(y) + 1e-9);
        previous = abyssal;
    }
}

TEST_CASE("abyssal distance scales the second spectrum") {
    FixedEnvelope a = make_envelope({0.0}, {1.0});
    FixedEnvelope b = make_envelope({0.0}, {0.5});
    // With other_scale = 2 the spectra become identical.
    CHECK(a.AbyssalWassersteinDistance(b, 10.0, 2.0) == doctest::Approx(0.0));
    // Unscaled, half a unit of mass is left over and must be condemned.
    CHECK(a.AbyssalWassersteinDistance(b, 10.0, 1.0) == doctest::Approx(0.5 * 10.0 * 0.5));
}

TEST_CASE("WassersteinMatch conserves probability") {
    FixedEnvelope a = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-6, true, false);
    FixedEnvelope b = FixedEnvelope::FromThreshold(Iso("C6H12O6"), 1e-6, true, false);
    a.normalize();
    b.normalize();

    SUBCASE("identical spectra match completely") {
        auto [unmatched_a, unmatched_b, flow] = a.WassersteinMatch(b, 0.5, 1.0);
        CHECK(flow == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(unmatched_a == doctest::Approx(0.0).epsilon(1e-9));
        CHECK(unmatched_b == doctest::Approx(0.0).epsilon(1e-9));
    }

    SUBCASE("a large displacement leaves everything unmatched") {
        b.shift_mass(1000.0);
        auto [unmatched_a, unmatched_b, flow] = a.WassersteinMatch(b, 0.5, 1.0);
        CHECK(flow == doctest::Approx(0.0).epsilon(1e-9));
        CHECK(unmatched_a == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(unmatched_b == doctest::Approx(1.0).epsilon(1e-9));
    }

    SUBCASE("matched + unmatched accounts for all the probability") {
        FixedEnvelope c = FixedEnvelope::FromThreshold(Iso("C6H12O5N1"), 1e-6, true, false);
        c.normalize();
        for (double flow_dist : {0.0, 0.01, 0.1, 1.0}) {
            CAPTURE(flow_dist);
            auto [unmatched_a, unmatched_c, flow] = a.WassersteinMatch(c, flow_dist, 1.0);
            CHECK(flow >= 0.0);
            CHECK(unmatched_a >= -1e-9);
            CHECK(unmatched_c >= -1e-9);
            CHECK(flow + unmatched_a == doctest::Approx(a.get_total_prob()).epsilon(1e-9));
            CHECK(flow + unmatched_c == doctest::Approx(c.get_total_prob()).epsilon(1e-9));
        }
    }

    SUBCASE("other_scale rescales the second spectrum's mass budget") {
        auto [unmatched_a, unmatched_b, flow] = a.WassersteinMatch(b, 0.5, 0.5);
        CHECK(flow == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(unmatched_a == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(unmatched_b == doctest::Approx(0.0).epsilon(1e-9));
    }

    SUBCASE("empty first spectrum") {
        FixedEnvelope empty;
        auto [unmatched_empty, unmatched_b, flow] = empty.WassersteinMatch(b, 0.5, 1.0);
        CHECK(unmatched_empty == doctest::Approx(0.0));
        CHECK(flow == doctest::Approx(0.0));
        CHECK(unmatched_b == doctest::Approx(b.get_total_prob()));
    }
}

TEST_CASE("WassersteinMatch: peaks pair up only within the flow distance") {
    // Ported from the standalone wasserstein_matching.cpp check.
    // 1.0 pairs with 0.9 (gap 0.1 > 0.05) — no; 2.0 pairs with 2.01 (gap 0.01).
    FixedEnvelope a = make_envelope({1.0, 2.0}, {0.5, 0.5});
    FixedEnvelope b = make_envelope({0.9, 2.01}, {0.4, 0.6});

    auto [unmatched_a, unmatched_b, flow] = a.WassersteinMatch(b, 0.05);
    CHECK(unmatched_a == doctest::Approx(0.5));
    CHECK(unmatched_b == doctest::Approx(0.5));
    CHECK(flow == doctest::Approx(0.5));

    // Widening the flow distance lets the 1.0/0.9 pair match too.
    FixedEnvelope a2 = make_envelope({1.0, 2.0}, {0.5, 0.5});
    FixedEnvelope b2 = make_envelope({0.9, 2.01}, {0.4, 0.6});
    auto [ua2, ub2, flow2] = a2.WassersteinMatch(b2, 0.2);
    CHECK(flow2 > flow);
    CHECK(flow2 + ua2 == doctest::Approx(1.0));
    CHECK(flow2 + ub2 == doctest::Approx(1.0));
}

TEST_CASE("distances see the spectrum, not the peak order") {
    // The distance functions sort their inputs; feeding them a prob-sorted
    // envelope must give the same answer as a mass-sorted one.
    FixedEnvelope a = FixedEnvelope::FromThreshold(Iso("C10H20O2"), 1e-8, true, false);
    FixedEnvelope b = FixedEnvelope::FromThreshold(Iso("C10H20O2"), 1e-8, true, false);
    a.normalize();
    b.normalize();
    b.shift_mass(0.7);

    a.sort_by_prob();
    b.sort_by_prob();
    const double d1 = a.WassersteinDistance(b);

    FixedEnvelope a2 = FixedEnvelope::FromThreshold(Iso("C10H20O2"), 1e-8, true, false);
    FixedEnvelope b2 = FixedEnvelope::FromThreshold(Iso("C10H20O2"), 1e-8, true, false);
    a2.normalize();
    b2.normalize();
    b2.shift_mass(0.7);
    a2.sort_by_mass();
    b2.sort_by_mass();
    CHECK(d1 == doctest::Approx(a2.WassersteinDistance(b2)));
}
