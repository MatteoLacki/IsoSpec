// isoMath.h/.cpp: the special functions and random variates the algorithms are
// built on.
//
// Each is checked against an independent reference — a textbook identity, a
// numerically integrated definition, or the statistics of a large sample.

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "doctest.h"
#include "isoMath.h"

using namespace IsoSpec;

TEST_CASE("minuslogFactorial equals -lgamma(n+1)") {
    // Both sides of the eager-fill boundary (ISOSPEC_LFACT_EAGER_FILL) and the
    // n < 2 short-circuit.
    const int points[] = {0, 1, 2, 3, 10, 100, 1000, 65535, 65536, 65537, 100000, 1000000};
    for (int n : points) {
        CAPTURE(n);
        CHECK(minuslogFactorial(n) == doctest::Approx(-std::lgamma(n + 1.0)).epsilon(1e-12));
    }
    CHECK(minuslogFactorial(0) == 0.0);
    CHECK(minuslogFactorial(1) == 0.0);

    // Repeated reads of a lazily-filled slot must be stable.
    const double first = minuslogFactorial(200000);
    CHECK(minuslogFactorial(200000) == first);
}

TEST_CASE("NormalCDF matches known quantiles") {
    // Standard normal reference values.
    CHECK(NormalCDF(0.0, 0.0, 1.0) == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(NormalCDF(1.0, 0.0, 1.0) == doctest::Approx(0.8413447).epsilon(1e-6));
    CHECK(NormalCDF(-1.0, 0.0, 1.0) == doctest::Approx(0.1586553).epsilon(1e-6));
    CHECK(NormalCDF(1.959964, 0.0, 1.0) == doctest::Approx(0.975).epsilon(1e-6));

    // Symmetry, monotonicity and the limits.
    for (double x = -4.0; x <= 4.0; x += 0.25) {
        CAPTURE(x);
        CHECK(NormalCDF(x, 0.0, 1.0) + NormalCDF(-x, 0.0, 1.0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(NormalCDF(x, 0.0, 1.0) < NormalCDF(x + 0.25, 0.0, 1.0));
    }
    CHECK(NormalCDF(-40.0, 0.0, 1.0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(NormalCDF(40.0, 0.0, 1.0) == doctest::Approx(1.0).epsilon(1e-9));

    // Location and scale.
    CHECK(NormalCDF(10.0, 10.0, 3.0) == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(NormalCDF(13.0, 10.0, 3.0) == doctest::Approx(NormalCDF(1.0, 0.0, 1.0)).epsilon(1e-6));
}

TEST_CASE("NormalPDF is the derivative of the CDF") {
    for (double x : {-2.0, -0.5, 0.0, 0.5, 2.0}) {
        CAPTURE(x);
        const double h = 1e-4;
        const double numeric = (NormalCDF(x + h, 0.0, 1.0) - NormalCDF(x - h, 0.0, 1.0)) / (2 * h);
        CHECK(NormalPDF(x) == doctest::Approx(numeric).epsilon(1e-4));
    }

    // Closed form, and the defaults.
    CHECK(NormalPDF(0.0) == doctest::Approx(1.0 / std::sqrt(2 * pi)));
    CHECK(NormalPDF(0.0, 0.0, 1.0) == NormalPDF(0.0));
    CHECK(NormalPDF(1.0, 1.0, 2.0) == doctest::Approx(1.0 / (2.0 * std::sqrt(2 * pi))));
    // Symmetric about the mean.
    CHECK(NormalPDF(3.0, 5.0, 2.0) == doctest::Approx(NormalPDF(7.0, 5.0, 2.0)));
}

TEST_CASE("NormalCDFInverse inverts NormalCDF") {
    // The rational approximation is documented as good to ~4.5e-4.
    for (double p : {0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999}) {
        CAPTURE(p);
        const double x = NormalCDFInverse(p);
        CHECK(NormalCDF(x, 0.0, 1.0) == doctest::Approx(p).epsilon(1e-3));
    }
    CHECK(NormalCDFInverse(0.5) == doctest::Approx(0.0).epsilon(1e-3));
    CHECK(NormalCDFInverse(0.975) == doctest::Approx(1.959964).epsilon(1e-3));
    // Antisymmetry about p = 0.5, and monotonicity.
    CHECK(NormalCDFInverse(0.1) == doctest::Approx(-NormalCDFInverse(0.9)).epsilon(1e-9));
    CHECK(NormalCDFInverse(0.2) < NormalCDFInverse(0.8));

    // Location/scale overload.
    CHECK(NormalCDFInverse(0.5, 10.0, 3.0) == doctest::Approx(10.0).epsilon(1e-3));
    CHECK(NormalCDFInverse(0.975, 10.0, 3.0) ==
          doctest::Approx(10.0 + 3.0 * NormalCDFInverse(0.975)).epsilon(1e-9));
}

TEST_CASE("LowerIncompleteGamma2 matches closed forms") {
    // LowerIncompleteGamma2(a, x) is the lower incomplete gamma of (a/2, x).
    // For a = 2 (i.e. s = 1) it is 1 - exp(-x).
    for (double x : {0.1, 0.5, 1.0, 2.0, 5.0}) {
        CAPTURE(x);
        CHECK(LowerIncompleteGamma2(2, x) == doctest::Approx(1.0 - std::exp(-x)).epsilon(1e-12));
        // For a = 1 (s = 1/2) it is sqrt(pi) * erf(sqrt(x)).
        CHECK(LowerIncompleteGamma2(1, x) ==
              doctest::Approx(std::sqrt(pi) * std::erf(std::sqrt(x))).epsilon(1e-12));
        // For a = 4 (s = 2) it is 1 - (1+x)exp(-x).
        CHECK(LowerIncompleteGamma2(4, x) ==
              doctest::Approx(1.0 - (1.0 + x) * std::exp(-x)).epsilon(1e-12));
    }

    // Monotone increasing in x, tending to gamma(a/2).
    for (int a : {1, 2, 3, 6}) {
        CAPTURE(a);
        double previous = -1.0;
        for (double x = 0.5; x < 20.0; x += 0.5) {
            const double v = LowerIncompleteGamma2(a, x);
            CHECK(v > previous);
            previous = v;
        }
        CHECK(previous == doctest::Approx(std::tgamma(a / 2.0)).epsilon(1e-4));
    }
}

TEST_CASE("InverseLowerIncompleteGamma2 inverts LowerIncompleteGamma2") {
    // Only for `a` large enough that the root falls inside the bisection
    // bracket the implementation uses; see the DOC case at the end of the file.
    for (int a : {4, 6, 8, 12}) {
        for (double target : {0.1, 0.5, 0.9}) {
            CAPTURE(a);
            CAPTURE(target);
            const double x = target * std::tgamma(a / 2.0);
            const double y = InverseLowerIncompleteGamma2(a, x);
            // The bisection stops at a relative width of 1e-3.
            CHECK(LowerIncompleteGamma2(a, y) == doctest::Approx(x).epsilon(1e-2));
        }
    }
}

TEST_CASE("InverseChiSquareCDF2 matches tabulated chi-square quantiles") {
    // Textbook chi-square quantiles (k degrees of freedom, p quantile).  Only
    // the regime where the underlying bisection brackets the root is checked
    // here — see the DOC case below for the rest.
    struct Ref { int k; double p; double q; };
    const Ref refs[] = {
        {4,  0.50,  3.357},
        {5,  0.50,  4.351},
        {5,  0.95, 11.070},
        {6,  0.95, 12.592},
        {8,  0.95, 15.507},
        {10, 0.95, 18.307},
        {10, 0.50,  9.342},
        {10, 0.99, 23.209},
        {20, 0.95, 31.410},
        {20, 0.50, 19.337},
    };
    for (const Ref& r : refs) {
        CAPTURE(r.k);
        CAPTURE(r.p);
        // The bisection inside stops at 0.1% relative width, so allow 1%.
        CHECK(InverseChiSquareCDF2(r.k, r.p) == doctest::Approx(r.q).epsilon(1e-2));
    }

    // Monotone in both arguments.
    for (int k : {5, 10, 20}) {
        CAPTURE(k);
        CHECK(InverseChiSquareCDF2(k, 0.5) < InverseChiSquareCDF2(k, 0.9));
        CHECK(InverseChiSquareCDF2(k, 0.9) < InverseChiSquareCDF2(k, 0.99));
    }
    for (double p : {0.5, 0.9, 0.99}) {
        CAPTURE(p);
        CHECK(InverseChiSquareCDF2(5, p) < InverseChiSquareCDF2(20, p));
    }
}

TEST_CASE("DOC: InverseLowerIncompleteGamma2 saturates for small `a`") {
    // The bisection in InverseLowerIncompleteGamma2 searches the interval
    // [0, tgamma(a)] — but tgamma(a) is a value in the function's *range*, not
    // a bound on its *argument*.  Whenever the root lies above tgamma(a) the
    // bisection cannot reach it and converges to the bracket's upper end
    // instead, silently returning tgamma(a).
    //
    // Consequence: InverseChiSquareCDF2(k, p) is wrong for small k — e.g.
    // chi-square(1) at p=0.95 is 3.841 but the call returns 2*tgamma(1) = 2 —
    // and, through Iso::saveMarginalLogSizeEstimates, the marginal size
    // estimates saturate too.  That only mis-ranks marginals (a heuristic:
    // it costs speed, not correctness of the computed distribution), which is
    // why this is documented rather than asserted.
    CHECK(InverseLowerIncompleteGamma2(1, 0.95 * std::tgamma(0.5)) ==
          doctest::Approx(std::tgamma(1)).epsilon(1e-3));
    CHECK(InverseChiSquareCDF2(1, 0.95) == doctest::Approx(2.0).epsilon(1e-3));
    CHECK(InverseChiSquareCDF2(2, 0.95) == doctest::Approx(2.0).epsilon(1e-3));
    // Saturated: the answer no longer even depends on p.
    CHECK(InverseChiSquareCDF2(2, 0.95) == doctest::Approx(InverseChiSquareCDF2(2, 0.99)));

    MESSAGE("InverseLowerIncompleteGamma2 bisects [0, tgamma(a)] and saturates "
            "when the root exceeds tgamma(a); fix by growing the upper bound "
            "until LowerIncompleteGamma2(a, hi) >= x before bisecting.");
}

TEST_CASE("rdvariate_beta_1_b samples Beta(1, b)") {
    std::mt19937 gen(12345);
    for (double b : {1.0, 2.0, 10.0}) {
        CAPTURE(b);
        const int n = 200000;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            const double x = rdvariate_beta_1_b(b, gen);
            REQUIRE(x >= 0.0);
            REQUIRE(x <= 1.0);
            sum += x;
        }
        // E[Beta(1,b)] = 1/(1+b).
        CHECK(sum / n == doctest::Approx(1.0 / (1.0 + b)).epsilon(2e-2));
    }
}

TEST_CASE("rdvariate_binom samples Binomial(n, p)") {
    std::mt19937 gen(6789);

    SUBCASE("degenerate success probabilities") {
        CHECK(rdvariate_binom(100, 1.0, gen) == 100);
        CHECK(rdvariate_binom(100, 2.0, gen) == 100);  // >= 1 is clamped
        CHECK(rdvariate_binom(100, 0.0, gen) == 0);
        CHECK(rdvariate_binom(0, 0.5, gen) == 0);
    }

    SUBCASE("mean and variance over many draws") {
        // Small n uses the inversion path, large n the BTRD algorithm (btrd.h).
        struct Case { std::size_t n; double p; };
        const Case cases[] = {
            {1, 0.5}, {10, 0.1}, {30, 0.5}, {1000, 0.01}, {1000, 0.5}, {1000000, 0.3},
        };
        for (const Case& c : cases) {
            CAPTURE(c.n);
            CAPTURE(c.p);
            const int draws = 20000;
            double sum = 0.0, sumsq = 0.0;
            for (int i = 0; i < draws; ++i) {
                const std::size_t k = rdvariate_binom(c.n, c.p, gen);
                REQUIRE(k <= c.n);
                sum += static_cast<double>(k);
                sumsq += static_cast<double>(k) * static_cast<double>(k);
            }
            const double mean = sum / draws;
            const double var = sumsq / draws - mean * mean;
            const double expected_mean = c.n * c.p;
            const double expected_var = c.n * c.p * (1 - c.p);
            // Sample mean is within 5 standard errors of the true mean.
            const double se = std::sqrt(expected_var / draws);
            CHECK(std::fabs(mean - expected_mean) < 5 * se + 1e-9);
            CHECK(var == doctest::Approx(expected_var).epsilon(0.1));
        }
    }

    SUBCASE("symmetry: Binom(n, p) and n - Binom(n, 1-p) agree in distribution") {
        const std::size_t n = 1000;
        const int draws = 20000;
        double sum_a = 0.0, sum_b = 0.0;
        for (int i = 0; i < draws; ++i) {
            sum_a += static_cast<double>(rdvariate_binom(n, 0.3, gen));
            sum_b += static_cast<double>(n - rdvariate_binom(n, 0.7, gen));
        }
        CHECK(sum_a / draws == doctest::Approx(sum_b / draws).epsilon(2e-2));
    }
}

TEST_CASE("the thread-local generators are usable directly") {
    // random_gen / stdunif back the default arguments of the variate helpers.
    for (int i = 0; i < 1000; ++i) {
        const double u = stdunif(random_gen);
        REQUIRE(u >= 0.0);
        REQUIRE(u < 1.0);
    }
    const double b = rdvariate_beta_1_b(3.0);
    CHECK(b >= 0.0);
    CHECK(b <= 1.0);
    CHECK(rdvariate_binom(10, 0.5) <= 10u);
}
