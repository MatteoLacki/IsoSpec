// isoMath.h/.cpp: the special functions and random variates the algorithms are
// built on.
//
// Each is checked against an independent reference — a textbook identity, a
// numerically integrated definition, or the statistics of a large sample.

#include <algorithm>
#include <cmath>
#include <limits>
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

    // Past the end of the table: computed, not cached.
    const int beyond = ISOSPEC_G_FACT_TABLE_SIZE + 5;
    CHECK(minuslogFactorial(beyond) == doctest::Approx(-std::lgamma(beyond + 1.0)).epsilon(1e-14));
}

TEST_CASE("lgamma_half and lgamma_stirling equal lgamma at integers and half-integers") {
    // Below 8 compiled-in tables, from 8 Stirling's series; every value within
    // a few ulp of the platform's lgamma, itself within a few ulp of exact.
    for (std::int64_t two_x = 2; two_x <= 4000; ++two_x) {
        CAPTURE(two_x);
        CHECK(lgamma_half(two_x) == doctest::Approx(std::lgamma(two_x * 0.5)).epsilon(1e-14));
    }
    const std::int64_t distant[] = {65535, 65536, 65537, 1999999, 2000000, 20000001,
                                std::int64_t{2} * ISOSPEC_G_FACT_TABLE_SIZE + 7};
    for (std::int64_t two_x : distant) {
        CAPTURE(two_x);
        CHECK(lgamma_half(two_x) == doctest::Approx(std::lgamma(two_x * 0.5)).epsilon(1e-14));
    }
    CHECK(lgamma_half(2) == 0.0);
    CHECK(lgamma_half(4) == 0.0);
    for (double x : {8.0, 8.25, 9.75, 1e3 + 0.1, 1e9, 1e15}) {
        CAPTURE(x);
        CHECK(lgamma_stirling(x) == doctest::Approx(std::lgamma(x)).epsilon(1e-14));
    }
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

namespace {
struct ChiSquareQuantile { int k; double x, q; };
// SciPy's chi-square quantiles (k up to 10^6, x from 1e-300 to 1 - 2^-53),
// written by gen_chi2_quantiles.py.
const ChiSquareQuantile chi2_quantiles[] = {
#include "chi2_quantiles.inc"
};
}  // namespace

TEST_CASE("InverseChiSquareCDF2 agrees with SciPy's chi-square quantiles") {
    // The callers only need about 0.1% (the quantile orders marginals by
    // estimated size); the Halley iteration gets to rounding, and this holds
    // it there. Rounding grows with k: the density's exponent has terms of
    // order k log k, so at k = 10^6 a libm differing by an ulp moves the
    // answer by ~1e-11 (glibc's worst row is 8.6e-12 there; below k = 1000,
    // 6e-14).
    REQUIRE(sizeof(chi2_quantiles) / sizeof(chi2_quantiles[0]) > 2000);
    double worst = 0.0;
    for (const ChiSquareQuantile& r : chi2_quantiles) {
        const double q = InverseChiSquareCDF2(r.k, r.x);
        const double error = std::fabs(q / r.q - 1.0);
        if (error > (r.k <= 1000 ? 1e-11 : 1e-9)) {
            CAPTURE(r.k);
            CAPTURE(r.x);
            CHECK(q == r.q);
        }
        worst = std::max(worst, error);
    }
    MESSAGE("worst relative error against SciPy: " << worst);
}

TEST_CASE("InverseChiSquareCDF2 matches closed forms") {
    for (double x : {1e-9, 0.01, 0.3, 0.5, 0.7, 0.99, 0.9999, 1.0 - 1e-12}) {
        CAPTURE(x);
        // k = 2: exponential with mean 2.
        CHECK(InverseChiSquareCDF2(2, x) == doctest::Approx(-2.0 * std::log1p(-x)).epsilon(1e-12));
        // k = 1: the square of a standard normal, P(chi2 <= q) = erf(sqrt(q/2)).
        CHECK(std::erf(std::sqrt(InverseChiSquareCDF2(1, x) / 2.0)) == doctest::Approx(x).epsilon(1e-12));
        // k = 4: P = 1 - (1 + q/2) exp(-q/2).
        const double h = InverseChiSquareCDF2(4, x) / 2.0;
        CHECK(-std::expm1(-h) - h * std::exp(-h) == doctest::Approx(x).epsilon(1e-12));
    }
}

TEST_CASE("InverseChiSquareCDF2 is strictly increasing in both arguments") {
    for (int k : {1, 2, 3, 4, 5, 8, 13, 40, 299, 300, 301, 5000}) {
        CAPTURE(k);
        double previous = 0.0;
        for (int i = 1; i < 2000; ++i) {
            const double x = i / 2000.0;
            const double q = InverseChiSquareCDF2(k, x);
            CAPTURE(x);
            REQUIRE(q > previous);
            previous = q;
        }
    }
    for (double x : {1e-6, 0.5, 0.99, 0.9999}) {
        CAPTURE(x);
        double previous = 0.0;
        for (int k = 1; k <= 1000; ++k) {
            const double q = InverseChiSquareCDF2(k, x);
            CAPTURE(k);
            REQUIRE(q > previous);
            previous = q;
        }
    }
}

TEST_CASE("InverseChiSquareCDF2 is defined for every argument") {
    // Its value becomes a sort key (Iso::saveMarginalLogSizeEstimates): it
    // must never be NaN, or std::sort's comparator stops being an ordering.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (int k : {-5, 0}) {
        CHECK(InverseChiSquareCDF2(k, 0.5) == 0.0);
        CHECK(InverseChiSquareCDF2(k, nan) == 0.0);
    }
    for (int k : {1, 2, 3, 8, 100, 1000000}) {
        CAPTURE(k);
        CHECK(InverseChiSquareCDF2(k, 0.0) == 0.0);
        CHECK(InverseChiSquareCDF2(k, -1.0) == 0.0);
        CHECK(InverseChiSquareCDF2(k, -inf) == 0.0);
        CHECK(InverseChiSquareCDF2(k, nan) == 0.0);
        // The quantile at 1 is infinite; the answer is the one just below.
        const double top = InverseChiSquareCDF2(k, 1.0 - 0x1p-53);
        CHECK(std::isfinite(top));
        CHECK(InverseChiSquareCDF2(k, 1.0) == top);
        CHECK(InverseChiSquareCDF2(k, 2.0) == top);
        CHECK(InverseChiSquareCDF2(k, inf) == top);
        // Down to the smallest doubles, where the quantile underflows.
        for (double x : {1e-300, 0x1p-1022, 0x1p-1074}) {
            const double q = InverseChiSquareCDF2(k, x);
            CHECK(std::isfinite(q));
            CHECK(q >= 0.0);
        }
    }
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
    // random_gen backs the default arguments of the variate helpers.
    for (int i = 0; i < 1000; ++i) {
        const double u = rdvariate_unif01(random_gen);
        REQUIRE(u >= 0.0);
        REQUIRE(u < 1.0);
    }
    const double b = rdvariate_beta_1_b(3.0);
    CHECK(b >= 0.0);
    CHECK(b <= 1.0);
    CHECK(rdvariate_binom(10, 0.5) <= 10u);
}
