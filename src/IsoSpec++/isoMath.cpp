/*
 *   This file has been released into public domain by John D. Cook
 *   and is used here with some slight modifications (which are hereby
 *   also released into public domain),
 *
 *   This file is part of IsoSpec.
 */


// NOLINT(legal/copyright)


#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <cstdlib>
#include <new>
#include "isoMath.h"
#include "platform.h"
#include "btrd.h"

namespace IsoSpec
{

/*
 * The body of lgamma_stirling() is the 8.0 <= x < 2**58 branch of fdlibm's
 * e_lgamma_r.c, under this notice:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 *
 * lnGamma(x) = (x-0.5)*(log(x)-1) + w, where w approximates
 * 0.5*ln(2*pi) - 0.5 + (Stirling's correction series in 1/x).
 * fdlibm switches to x*(log(x)-1) from 2**58 up; there w is below an
 * ulp of the result, so the one formula serves every argument here.
 * Against correctly rounded values: at most 2 ulp off, on every integer and
 * half-integer from 8 to 10**7.
 */
double lgamma_stirling(double x) noexcept
{
    ISOSPEC_IMPOSSIBLE(x < 8.0);
    constexpr double w0 =  4.18938533204672725052e-01;
    constexpr double w1 =  8.33333333333329678849e-02;
    constexpr double w2 = -2.77777777728775536470e-03;
    constexpr double w3 =  7.93650558643019558500e-04;
    constexpr double w4 = -5.95187557450339963135e-04;
    constexpr double w5 =  8.36339918996282139126e-04;
    constexpr double w6 = -1.63092934096575273989e-03;
    const double t = std::log(x);
    const double z = 1.0 / x;
    const double y = z * z;
    const double w = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))));
    return (x - 0.5) * (t - 1.0) + w;
}

double lgamma_half(std::int64_t two_x) noexcept
{
    ISOSPEC_IMPOSSIBLE(two_x < 2);
    if (two_x % 2 == 0)
    {
        const std::int64_t n = two_x / 2 - 1;      // lnGamma(n+1) = ln(n!)
        if (n < ISOSPEC_G_FACT_TABLE_SIZE)
            return 0.0 - minuslogFactorial(static_cast<int>(n));    // +0.0, not -0.0, at x = 1, 2
        return lgamma_stirling(static_cast<double>(n) + 1.0);
    }
    // lnGamma(1.5), lnGamma(2.5), ..., lnGamma(7.5), correctly rounded.
    static constexpr double small_half_integers[7] = {
        -0x1.eeb95b094c191p-4, 0x1.2383e809a67e8p-2, 0x1.3373018970a36p+0, 0x1.3a140a3a623cbp+1,
        0x1.fa99a5e94985ap+1, 0x1.6a676ab2aa095p+2, 0x1.e233060e41f7fp+2};
    if (two_x < 16)
        return small_half_integers[(two_x - 3) / 2];
    return lgamma_stirling(static_cast<double>(two_x) * 0.5);
}


double* alloc_lfact_table()
{
    double* ret;
# if ISOSPEC_GOT_MMAN
    ret = reinterpret_cast<double*>(mmap(nullptr, sizeof(double)*ISOSPEC_G_FACT_TABLE_SIZE, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
    if(ret == MAP_FAILED)
        throw std::bad_alloc();
#else
    ret = reinterpret_cast<double*>(calloc(ISOSPEC_G_FACT_TABLE_SIZE, sizeof(double)));
    if(ret == nullptr)
        throw std::bad_alloc();
#endif

    // Eager-fill the bottom of the table.  Done before any thread can call
    // minuslogFactorial(), so the eager region is immutable post-init and the
    // hot path reads it with no synchronisation.  The lazy tail (if any) uses
    // std::atomic_ref to handle racing fills safely.
    constexpr std::size_t eager_end =
        static_cast<std::size_t>(ISOSPEC_G_FACT_TABLE_SIZE) <
                static_cast<std::size_t>(ISOSPEC_LFACT_EAGER_FILL)
            ? static_cast<std::size_t>(ISOSPEC_G_FACT_TABLE_SIZE)
            : static_cast<std::size_t>(ISOSPEC_LFACT_EAGER_FILL);
    // Below 7, n+1 is under Stirling's crossover at 8: -ln(n!), correctly
    // rounded. Written even when ISOSPEC_LFACT_EAGER_FILL is smaller, which is
    // harmless: the lazy path only computes a slot it finds still 0.0.
    static_assert(ISOSPEC_G_FACT_TABLE_SIZE >= 7, "log-factorial table must hold the compiled-in values");
    static constexpr double small_minuslogfactorials[7] = {
        0.0, 0.0, -0x1.62e42fefa39efp-1, -0x1.cab0bfa2a2002p+0, -0x1.96ca77c922cf9p+1,
        -0x1.326643c4479c9p+2, -0x1.a51273acf01cap+2};
    for (std::size_t n = 0; n < 7; ++n)
        ret[n] = small_minuslogfactorials[n];
    for (std::size_t n = 7; n < eager_end; ++n)
        ret[n] = -lgamma_stirling(static_cast<double>(n + 1));

    // Deliberately never freed. An atexit() handler releasing it would run
    // while threads that outlive main() -- detached, or still winding down --
    // may be reading it, a use-after-free at exit; the OS reclaims it at
    // process end anyway. The cost: a library that is dlclose()d and loaded
    // again leaks one table per reload, mostly untouched address space past
    // the eager prefix.
    return ret;
}

double* g_lfact_table = alloc_lfact_table();


double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p < 0.5)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    else
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
}

double NormalCDFInverse(double p, double mean, double stdev)
{
    return mean + stdev * NormalCDFInverse(p);
}

double NormalCDF(double x, double mean, double stdev)
{
    x = (x-mean)/stdev * 0.7071067811865476;

    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

double NormalPDF(double x, double mean, double stdev)
{
    double two_variance = stdev * stdev * 2.0;
    double delta = x-mean;
    return exp( -delta*delta / two_variance )      /     sqrt( two_variance * pi );
}

const double sqrt_pi = 1.772453850905516027298167483341145182798;

namespace
{
// The regularised incomplete gamma functions at a = k/2, k >= 1: their values
// P(a, s) and Q(a, s) = 1 - P(a, s), each computed directly so that neither
// loses its relative precision to cancellation in its tail, and the density
// dP/ds = s^(a-1) e^(-s) / Gamma(a), which every Halley step needs anyway.
// `lgamma_a` is lnGamma(a). Region by region, as Boost's gamma_p does it:
//   s <  a+1:  P = s^a e^-s / Gamma(a+1) * sum_n s^n / ((a+1)...(a+n)),
//              whose prefix is density * s / a and whose terms fall at once;
//   s >= a+1:  Q's finite sum, for k = 2m and k = 2m+1
//              Q(m, s)     = e^-s * sum_{j<m} s^j / j!
//              Q(m+1/2, s) = erfc(sqrt(s)) + e^-s * sum_{j<m} s^(j+1/2) / Gamma(j+3/2),
//              summed downward from its last term, which is the density:
//              each term is the next one times (j+a-m)/s < 1, so nothing
//              overflows whatever k, and the sum stops once the terms are
//              below rounding.
struct GammaHalf
{
    double P, Q, density;
};

GammaHalf gamma_half(int k, double s, double lgamma_a)
{
    const double a = 0.5 * k;
    const double density = exp((a - 1.0) * log(s) - s - lgamma_a);

    if(s < a + 1.0)
    {
        double t = 1.0, series = 1.0;
        for(double d = a + 1.0; t > series * 0x1p-54; d += 1.0)
        {
            t *= s / d;
            series += t;
        }
        const double P = density * (s / a) * series;
        return {P, 1.0 - P, density};
    }

    double sum = 0.0, t = density;
    for(int i = 0; i < k / 2 && t > sum * 0x1p-54; i++)
    {
        sum += t;
        t *= (a - 1.0 - i) / s;
    }
    const double Q = (k % 2 == 0) ? sum : erfc(sqrt(s)) + sum;
    return {1.0 - Q, Q, density};
}
}  // namespace

double InverseChiSquareCDF2(int k, double x)
{
    // q = 2s with P(k/2, s) = x. The initial guess is the larger of
    //   * Wilson and Hilferty (1931): (q/k)^(1/3) is close to normal with mean
    //     1 - 2/(9k) and variance 2/(9k), so q ~ k (1 - 2/(9k) + z sqrt(2/(9k)))^3
    //     with z the normal quantile of x -- within a fraction of a percent
    //     away from the far tails;
    //   * the leading term of P's series, P(a, s) <= s^a / Gamma(a+1), whose
    //     inverse s = (x Gamma(a+1))^(1/a) is a lower bound on the root
    //     everywhere and accurate in the lower tail, where Wilson-Hilferty's
    //     cube can go negative.
    // Then Halley's iteration (as in DiDonato and Morris 1986, and Boost's
    // gamma_p_inv): f'' / f' = (a-1)/s - 1 for f = P(a, s) - x, so each step
    // costs one evaluation of gamma_half, and from this guess two (rarely
    // more; at most 10 in the far tails of k < 10) reach full precision. The callers need about 0.1% (the quantile
    // only orders marginals by estimated size), but the steps past that cost
    // little and make the function testable against reference values. No
    // tgamma or pow: one exp and one log per step, lnGamma from lgamma_half.
    //
    // Defined everywhere, so that the sort keys built from it are never NaN:
    // 0 for k <= 0 (no freedom) and for x <= 0 or NaN; x >= 1, whose quantile
    // is infinite, is answered at the largest double below 1.
    if(k <= 0 || !(x > 0.0))
        return 0.0;
    x = (std::min)(x, 1.0 - 0x1p-53);
    if(k == 2)
        return -2.0 * log1p(-x);   // exponential with mean 2: exact

    const double a = 0.5 * k;
    // lnGamma(1/2) = ln(sqrt(pi)) is below lgamma_half's domain.
    const double lgamma_a = k == 1 ? 0.5723649429247001 : lgamma_half(k);

    double s_wh;
    if(k == 1 && x >= 0.5)
    {
        // The square of a standard normal: q = z^2 for z the normal quantile
        // at (1+x)/2, taken from its tail (1-x)/2, which keeps its precision;
        // NormalCDFInverse's 4.5e-4 makes it a guess. Below the median the
        // series bound is the better one.
        const double z = NormalCDFInverse(0.5 * (1.0 - x));
        s_wh = 0.5 * z * z;
    }
    else
    {
        const double v = 2.0 / (9.0 * k);
        const double h = 1.0 - v + NormalCDFInverse(x) * sqrt(v);
        s_wh = h > 0.0 ? a * h * h * h : 0.0;
    }
    const double lower_bound = exp((log(x) + lgamma_half(static_cast<std::int64_t>(k) + 2)) / a);
    double s = (std::max)(s_wh, lower_bound);
    if(!(s > 0.0))
        return 0.0;   // below the smallest double

    // Solved for P below the median and for Q above it, in whichever tail
    // keeps its relative precision; both are increasing in s.
    const bool upper = x > 0.5;
    const double target = upper ? 1.0 - x : x;
    double lo = lower_bound, hi = std::numeric_limits<double>::infinity();
    for(int iteration = 0; iteration < 100; iteration++)
    {
        const GammaHalf g = gamma_half(k, s, lgamma_a);
        const double f = upper ? target - g.Q : g.P - target;
        if(f == 0.0 || !(g.density > 0.0))
            break;
        if(f < 0.0)
            lo = s;
        else
            hi = s;
        const double newton = f / g.density;
        const double curvature = 1.0 - 0.5 * newton * ((a - 1.0) / s - 1.0);
        const double halley = curvature > 0.5 ? newton / curvature : newton;
        // Converged once the step is below 2^-20 relative: Halley's error
        // after a step of d is of order d^3, so taking that step puts s
        // within rounding of the root (checked over k = 1..1000 and x from
        // 1e-15 to 1 - 1e-15 against SciPy: worst 4e-14, one evaluation
        // fewer than a 2^-40 test, the same answers). Checked before the
        // bracket, which a step below one ulp of s cannot enter.
        if(fabs(halley) <= s * 0x1p-20)
        {
            s -= halley;
            break;
        }
        // A step leaving the bracket falls back to Newton's, then to
        // bisection (or doubling, while the bracket has no upper end).
        double next = s - halley;
        if(!(next > lo && next < hi))
            next = s - newton;
        if(!(next > lo && next < hi))
            next = std::isfinite(hi) ? 0.5 * (lo + hi) : 2.0 * (std::max)(s, lo);
        s = next;
    }
    return 2.0 * s;
}


thread_local std::mt19937 random_gen(std::random_device{}());

size_t rdvariate_binom(size_t tries, double succ_prob, std::mt19937& rgen)
{
    if (succ_prob >= 1.0)
        return tries;
    return IsoSpec::boost_binomial_distribution_variate(tries, succ_prob, rgen);
}



}  // namespace IsoSpec

