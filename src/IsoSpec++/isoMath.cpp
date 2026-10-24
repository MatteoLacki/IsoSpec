/*
 *   This file has been released into public domain by John D. Cook
 *   and is used here with some slight modifications (which are hereby
 *   also released into public domain),
 *
 *   This file is part of IsoSpec.
 */


// NOLINT(legal/copyright)


#include <cmath>
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

double LowerIncompleteGamma2(int a, double x)
{
    double base;
    double exp_minus_x = exp(-x);
    double current_s;
    if(a % 2 == 0)
    {
        base = 1 - exp_minus_x;
        current_s = 1.0;
        a--;
    }
    else
    {
        base = sqrt_pi * erf(sqrt(x));
        current_s = 0.5;
    }

    a = a/2;
    for(; a; a--)
    {
        base = base * current_s - pow(x, current_s) * exp_minus_x;
        current_s += 1.0;
    }

    return base;
}

double InverseLowerIncompleteGamma2(int a, double x)
{
    double l = 0.0;
    double p = tgamma(a);
    double s;

    do {
        s = (l+p) / 2.0;
        double v = LowerIncompleteGamma2(a, s);
        if (x < v)
            p = s;
        else
            l = s;
    } while((p-l)*1000.0 > p);

    return s;
}

thread_local std::mt19937 random_gen(std::random_device{}());

size_t rdvariate_binom(size_t tries, double succ_prob, std::mt19937& rgen)
{
    if (succ_prob >= 1.0)
        return tries;
    return IsoSpec::boost_binomial_distribution_variate(tries, succ_prob, rgen);
}



}  // namespace IsoSpec

