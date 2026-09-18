/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */

#pragma once

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>

#if !defined(ISOSPEC_G_FACT_TABLE_SIZE)
// 10M should be enough for anyone, right?
// Actually, yes. If anyone tries to input a molecule that has more than 10M atoms,
// he deserves to get an exception thrown in his face. OpenMS guys don't want to alloc
// a table of 10M to memoize the necessary values though, use something smaller for them.
  #if ISOSPEC_BUILDING_OPENMS
    #define ISOSPEC_G_FACT_TABLE_SIZE 1024
  #else
    #define ISOSPEC_G_FACT_TABLE_SIZE 1024*1024*10
  #endif
#endif

// Prefix of g_lfact_table that is filled eagerly at startup.  Accessing those
// entries needs no synchronisation (immutable post-init); the lazy tail above
// goes through std::atomic_ref to keep the racing fill thread-safe.
#if !defined(ISOSPEC_LFACT_EAGER_FILL)
  #if ISOSPEC_BUILDING_OPENMS
    #define ISOSPEC_LFACT_EAGER_FILL ISOSPEC_G_FACT_TABLE_SIZE
  #else
    #define ISOSPEC_LFACT_EAGER_FILL 65536
  #endif
#endif

namespace IsoSpec
{

extern double* g_lfact_table;

static inline double minuslogFactorial(int n)
{
    if (n < 2)
        return 0.0;
    if (static_cast<std::size_t>(n) < static_cast<std::size_t>(ISOSPEC_LFACT_EAGER_FILL))
        return g_lfact_table[n];                              // immutable after startup

    #if ISOSPEC_BUILDING_OPENMS
    if (n >= ISOSPEC_G_FACT_TABLE_SIZE)
        return -lgamma(n+1);
    #endif

    std::atomic_ref<double> slot(g_lfact_table[n]);
    double v = slot.load(std::memory_order_relaxed);
    if (v == 0.0) {
        v = -lgamma(n+1);
        slot.store(v, std::memory_order_relaxed);
    }
    return v;
}

const double pi = 3.14159265358979323846264338328;
const double logpi = 1.144729885849400174143427351353058711647294812915311571513623071472137769884826079783623270275489708;

double NormalCDFInverse(double p);
double NormalCDFInverse(double p, double mean, double stdev);
double NormalCDF(double x, double mean, double stdev);
double NormalPDF(double x, double mean = 0.0, double stdev = 1.0);

// Returns lower incomplete gamma function of a/2, x, where a is int and > 0.
double LowerIncompleteGamma2(int a, double x);

// Returns y such that LowerIncompleteGamma2(a, y) == x. Approximately.
double InverseLowerIncompleteGamma2(int a, double x);

// Computes the inverse Cumulative Distribution Funcion of the Chi-Square distribution with k degrees of freedom
inline double InverseChiSquareCDF2(int k, double x)
{
    return InverseLowerIncompleteGamma2(k, x*tgamma(static_cast<double>(k)/2.0)) * 2.0;
}

extern thread_local std::mt19937 random_gen;

//! Uniform draw from [0,1), 53-bit resolution.
/*! Not std::uniform_real_distribution: libstdc++ routes that through
    std::generate_canonical, which takes std::log2 of a long double at every
    call -- __ieee754_logl alone was 21% of FromStochastic's runtime profile.
    Two raw 32-bit draws glued into a 53-bit mantissa produce the same
    distribution directly. Measured: 1.51x end-to-end on FromStochastic
    (C23832H37816N6528O7031S170, 3e7 molecules, Piledriver). */
inline double rdvariate_unif01(std::mt19937& rgen = random_gen)
{
    uint64_t hi = rgen();
    uint64_t lo = rgen();
    return static_cast<double>(((hi << 32) | lo) >> 11) * 0x1.0p-53;
}

inline double rdvariate_beta_1_b(double b, std::mt19937& rgen = random_gen)
{
    return 1.0 - pow(rdvariate_unif01(rgen), 1.0/b);
}


size_t rdvariate_binom(size_t tries, double succ_prob, std::mt19937& rgen = random_gen);




}  // namespace IsoSpec
