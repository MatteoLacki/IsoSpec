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

// The library never calls the C library's lgamma(). It reports the sign of
// Gamma(x) through the global `signgam`, so every call writes shared state --
// a data race between threads -- and lgamma_r(), which avoids that, is missing
// from MSVC's and MinGW's runtimes (MinGW's own lgamma, linked even under UCRT,
// writes signgam too). Every argument the library needs is an integer or a
// half-integer, so lnGamma is served instead by: exact values below 8, compiled
// in; Stirling's series from 8 up; and, for integers, the log-factorial cache
// below. Results are within 2 ulp of correctly rounded (glibc >= 2.43 rounds
// correctly), at about half the cost: 4.75 vs 9.06 ns per value on x86-64-v3.

// lnGamma(x) for x >= 8, by Stirling's series. Calls only log(), which writes
// no shared state for such arguments.
double lgamma_stirling(double x) noexcept;

extern double* g_lfact_table;

// -ln(n!) = -lnGamma(n+1), from a cache: entries below 7 are compiled in, the
// rest up to ISOSPEC_LFACT_EAGER_FILL are filled at startup, the remainder of
// the table on first use; beyond the table it is computed every time.
static inline double minuslogFactorial(int n)
{
    if (n < 2)
        return 0.0;
    if (static_cast<std::size_t>(n) < static_cast<std::size_t>(ISOSPEC_LFACT_EAGER_FILL))
        return g_lfact_table[n];                              // immutable after startup

    if (n >= ISOSPEC_G_FACT_TABLE_SIZE)
        return -lgamma_stirling(n + 1.0);

    std::atomic_ref<double> slot(g_lfact_table[n]);
    double v = slot.load(std::memory_order_relaxed);
    if (v == 0.0) {
        v = -lgamma_stirling(n + 1.0);
        slot.store(v, std::memory_order_relaxed);
    }
    return v;
}

// lnGamma(two_x / 2) for two_x >= 2, i.e. at an integer or half-integer. The
// argument is passed doubled so that no other argument is expressible. Integers
// come from the log-factorial cache; half-integers below 8 from a small table,
// above it from Stirling. Not on any hot path (setup-time size estimates), so
// half-integers get no cache of their own.
double lgamma_half(std::int64_t two_x) noexcept;

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

//! A generator for reproducible sampling: the same seed gives the same stream
//! on every platform and through every binding (C ABI, Python, Rust), since
//! std::mt19937 and std::seed_seq are fully specified by the standard. Both
//! halves of the seed are used.
inline std::mt19937 seeded_rng(std::uint64_t seed)
{
    std::seed_seq seq{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)};
    return std::mt19937(seq);
}

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

//! Pivot positions for the quickselects (trimming FromTotalProb's envelope,
//! quickselect() in misc.cpp): a splitmix64 stream from a fixed start, so a
//! selection, and the order it leaves the array in, is the same on every run.
//! Drawing from random_gen instead (seeded from std::random_device) made
//! FromTotalProb's output order, and which of several equally probable
//! configurations survive the cut, differ between runs of the same program.
//! Pseudo-random rather than e.g. middle-of-range so the expected cost stays
//! linear whatever order the generators emit configurations in.
class PivotSequence
{
    std::uint64_t state = 0;

 public:
    //! A position in [0, len), len > 0.
    std::size_t next(std::size_t len) noexcept
    {
        std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z ^= z >> 31;
        return static_cast<std::size_t>(z % len);
    }
};




}  // namespace IsoSpec
