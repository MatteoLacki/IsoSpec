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

#include <cstddef>
#include "platform.h"
#include "iso_simd.h"

// ---------------------------------------------------------------------------
// Runtime ISA dispatch.
//
// A single binary can carry several builds of the two batched fill kernels --
// one per instruction-set level -- and pick the best one the CPU actually
// supports, at run time. This matters because the distributed builds have no
// -march at all (the PyPI wheel is compiled at the x86-64 baseline, so its
// vectors are 2 doubles wide), while the machines running them are mostly
// AVX2-or-better.
//
// Why it is built this way, and not the two obvious easier ways:
//
//   * __attribute__((target(...))) / target_clones() cannot do it.
//     std::experimental::simd fixes its ABI when the templates are instantiated,
//     from the -march preprocessor macros; a target attribute applied afterwards
//     changes only the *encoding*. Measured: fixed_size_simd<double,8> inside a
//     target("avx512f") function on a baseline TU emits nothing but xmm, and
//     native_simd<double> under target_clones still reports 2 lanes. So the
//     levels have to be separate translation units, each with its own -march.
//
//   * Those TUs must not include the library's own headers. A TU compiled at
//     -march=x86-64-v4 emits VEX/EVEX code into the COMDAT copies of every
//     inline and template function it instantiates, and the linker then picks
//     one copy of each for the whole program -- by link order. Measured: the
//     same two objects linked in the opposite order gave a surviving copy of a
//     shared inline function using %zmm, reachable from baseline code paths,
//     i.e. SIGILL on any CPU without AVX-512.
//
// Hence the interface below is deliberately POD-only: raw pointers and doubles,
// no IsoSpec types, nothing with a definition that could be duplicated. The
// kernel bodies (isa_kernel_impl.h) are static, and only the tables below have
// external linkage. libstdc++'s simd internals are all always_inline, so they
// leave no mergeable symbols behind either -- verified with nm on both TUs.
// ---------------------------------------------------------------------------

namespace IsoSpec
{

//! Batch-emit configurations of the current marginal-0 run into out_masses/out_probs.
/*!
    Consumes as much of the run as fills whole vectors, leaving the remaining
    fewer-than-W configurations to the caller's scalar drain. Advances
    *lProbs_ptr past what it consumed and returns how many it emitted.

    \param lProbs_ptr  In/out cursor; on entry it points at the last-emitted
                       configuration, per the generators' "advance, then emit"
                       convention.
    \param lProbs_start Index 0 of marginal 0's log-probability array.
    \param marginal_masses, marginal_probs  Marginal 0's tables, indexed as lProbs is.
    \param lcfmsv      The cutoff the run ends at.
    \param partial_mass, partial_prob  The higher marginals' contribution, constant
                       across the run.
*/
typedef std::size_t (*isa_fill_run_fn)(const double** lProbs_ptr,
                                       const double* lProbs_start,
                                       const double* marginal_masses,
                                       const double* marginal_probs,
                                       double lcfmsv,
                                       double partial_mass,
                                       double partial_prob,
                                       double* out_masses,
                                       double* out_probs);

//! As above, but scattering into a dense bin accumulator instead of storing.
/*! A configuration of mass m goes to acc[floor((m + hwmm) * inv_bin_width)],
    which the caller guarantees is one of the no_bins slots of acc. no_bins is
    only for a debug-build check, and is passed last and as a double (exact:
    far below 2^53). As a seventh integer argument it would go on the stack,
    and gcc then gave the caller a frame pointer, its spills losing the fast
    rsp-relative store forwarding of recent x86 cores (Binned()'s layered fill
    2-3% slower on x86-64-v3); anywhere but last it moves hwmm and
    inv_bin_width to other registers, after which clang interleaved the
    scatter's lanes differently (+1% on full enumerations into wide bins). */
typedef std::size_t (*isa_bin_run_fn)(const double** lProbs_ptr,
                                      const double* lProbs_start,
                                      const double* marginal_masses,
                                      const double* marginal_probs,
                                      double lcfmsv,
                                      double partial_mass,
                                      double partial_prob,
                                      double* acc,
                                      double hwmm,
                                      double inv_bin_width,
                                      double no_bins);

//! One instruction-set level's build of the kernels.
struct IsaKernels {
    isa_fill_run_fn fill_run;
    isa_bin_run_fn  bin_run;
    const char*     name;   /*!< For diagnostics and the tests. */
    std::size_t     lanes;  /*!< Doubles per vector; 1 when there is no SIMD at all. */
};

//! ISA levels, ordered by capability. Only ever grows at the top.
/*! The numeric order IS the capability order -- best_available_level() picks the
    highest usable value -- so a new level in the middle renumbers the ones above
    it. That is an ABI break for anything already compiled against this header,
    and is why ISA_LEVEL_AVX was inserted before 2.4.0 shipped rather than after.
    Any level added from now on goes on the top, or the ordering stops being
    expressible as `<`.

    Note there is deliberately no x86-64-v2 level: v2 is SSE4.2 and POPCNT, which
    do not widen a vector of doubles at all, so it would be identical to the
    baseline in every way that matters here. */
enum IsaLevel {
    ISA_LEVEL_BASELINE = 0,  /*!< Whatever the unadorned target gives: SSE2 on x86-64, NEON on AArch64. */
    ISA_LEVEL_AVX      = 1,  /*!< AVX without AVX2: 4 doubles, no FMA. Sandy/Ivy Bridge, Bulldozer-Excavator. */
    ISA_LEVEL_V3       = 2,  /*!< x86-64-v3: AVX2 + FMA, 4 doubles. */
    ISA_LEVEL_V4       = 3,  /*!< x86-64-v4: AVX-512F/BW/DQ/VL, 8 doubles. */
    ISA_LEVEL_COUNT    = 4
};

//! The kernels selected for this machine. Resolved once, on first call.
const IsaKernels& isa_kernels();

//! The level those kernels came from.
IsaLevel current_isa_level();

//! Highest level this build contains *and* this CPU supports.
IsaLevel max_supported_isa_level();

//! Is this level present in the binary and usable on this CPU?
bool isa_level_available(IsaLevel level);

//! Force a level, for testing and for working around a bad CPU.
/*! Returns false and changes nothing if the level is not available. Also
    settable through the ISOSPEC_ISA_LEVEL environment variable
    ("baseline", "v3", "v4"), read when the selection is first resolved. */
bool set_isa_level(IsaLevel level);

const char* isa_level_name(IsaLevel level);

//! What the batched kernels actually running on this machine are vectorised to.
/*!
    A stable lowercase token, never NULL and never to be freed: one of

      - "scalar"  -- no vectorisation at all. Either the standard library has no
                     <experimental/simd> (MSVC, and Apple Clang's libc++), or the
                     build disabled the SIMD path; every configuration goes
                     through the scalar drain and the kernels are do-nothing stubs.
      - "sse2"    -- x86-64's guaranteed baseline, 2 doubles per vector.
      - "avx"     -- 4 doubles, no FMA (a -march build targeting pre-Haswell x86).
      - "avx2"    -- x86-64-v3: AVX2+FMA, 4 doubles.
      - "avx512"  -- x86-64-v4: AVX-512F/BW/DQ/VL, 8 doubles.
      - "neon"    -- AArch64's guaranteed baseline, 2 doubles.
      - "simd"    -- vectorised, on an architecture this function has no name for.

    This is a finer-grained question than current_isa_level(), and deliberately
    so: that one names the dispatch bucket, which says nothing about whether the
    bucket's kernels vectorise. A build whose standard library lacks the SIMD
    header still reports level "baseline" while running entirely scalar -- the
    case this function exists to make visible.

    New tokens may be added as levels are; treat an unrecognised one as "some
    vector unit", not as an error.
*/
const char* active_simd_level();

}  // namespace IsoSpec
