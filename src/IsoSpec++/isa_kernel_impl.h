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

// The bodies of the two batched fill kernels, in a form meant to be included
// once per instruction-set level, by a translation unit compiled with that
// level's -march. Define ISOSPEC_ISA_NAMESPACE to the namespace to put them in
// before including this; everything lands there, so several levels can coexist
// in one binary.
//
// Nothing here may be given external linkage except the kernels table at the
// bottom, and nothing outside <experimental/simd> and this project's POD-only
// headers may be included -- see the long note in isa_kernels.h for what goes
// wrong otherwise. This file therefore has NO include guard: it is included
// several times on purpose, with a different namespace each time.

#ifndef ISOSPEC_ISA_NAMESPACE
#error "define ISOSPEC_ISA_NAMESPACE before including isa_kernel_impl.h"
#endif

#include <cstddef>
#include "platform.h"      // ISOSPEC_HAS_SIMD, simd_double, simd_ns -- macros and typedefs only
#include "isa_kernels.h"   // the POD interface these implement

namespace IsoSpec
{
namespace ISOSPEC_ISA_NAMESPACE
{

#if ISOSPEC_HAS_SIMD

// Both kernels below walk the run the same way. Each step takes the W positions
// after the cursor at once: one compare decides that all W are above the cutoff
// -- sound because lProbs descends, so the highest index guards the rest, and
// safe to read because of the ISOSPEC_LPROB_GUARDIANS sentinels past the end --
// and then the marginal's mass and probability tables are read as vectors and
// combined with the higher marginals' (loop-invariant) partial mass and product.
//
// Because the compare only succeeds when the *last* lane is above the cutoff,
// every lane produced is a real configuration, so the reads from the mass and
// probability tables, which carry no sentinels, stay in bounds.
//
// The output is packed contiguously and the caller's scalar drain moves the
// cursor by a non-multiple of W, so neither the loads nor the stores can be
// assumed vector-aligned: element_aligned throughout.

static std::size_t fill_run(const double** lProbs_ptr,
                            const double* lProbs_start,
                            const double* marginal_masses,
                            const double* marginal_probs,
                            double lcfmsv,
                            double partial_mass,
                            double partial_prob,
                            double* out_masses,
                            double* out_probs)
{
    constexpr std::size_t W = simd_double::size();
    const double* p = *lProbs_ptr;
    std::size_t emitted = 0;

    const simd_double vpartial_prob(partial_prob);
    const simd_double vpartial_mass(partial_mass);

    simd_double masses;
    simd_double probs;

    while(ISOSPEC_LIKELY(*(p + W) >= lcfmsv))
    {
        const std::size_t offset = static_cast<std::size_t>(p - lProbs_start) + 1;

        probs.copy_from(marginal_probs + offset, simd_ns::element_aligned);
        probs *= vpartial_prob;
        masses.copy_from(marginal_masses + offset, simd_ns::element_aligned);
        masses += vpartial_mass;

        masses.copy_to(out_masses + emitted, simd_ns::element_aligned);
        probs.copy_to(out_probs + emitted, simd_ns::element_aligned);

        emitted += W;
        p += W;
    }

    *lProbs_ptr = p;
    return emitted;
}

static std::size_t bin_run(const double** lProbs_ptr,
                           const double* lProbs_start,
                           const double* marginal_masses,
                           const double* marginal_probs,
                           double lcfmsv,
                           double partial_mass,
                           double partial_prob,
                           double* acc,
                           double hwmm,
                           double inv_bin_width,
                           double no_bins)
{
    constexpr std::size_t W = simd_double::size();
    const double* p = *lProbs_ptr;
    std::size_t emitted = 0;

    const simd_double vpartial_prob(partial_prob);
    const simd_double vpartial_mass(partial_mass);
    const simd_double vhwmm(hwmm);
    const simd_double vinv_bin_width(inv_bin_width);
    (void)no_bins;  // only for the debug check below

    simd_double masses;
    simd_double probs;

    while(ISOSPEC_LIKELY(*(p + W) >= lcfmsv))
    {
        const std::size_t offset = static_cast<std::size_t>(p - lProbs_start) + 1;

        probs.copy_from(marginal_probs + offset, simd_ns::element_aligned);
        probs *= vpartial_prob;
        masses.copy_from(marginal_masses + offset, simd_ns::element_aligned);
        masses += vpartial_mass;

        // The vector form of bin_slot(): the same arithmetic in the same order,
        // so it lands in the same bin and, like it, stays inside the allocation
        // (no floor: the values are >= 0, so the conversion's truncation is it).
        // The scatter that consumes it cannot be vectorised -- the indices repeat
        // and collide, and resolving that would need AVX-512CD conflict detection
        // -- but everything ahead of it is, which is where the work is.
        const simd_double slots = (masses + vhwmm) * vinv_bin_width;

        // Hand the probabilities to the scatter through a buffer rather than by
        // subscripting the vector. Clang otherwise sees every lane extracted and
        // un-vectorises the multiply above into W scalar vmulsd's, each with its
        // own 8-byte load, putting a load->mul->add->store chain on every lane:
        // measured 2.03x slower on Alder Lake, 1.36x on Piledriver. The spill
        // itself is free -- these values are consumed one lane at a time either
        // way, so the round trip through memory was always going to happen; this
        // only fixes where the multiply lands relative to it. GCC is unaffected
        // (it kept the packed form regardless) and the results are bit-identical.
        alignas(DOUBLE_SIMD_ALIGNMENT) double pbuf[W];
        probs.copy_to(pbuf, simd_ns::vector_aligned);

        // The check is compiled out of release builds altogether: even as a
        // named slot variable around the same arithmetic, it changed how clang
        // interleaves these lanes.
        for(std::size_t ii = 0; ii < W; ii++)
        {
#ifdef ISOSPEC_DEBUG
            ISOSPEC_IMPOSSIBLE(!(slots[ii] >= 0.0 && slots[ii] < no_bins));
#endif
            acc[static_cast<std::ptrdiff_t>(slots[ii])] += pbuf[ii];
        }

        emitted += W;
        p += W;
    }

    *lProbs_ptr = p;
    return emitted;
}

constexpr std::size_t kernel_lanes = simd_double::size();

#else  // no <experimental/simd> on this toolchain

// Emit nothing and let the caller's scalar drain do all the work: that drain is
// a complete implementation on its own, so this needs no scalar duplicate of
// the loops above.

static std::size_t fill_run(const double**, const double*, const double*, const double*,
                            double, double, double, double*, double*)
{
    return 0;
}

static std::size_t bin_run(const double**, const double*, const double*, const double*,
                           double, double, double, double*, double, double, double)
{
    return 0;
}

constexpr std::size_t kernel_lanes = 1;

#endif  // ISOSPEC_HAS_SIMD

// The one symbol with external linkage. `extern` is load-bearing: a const object
// at namespace scope would otherwise have internal linkage and the whole level
// would be optimised out of the object file.
extern const IsaKernels kernels;
extern const IsaKernels kernels = { &fill_run, &bin_run, ISOSPEC_ISA_NAME, kernel_lanes };

}  // namespace ISOSPEC_ISA_NAMESPACE
}  // namespace IsoSpec
