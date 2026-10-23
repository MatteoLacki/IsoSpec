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

#include <atomic>
#include <cstdlib>
#include <cstring>
#include "isa_kernels.h"

// Compiled at the baseline, always. It must not be built for a raised -march:
// it runs before anything has established what the CPU can do.

#if !ISOSPEC_ISA_DISPATCH
// No dispatch: the kernels are compiled straight into this translation unit at
// whatever -march the build used, and the "selection" is a constant. This is
// what a plain unity build, and the R package, get.
#define ISOSPEC_ISA_NAMESPACE isa_native
#define ISOSPEC_ISA_NAME "native"
#include "isa_kernel_impl.h"
#endif

namespace IsoSpec
{

#if ISOSPEC_ISA_DISPATCH

// Defined by the isa_kernels_*.cpp translation units, each compiled with its own
// -march. Only ever referenced through these declarations: nothing from those
// TUs is inlined here, which is the whole point.
namespace isa_baseline { extern const IsaKernels kernels; }
#if ISOSPEC_ISA_BUILD_X86
namespace isa_avx { extern const IsaKernels kernels; }
namespace isa_v3 { extern const IsaKernels kernels; }
namespace isa_v4 { extern const IsaKernels kernels; }
#endif

namespace
{

constexpr const IsaKernels* level_table(IsaLevel level)
{
    switch(level)
    {
        case ISA_LEVEL_BASELINE: return &isa_baseline::kernels;
#if ISOSPEC_ISA_BUILD_X86
        case ISA_LEVEL_AVX:      return &isa_avx::kernels;
        case ISA_LEVEL_V3:       return &isa_v3::kernels;
        case ISA_LEVEL_V4:       return &isa_v4::kernels;
#else
        case ISA_LEVEL_AVX:      return nullptr;
        case ISA_LEVEL_V3:       return nullptr;
        case ISA_LEVEL_V4:       return nullptr;
#endif
        case ISA_LEVEL_COUNT:
        default:                 return nullptr;
    }
}

// Whether the CPU (and the OS, which must have enabled the register state) can
// run a level. __builtin_cpu_supports covers both: libgcc's initialiser checks
// the CPUID feature bits and XCR0, so an AVX-capable CPU under a kernel that
// does not save the wide registers correctly reports false.
bool cpu_supports(IsaLevel level)
{
    if(level == ISA_LEVEL_BASELINE)
        return true;  // architecturally guaranteed, by definition of "baseline"

#if ISOSPEC_ISA_BUILD_X86 && (defined(__GNUC__) || defined(__clang__))
    switch(level)
    {
        case ISA_LEVEL_AVX:
            return __builtin_cpu_supports("avx");
        case ISA_LEVEL_V3:
            return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
        case ISA_LEVEL_V4:
            // x86-64-v4 is F+BW+DQ+VL together; a CPU with only some of those
            // cannot run code compiled for the set.
            return __builtin_cpu_supports("avx512f")
                && __builtin_cpu_supports("avx512bw")
                && __builtin_cpu_supports("avx512dq")
                && __builtin_cpu_supports("avx512vl");
        default:
            return false;
    }
#else
    // No way to ask: stay on the level that is always safe.
    return false;
#endif
}

// A level is usable only if it is both present in this binary and runnable here.
// "Present" is not a given: on a non-x86 build, and on a toolchain too old for
// the AVX-512 flags, the corresponding TU compiles to nothing and its table has
// no kernels.
bool level_usable(IsaLevel level)
{
    const IsaKernels* t = level_table(level);
    return t != nullptr && t->fill_run != nullptr && cpu_supports(level);
}

IsaLevel level_from_name(const char* s)
{
    if(strcmp(s, "baseline") == 0) return ISA_LEVEL_BASELINE;
    if(strcmp(s, "avx") == 0)      return ISA_LEVEL_AVX;
    if(strcmp(s, "v3") == 0)       return ISA_LEVEL_V3;
    if(strcmp(s, "v4") == 0)       return ISA_LEVEL_V4;
    return ISA_LEVEL_COUNT;        // unrecognised
}

IsaLevel best_available_level()
{
    IsaLevel best = ISA_LEVEL_BASELINE;
    for(int l = ISA_LEVEL_BASELINE; l < ISA_LEVEL_COUNT; l++)
        if(level_usable(static_cast<IsaLevel>(l)))
            best = static_cast<IsaLevel>(l);
    return best;
}

// The selection. Resolved on first use rather than during static
// initialisation, so it cannot lose a race with the initialisation of
// anything it depends on -- but deliberately *not* via a function-local
// static: that would need a runtime guard (a lock, on most ABIs) to stay
// safe if two threads race to initialize it for the first time, which is
// unnecessary here. Computing the answer is a pure, deterministic function
// of CPU features and an environment variable, so letting two racing
// threads compute it twice is harmless; only the publication needs care.
struct Selection {
    IsaLevel level;
    const IsaKernels* kernels;
};

// One immutable record per level, so that publishing a selection is a single
// pointer store. With the level and the kernels in two separate atomics, a
// reader racing set_isa_level() could pair one call's level with another's
// kernels. constinit: the addresses of the extern tables are link-time
// constants, so this is never part of dynamic initialisation -- a reader
// during another TU's static initialisation still sees it filled in.
constinit const Selection level_records[ISA_LEVEL_COUNT] = {
    {ISA_LEVEL_BASELINE, level_table(ISA_LEVEL_BASELINE)},
    {ISA_LEVEL_AVX,      level_table(ISA_LEVEL_AVX)},
    {ISA_LEVEL_V3,       level_table(ISA_LEVEL_V3)},
    {ISA_LEVEL_V4,       level_table(ISA_LEVEL_V4)},
};

const Selection* compute_selection()
{
    IsaLevel level = best_available_level();

    // An override lets a user work around a CPU whose wide units are a trap
    // (or slower), and lets the test suite exercise every level the host can
    // actually run. Ignored, rather than rejected, if it names something
    // unavailable -- a stray environment variable must not stop the library.
    //
    // getenv() is safe against concurrent getenv(); a host program calling
    // setenv() on another thread at the same moment is POSIX's problem, not
    // one this library can solve.
    const char* env = getenv("ISOSPEC_ISA_LEVEL");
    if(env != nullptr)
    {
        const IsaLevel wanted = level_from_name(env);
        if(wanted != ISA_LEVEL_COUNT && level_usable(wanted))
            level = wanted;
    }

    return &level_records[level];
}

std::atomic<const Selection*> g_selection{nullptr};

const Selection& selection()
{
    const Selection* s = g_selection.load(std::memory_order_acquire);
    if(s != nullptr)
        return *s;

    // Publish only if nothing has been published meanwhile: an explicit
    // set_isa_level() that landed while this thread was computing the default
    // must win over it, not be silently overwritten. On failure `expected`
    // holds whatever got there first.
    const Selection* expected = nullptr;
    const Selection* computed = compute_selection();
    if(g_selection.compare_exchange_strong(expected, computed, std::memory_order_acq_rel, std::memory_order_acquire))
        return *computed;
    return *expected;
}

}  // anonymous namespace

const IsaKernels& isa_kernels()
{
    return *selection().kernels;
}

IsaLevel current_isa_level()
{
    return selection().level;
}

IsaLevel max_supported_isa_level()
{
    return best_available_level();
}

bool isa_level_available(IsaLevel level)
{
    return level_usable(level);
}

bool set_isa_level(IsaLevel level)
{
    if(!level_usable(level))
        return false;
    g_selection.store(&level_records[level], std::memory_order_release);
    return true;
}

#else  // !ISOSPEC_ISA_DISPATCH

const IsaKernels& isa_kernels()
{
    return isa_native::kernels;
}

IsaLevel current_isa_level() { return ISA_LEVEL_BASELINE; }

IsaLevel max_supported_isa_level() { return ISA_LEVEL_BASELINE; }

bool isa_level_available(IsaLevel level) { return level == ISA_LEVEL_BASELINE; }

bool set_isa_level(IsaLevel level) { return level == ISA_LEVEL_BASELINE; }

#endif  // ISOSPEC_ISA_DISPATCH

namespace
{

// What this translation unit's own -march gives. It is compiled with the
// build's ambient flags -- exactly the ones isa_kernels_baseline.cpp gets, since
// IsaDispatch.cmake raises -march only for the raised levels' sources -- so this is
// the honest answer for the baseline level too, including the `-march=native`
// builds where "baseline" is in fact as wide as the host allows.
const char* ambient_simd_name()
{
#if defined(__AVX512F__)
    return "avx512";
#elif defined(__AVX2__)
    return "avx2";
#elif defined(__AVX__)
    return "avx";
#elif ISOSPEC_ISA_BUILD_X86
    return "sse2";
#elif defined(__aarch64__) || defined(__ARM_NEON) || defined(_M_ARM64)
    return "neon";
#else
    return "simd";
#endif
}

}  // anonymous namespace

const char* active_simd_level()
{
    const IsaKernels& k = isa_kernels();

    // The stub tables: a level that was not built (fill_run == nullptr), or one
    // built without any SIMD to build with, which publishes a kernel that emits
    // nothing and leaves the whole run to the caller's scalar drain.
    if(k.fill_run == nullptr || k.lanes <= 1)
        return "scalar";

#if ISOSPEC_ISA_DISPATCH
    // The raised levels are known exactly: their sources are compiled with the
    // flags naming them, and refuse to publish a table otherwise.
    switch(current_isa_level())
    {
        case ISA_LEVEL_AVX: return "avx";
        case ISA_LEVEL_V3: return "avx2";
        case ISA_LEVEL_V4: return "avx512";
        default:           break;  // baseline: whatever the ambient flags gave
    }
#endif

    return ambient_simd_name();
}

const char* isa_level_name(IsaLevel level)
{
    switch(level)
    {
        case ISA_LEVEL_BASELINE: return "baseline";
        case ISA_LEVEL_AVX:      return "avx";
        case ISA_LEVEL_V3:       return "v3";
        case ISA_LEVEL_V4:       return "v4";
        case ISA_LEVEL_COUNT:
        default:                 return "unknown";
    }
}

}  // namespace IsoSpec
