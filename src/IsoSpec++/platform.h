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

#include "platform_incl.h"

#if defined(__unix__) || defined(__unix) || \
        (defined(__APPLE__) && defined(__MACH__))
#define ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY true
#define ISOSPEC_TEST_WE_ARE_ON_WINDOWS false /* CYGWIN doesn't really count as Windows for our purposes, we'll be using UNIX API anyway */
#define ISOSPEC_TEST_GOT_SYSTEM_MMAN true
#define ISOSPEC_TEST_GOT_MMAN true
#elif defined(__MINGW32__) || defined(_WIN32)
#define ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY false
#define ISOSPEC_TEST_WE_ARE_ON_WINDOWS true
#define ISOSPEC_TEST_GOT_SYSTEM_MMAN false
#define ISOSPEC_TEST_GOT_MMAN true
#else
#define ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY false /* Well, probably... */
#define ISOSPEC_TEST_WE_ARE_ON_WINDOWS false
#define ISOSPEC_TEST_GOT_SYSTEM_MMAN false
#define ISOSPEC_TEST_GOT_MMAN false
#endif

#if !defined(ISOSPEC_WE_ARE_ON_UNIX_YAY)
#define ISOSPEC_WE_ARE_ON_UNIX_YAY ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY
#endif

#if !defined(ISOSPEC_WE_ARE_ON_WINDOWS)
#define ISOSPEC_WE_ARE_ON_WINDOWS ISOSPEC_TEST_WE_ARE_ON_WINDOWS
#endif

#if !defined(ISOSPEC_GOT_SYSTEM_MMAN)
#define ISOSPEC_GOT_SYSTEM_MMAN ISOSPEC_TEST_GOT_SYSTEM_MMAN
#endif

#if !defined(ISOSPEC_GOT_MMAN)
#define ISOSPEC_GOT_MMAN ISOSPEC_TEST_GOT_MMAN
#endif


// Note: __GNUC__ is defined by clang and gcc
#ifdef __GNUC__
#define ISOSPEC_IMPOSSIBLE(condition) if(condition) __builtin_unreachable();
#define ISOSPEC_LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define ISOSPEC_UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)
// For aggressive inlining
#define ISOSPEC_FORCE_INLINE __attribute__ ((always_inline)) inline
#elif defined _MSC_VER
#define ISOSPEC_IMPOSSIBLE(condition) __assume(!(condition));
#define ISOSPEC_LIKELY(condition) condition
#define ISOSPEC_UNLIKELY(condition) condition
#define ISOSPEC_FORCE_INLINE __forceinline
#else
#define ISOSPEC_IMPOSSIBLE(condition)
#define ISOSPEC_LIKELY(condition) condition
#define ISOSPEC_UNLIKELY(condition) condition
#define ISOSPEC_FORCE_INLINE inline
#endif

#ifdef ISOSPEC_DEBUG
#undef ISOSPEC_IMPOSSIBLE
#include <cassert>
#define ISOSPEC_IMPOSSIBLE(condition) assert(!(condition));
#endif /* ISOSPEC_DEBUG */


#if ISOSPEC_GOT_MMAN
    #if ISOSPEC_GOT_SYSTEM_MMAN
        #include <sys/mman.h>
    #else
        #include "mman.h"
    #endif
#else
    #include <stdlib.h>     /* malloc, free, rand */
#endif


#if defined(OPENMS_DLLAPI) /* IsoSpec is being built as a part of OpenMS: use their visibility macros */
#define ISOSPEC_EXPORT_SYMBOL OPENMS_DLLAPI
#else /* it's a can of worms we don't yet want to open ourselves though... */
#define ISOSPEC_EXPORT_SYMBOL
#endif

#if (defined(_WIN32) || defined(_WIN64)) && !defined(MXE)
#define ISOSPEC_C_API __declspec(dllexport)
#else
#define ISOSPEC_C_API
#endif

#if !defined(__cpp_if_constexpr)
#define constexpr_if if
#define ISOSPEC_MAYBE_UNUSED
#else
#define constexpr_if if constexpr
#define ISOSPEC_MAYBE_UNUSED [[maybe_unused]]
#endif


// C++-only from here on: cwrapper.h includes this header, and that one must
// stay includable from a C translation unit (it is the C ABI's header).
#ifdef __cplusplus

#include <cstddef>  // std::size_t, used below whether or not SIMD is available

#ifdef __has_include
    // Deliberately NOT preferring <simd>: the standardized std::simd (P1928,
    // landing piecemeal in libstdc++/libc++) does not expose the experimental
    // TS's native_simd alias, so __has_include(<simd>) can be true while the
    // header provides none of the API this codebase uses. GCC 16's <simd> is
    // exactly this case (compiles under -std=c++2c, not under our -std=c++20)
    // and silently broke the macOS build when given priority here. Revisit
    // once std::simd's actual shipped API is verified against our usage.
    #if __has_include(<experimental/simd>)
        #define ISOSPEC_HAS_SIMD 1
        #include <experimental/simd>
        namespace simd_ns = std::experimental;
        using simd_double = simd_ns::native_simd<double>;
        constexpr std::size_t DOUBLE_SIMD_ALIGNMENT =
                alignof(simd_ns::native_simd<double>);
    #endif
#endif

#if !defined(ISOSPEC_HAS_SIMD)
    #define ISOSPEC_HAS_SIMD 0

    // Hopefully <simd> will be available everywhere before someone
    // comes up with a CPU that has 256-byte SIMD registers.
    constexpr std::size_t DOUBLE_SIMD_ALIGNMENT = 128;
#endif

#endif  /* __cplusplus */