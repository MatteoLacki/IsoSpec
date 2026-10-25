//! The instruction-set level the library's SIMD kernels run at.
//!
//! The kernels are compiled for several x86-64 levels and the best one the
//! CPU supports is chosen at run time (it can also be forced with the
//! `ISOSPEC_ISA_LEVEL` environment variable: `baseline`, `avx`, `v3`, `v4`).
//! The choice is process-wide. Off x86-64, and with compilers lacking
//! `<experimental/simd>` (MSVC, Apple clang), only the baseline exists.

use crate::ffi::bridge as sys;

/// An instruction-set level of the kernels, in increasing capability.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum IsaLevel {
    /// What the target gives without extra flags: SSE2 on x86-64, NEON on aarch64.
    Baseline,
    /// AVX without AVX2: 4 doubles, no FMA.
    Avx,
    /// x86-64-v3: AVX2 and FMA.
    V3,
    /// x86-64-v4: AVX-512 (F, BW, DQ, VL).
    V4,
}

impl IsaLevel {
    /// All levels, lowest first.
    pub const ALL: [IsaLevel; 4] = [
        IsaLevel::Baseline,
        IsaLevel::Avx,
        IsaLevel::V3,
        IsaLevel::V4,
    ];

    fn from_raw(level: i32) -> IsaLevel {
        IsaLevel::ALL[level.clamp(0, 3) as usize]
    }

    fn raw(self) -> i32 {
        self as i32
    }

    /// The level's name: `"baseline"`, `"avx"`, `"v3"` or `"v4"`.
    pub fn name(self) -> &'static str {
        sys::isa_level_name(self.raw())
    }

    /// Whether this build has the level and the CPU can run it.
    pub fn is_available(self) -> bool {
        sys::isa_level_available(self.raw())
    }
}

/// The level in use.
pub fn current() -> IsaLevel {
    IsaLevel::from_raw(sys::isa_current_level())
}

/// The highest level available.
pub fn max_supported() -> IsaLevel {
    IsaLevel::from_raw(sys::isa_max_supported_level())
}

/// Switches the whole process to `level`; returns `false`, changing nothing,
/// if it is not available. For testing and benchmarking: the default is
/// already the best available.
pub fn set(level: IsaLevel) -> bool {
    sys::isa_set_level(level.raw())
}

/// What the kernels actually run as: `"scalar"` (no SIMD at all), `"sse2"`,
/// `"neon"`, `"avx"`, `"avx2"`, `"avx512"`, ...
pub fn simd_level() -> &'static str {
    sys::simd_level()
}
