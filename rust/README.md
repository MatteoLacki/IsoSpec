# isospec

Fine-structure isotopic distributions of molecules, from Rust: a safe API
over [IsoSpec](https://github.com/MatteoLacki/IsoSpec), the C++ library also
behind the `IsoSpecPy` (PyPI) and `IsoSpecR` (CRAN) packages.

Given a molecule, IsoSpec computes the masses and probabilities of its
isotopologues: all of those above a probability threshold, or the smallest
set covering a chosen share of the total probability, which it finds without
enumerating the (often astronomically many) others.

```rust
use isospec::Iso;

fn main() -> Result<(), isospec::Error> {
    let iso = Iso::from_formula("C254H377N65O75S6")?; // insulin
    let mut env = iso.envelope().total_prob(0.999)?;    // smallest set covering 99.9%
    env.sort_by_mass();
    for (mass, prob) in env.iter() {
        println!("{mass:.6}\t{prob:.6e}");
    }
    Ok(())
}
```

## Building

The C++ library is compiled from source by the build script, from sources
shipped in the crate: no system library, CMake or pre-built binary is
needed, only a C++20 compiler (tested: recent gcc and clang on Linux, Apple
clang, MSVC 2022).

| Platform | Kernels |
|---|---|
| Linux x86-64, gcc or clang with libstdc++ | SIMD; the best of SSE2, AVX, AVX2 and AVX-512 chosen at run time |
| Linux aarch64 | SIMD (NEON) |
| macOS (Apple clang) | scalar: libc++ has no `<experimental/simd>` |
| Windows (MSVC) | scalar: likewise |

"Scalar" is the same algorithms without the vectorised inner loops: correct
everywhere, somewhat slower on the largest enumerations. `isospec::simd_level()`
says what a build runs; `isospec::isa` reports the levels, and the
`ISOSPEC_ISA_LEVEL` environment variable (`baseline`, `avx`, `v3`, `v4`)
forces one, for testing.

The C++ is compiled with the optimisation level of the Cargo profile, so
debug builds of your program run IsoSpec unoptimised. To keep it fast there:

```toml
[profile.dev.package.isospec]
opt-level = 3
```

## What is there

- **Molecules** (`Iso`): from formulas, peptide sequences with
  `[UNIMOD:<id>]` modifications (an embedded Unimod table, or your own CSV),
  plain FASTA sequences, element `Composition`s, or raw isotope masses and
  abundances for labelled or otherwise non-natural elements. Monoisotopic,
  average, mode, lightest and heaviest masses, variance, configurations.
- **Envelopes** (`Envelope`): threshold, total-probability, simulated
  (stochastic, seedable) and binned distributions, or from your own arrays;
  sorting, scaling, normalising, mass shifts, resampling, binning, sums,
  convolutions, linear combinations, statistics, and Wasserstein distances
  (plain, oriented, abyssal, matching).
- **Generators**: `ThresholdGenerator`, `LayeredGenerator`,
  `OrderedGenerator` (strictly decreasing probability) and
  `StochasticGenerator` stream isotopologues without storing them, as
  iterators (`for peak in generator`) or step by step with configurations
  (`generator.advance()`, `generator.conf()`).
- **Elements and compositions**: the built-in element table with isotope
  masses and abundances; compositions with arithmetic, parsing and printing.

Seeded sampling gives the same result as IsoSpecPy and the C++ library for
the same seed, on every platform.

See `rust/examples/` in the repository for runnable programs.

## Versions

The crate follows semver on its own version numbers; the C++ library it
bundles is versioned separately, as `isospec::LIBRARY_VERSION`, and each
crate release names the library version it ships (see `CHANGELOG.md`).

The minimum supported Rust version is 1.88 (that of `cxx`). It is raised
only in a release that bumps the crate's minor version (before 1.0, a
release that may also change the API), never in a patch release.

## Threads

All types can move between threads. `Iso`, `Envelope` and `UnimodTable` are
`Sync`: their `&self` methods may run concurrently. The library itself is
single-threaded per call; run independent molecules on separate threads
(rayon over peptides, for instance).

## Errors

Invalid input is an `Error`, whose `ErrorKind` names what went wrong; it is
never undefined behaviour. The C++ library's exceptions are caught at the
boundary, the arguments it cannot check itself are checked before the call,
and the API is fuzzed (with the C++ under AddressSanitizer).

## Licence and citation

BSD 2-clause. If IsoSpec is useful in your research, please cite:

- Łącki, M. K.; Valkenborg, D.; Startek, M. P. *IsoSpec2: Ultrafast Fine
  Structure Calculator.* Analytical Chemistry **2020**, *92* (14), 9472–9475.
  <https://doi.org/10.1021/acs.analchem.0c00959>
- Łącki, M. K.; Startek, M.; Valkenborg, D.; Gambin, A. *IsoSpec: Hyperfast
  Fine Structure Calculator.* Analytical Chemistry **2017**, *89* (6),
  3272–3277. <https://doi.org/10.1021/acs.analchem.6b01459>
