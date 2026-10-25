# The Rust crate (`isospec`)

A safe Rust API over the C++ library, published as one crate, `isospec`
(0.x, edition 2024, MSRV 1.88, the floor of cxx 1.0.205). It compiles the C++
from source in its build script: no system library, no CMake.

## Layout

`Cargo.toml` sits at the repository root and everything Rust under `rust/`.
The root position is the point: `cargo package` only ships files inside the
crate directory, and from the root the `include` list can name
`src/IsoSpec++/*.h` and `*.cpp` directly, so the published crate builds from
exactly the sources in the repository, with no copy to keep in sync (the
alternative, a `rust/` crate with a script copying the C++ in before every
package, was rejected for that reason). Consequences:

- `[lib] path`, `build`, and explicit `[[test]]`/`[[example]]` entries with
  `autotests = false` etc.: the repository's `tests/` (C++, Python, R) and
  `src/` are not cargo's.
- An empty `[workspace]`, so the crate is never taken for a member of an
  enclosing workspace, and a crate further down (`rust/fuzz`, which has its
  own `[workspace]`) is not made a member of this one.
- `.gitignore` has `/target/` and `/Cargo.lock`.

## Binding: cxx to the C++ classes, not the C ABI

Python goes through `cwrapper.h`; Rust binds the C++ classes through cxx,
which catches exceptions at the boundary (a C ABI would have to encode every
one as NULL/NaN) and lets the safe types hold `UniquePtr`s to the real
objects. What cxx cannot bind directly -- constructors, overloads, default
arguments, class templates (the generators are templates behind `using`
aliases), arrays of pointers (`LinearCombination`) -- goes through
`rust/shim/shim.{h,cc}`: thin free functions and wrapper structs, nothing
that computes.

- **Exceptions.** cxx's default translation keeps only `what()`. The shim
  defines `rust::behavior::trycatch`, catching the standard classes most
  derived first and prefixing the message with a letter (`I` invalid_argument
  or domain_error, `R` out_of_range, `L` length_error, `O` overflow_error, `M`
  bad_alloc, `X` anything else); `error.rs` turns it back into `ErrorKind`.
  Every bridge function that can throw is declared `Result`: one that is not,
  and throws, calls `std::terminate`.
- **Checks the C++ cannot make.** It takes raw pointers plus counts, so the
  Rust side checks that slices are as long as the counts say (isotope masses
  vs probabilities, `from_arrays`, configuration buffers) and that element
  indices are in range before calling. Values are validated by the C++ (the
  hardened constructors and envelope operations), so their errors come back
  as the C++ wrote them.
- **Null slices.** An empty envelope's arrays may be null pointers, from which
  a `rust::Slice` (and a Rust slice) may not be built: the shim returns an
  empty slice instead.
- **Elements** are a row index into the C++ element table (`elem_table_*`),
  only ever obtained from the library, so always valid.
- **Unimod tables**: the embedded one and those loaded by path live for the
  whole process in the C++ (the path cache is never evicted), so they are
  `&'static`; CSV-parsed ones are owned.
- **The stochastic generator** keeps a reference to the `std::mt19937` it is
  given, so the shim's wrapper owns a seeded one next to it (declared first:
  constructed before, destroyed after). Unseeded, it is given the library's
  `random_gen`, the thread-local generator of whichever thread is sampling.

## Keeping 0.x changes non-breaking

- `links = "isospec"`: two semver-incompatible versions in one build would
  both link the C++ symbols from static archives, and the linker silently
  keeps one definition of each (one version's Rust calling the other's C++).
  Cargo refuses such a graph instead. It had to be there from the first
  release: a published version without it can always coexist with later ones.
- Option structs (`ThresholdOptions`, `LayeredOptions`, `StochasticOptions`,
  `PeptideOptions`) and `ElementSpec` are `#[non_exhaustive]`: built from
  `default()`/`new` plus `with_*` setters (public fields stay readable and
  assignable), so a new C++ parameter is a new field, not a breaking change.
  Output structs (`Peak`, `Isotope`, `UnimodEntry`, `Match`,
  `ExpandedComposition`) likewise, and the enums that can grow (`IsaLevel`,
  `Masses`, `ErrorKind`); `IsaLevel::ALL` is a slice for the same reason.
  `Threshold` stays exhaustive: absolute or relative is all a threshold can be.
- `isa::set` is `#[doc(hidden)]`, for the crate's tests and benchmarks; not
  part of the API (`ISOSPEC_ISA_LEVEL` is the supported way to force a level).

## Ownership and threads

- Generators take their `Iso` by value (the C++ constructors move out of it).
- Envelopes are built from `&Iso` (the C++ copies it).
- `Iso`, `Envelope`, `UnimodTable`: `Send + Sync`. AGENTS.md's Threading
  section is the basis: const members may run concurrently. The safe API puts
  every non-const C++ call (sorting, `get_total_prob`, the statistics, the
  Wasserstein distances, which sort both arguments) behind `&mut`.
- Generators: `Send`, not `Sync`.
- `tests/threads.rs` and the TSan run below are the evidence.

## Generators: cursor and batched iterator

`advance()` steps one isotopologue per C++ call and is the only path that
can report the configuration (`conf()`, valid until the next step: enforced
by the borrow). Iteration (`IntoIterator`, `Peaks`) is separate and batched:
one call fills 256 masses, probabilities and log-probabilities.

Measured (best of 7, threshold enumerations of 28k-38k peaks, release-wheel
flags): one call per peak cost 15-20 ns/peak, against 2.2 for the same C++
loop. Most of it was store-forwarding stalls: a peak written as 8-byte
stores (by the C++, or into `self.current`) and immediately read back as
16-byte loads. The batch removes the per-peak handover. (Separately, the C++
asserts were live in release builds, because cc never defines `NDEBUG`; it is
now defined when Cargo's debug assertions are off, which brought one
envelope case from 9.4 to 6.0 ns/peak, level with C++.) `__restrict` on the fill loop's outputs
(without it every store may alias the generator's state, which is then
reloaded) and an inlined `next` with a cold refill bring it to:

| host | Rust iterator | C++ loop, all three fields used | envelope Rust / C++ |
|---|---|---|---|
| nighthaven (x86-64-v3) | 2.0-3.3 ns/peak | 0.3-0.7 | 3.8-4.7 / 4.5-5.4 (parity, noise both ways) |
| solace (x86-64-v4, WSL) | 2.8-3.2 | 0.8-1.2 | 6.0 / 6.1 |
| spot (M1, macOS, scalar kernels) | 2.9 | 1.1 | 2.9 / 3.0 |

The remaining ~2 ns is writing the peaks to memory and reading them back,
which a hand-written C++ loop (vectorised by the compiler) never does.
Envelopes, which do their work entirely in C++, are at parity; the ordered
generator, dominated by its priority queue, likewise.

## Not exposed, deliberately

- The raw-pointer constructors and `release_*` of `FixedEnvelope`
  (`to_vecs()` copies instead).
- `tabSize`/`hashSize` (64 is the measured best; `hashSize` does nothing).
- The `Marginal` classes and generator internals: custom elements need only
  `Iso::from_isotopes`/`add_element`.
- `quickselect` and the rest of `misc.h`; Python-only helpers (plotting).

## Testing

- `rust/tests/`: every constructor, accessor and operation, with expected
  values from IsoSpecPy 2.5.0; seeded sampling gives Python's exact counts.
  Every `Result` function has a failing case.
- `rust/src/guarantees.rs`: `compile_fail` doctests for the borrow-checked
  guarantees, each with a twin that compiles.
- Run locally (not in CI): the suite under ASan and under TSan with the C++
  instrumented too (`RUSTFLAGS=-Zsanitizer=... CXX=clang++
  CXXFLAGS=-fsanitize=... cargo +nightly test -Zbuild-std --target
  x86_64-unknown-linux-gnu --tests`), clean; coverage with `cargo llvm-cov`;
  `rust/fuzz` targets under ASan.
- CI (`rust.yml`): debug and release tests on Linux x86-64 and aarch64, macOS,
  Windows; `cargo package` (which builds the packaged crate), `cargo publish
  --dry-run`; fmt, clippy, docs as docs.rs builds them, and the MSRV.

## docs.rs

`build.rs` compiles no C++ when `DOCS_RS` is set (rustdoc does not link, and
docs.rs's compiler may lack `<experimental/simd>`); the library version still
comes through (`ISOSPEC_LIBRARY_VERSION`, read from `VERSION`).
