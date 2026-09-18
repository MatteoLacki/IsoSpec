# IsoSpec

IsoSpec is a fine-structure isotopic distribution calculator for chemical formulas — given a molecule, it computes the masses and probabilities of its isotopologues. It is fast enough to handle large molecules (proteins, oligonucleotides) where the full distribution has astronomically many configurations, and it returns provably-optimal subsets that cover a chosen fraction of the total probability mass.

IsoSpec is primarily used as a library by mass spectrometry software. It is implemented in C++ (`src/IsoSpec++/`) and shipped with first-class bindings for Python (`IsoSpecPy`, on PyPI) and R (`IsoSpecR`, on CRAN).

## Installation

### Python

```bash
pip install IsoSpecPy
```

The wheel bundles the C++ library — no separate native install required. Compatible with CPython and PyPy on Linux, macOS, Windows, and Cygwin/MinGW. Python ≥ 3.6.

To build from source: `pip install .` from a checkout. A C++20 compiler is required.

### R

```r
install.packages("IsoSpecR")
```

From source: `cd src && R CMD build IsoSpecR && R CMD INSTALL IsoSpecR_*.tar.gz`.

### C++ library

```bash
cd src/IsoSpec++
make            # produces libIsoSpec++.so
```

Or via CMake from the repo root:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
sudo cmake --install build
```

Requires a C++20 compiler. CMake produces both shared and static libraries and installs headers under `${CMAKE_INSTALL_INCLUDEDIR}/IsoSpec++`.

## Quick start

The most common usage is `IsoTotalProb`: ask for the smallest set of isotopologues whose summed probability covers a chosen fraction of the spectrum (typically 0.99–0.9999). The result is a materialized envelope — arrays of masses and probabilities you can read directly.

### Python

```python
import IsoSpecPy

# Isotopologues covering at least 99.9% of the probability mass of water.
iso = IsoSpecPy.IsoTotalProb(formula="H2O1", prob_to_cover=0.999)

for mass, prob in zip(iso.masses, iso.probs):
    print(mass, prob)

# From an amino-acid sequence:
iso = IsoSpecPy.IsoTotalProb(peptide_sequence="AAAPPGQAAC", prob_to_cover=0.999)
print(list(zip(iso.masses, iso.probs)))

# ...optionally carrying [UNIMOD:<id>] modification tags (the notation this
# monorepo's SAGE search-engine fork writes: [UNIMOD:<id>]-SEQUENCE prefix
# for an N-terminal mod, X[UNIMOD:<id>] inline for an internal one,
# SEQUENCE-[UNIMOD:<id>] suffix for a C-terminal one):
iso = IsoSpecPy.IsoTotalProb(peptide_sequence="PEPTC[UNIMOD:4]DEK", prob_to_cover=0.999)
print(list(zip(iso.masses, iso.probs)))

# Mind the water: peptide_sequence= yields the *residue* composition, without
# the terminal H2O of a free peptide, so the masses above are 18.0106 Da
# lighter than the neutral molecule. Add it explicitly for the real peptide
# (the C++ Iso::FromFASTAWithMods below adds it by default -- the two APIs
# differ here, and this is the older Python behaviour kept for compatibility):
iso = IsoSpecPy.IsoTotalProb(peptide_sequence="PEPTC[UNIMOD:4]DEK", formula="H2O1",
                             prob_to_cover=0.999)

# Or get just the composition (element -> atom count), no envelope:
composition = IsoSpecPy.ParsePeptideSequence("PEPTC[UNIMOD:4]DEK")
print(dict(composition))  # {'H': 60, 'C': 39, 'N': 10, 'O': 16, 'S': 1}
```

Resolves against a table of ~980 Unimod entries embedded in the library (isotope-labelled and glycan/derivatization "brick" modifications are excluded -- they aren't representable as a pure elemental delta). Pass `unimod_db_path="/path/to/table.csv"` to either function to resolve against a different table instead (same `id,name,mono_mass,composition` CSV shape as the embedded default).

See `Examples/Python/` for radiolabelling, custom elements, binned spectra, and FASTA/Unimod modifications.

### C++

```cpp
#include "IsoSpec++/isoSpec++.h"
#include "IsoSpec++/fixedEnvelopes.h"

using namespace IsoSpec;

int main() {
    FixedEnvelope iso = FixedEnvelope::FromTotalProb("H2O1", 0.999, true, true);

    for (size_t i = 0; i < iso.confs_no(); ++i) {
        std::cout << iso.mass(i) << '\t' << iso.prob(i) << '\n';
    }
}
```

From a peptide sequence, optionally carrying `[UNIMOD:<id>]` modification tags (`fasta_mods.h`; same notation as the Python/R examples above):

```cpp
#include "IsoSpec++/isoSpec++.h"
#include "IsoSpec++/fasta_mods.h"
#include "IsoSpec++/fixedEnvelopes.h"

using namespace IsoSpec;

int main() {
    FixedEnvelope iso = FixedEnvelope::FromTotalProb(
        Iso::FromFASTAWithMods("PEPTC[UNIMOD:4]DEK"), 0.999, true, true);

    for (size_t i = 0; i < iso.confs_no(); ++i) {
        std::cout << iso.mass(i) << '\t' << iso.prob(i) << '\n';
    }
}
```

`Iso::FromFASTAWithMods` takes an optional trailing `unimod_db_path` argument to resolve against a different Unimod table instead of the one embedded at compile time. For just the composition, no envelope, call `parse_fasta_with_mods`/`parse_fasta_with_mods_full` directly (or the allocation-free `_into` forms, for a hot loop over many sequences) — see `fasta_mods.h`.

Quickest way to build:

```bash
clang++ -std=c++20 water.cpp src/IsoSpec++/unity-build.cpp -o water
```

(`unity-build.cpp` is a single translation unit that `#include`s every source file — see `Examples/C++/COMPILING`.)

### R

```r
library(IsoSpecR)
water <- c(H = 2, O = 1)
IsoSpecify(molecule = water, stopCondition = 0.999)
```

`IsoSpecify` only ever takes a raw named-integer-vector `molecule` -- it has no sequence-string parameter -- so a `[UNIMOD:<id>]`-annotated peptide sequence takes two calls instead of one: `RParsePeptideSequence()` to resolve the sequence into a composition, then `IsoSpecify()` as above.

```r
composition <- RParsePeptideSequence("PEPTC[UNIMOD:4]DEK")
print(composition)  #  H  C  N  O  S
                     # 60 39 10 16  1
IsoSpecify(molecule = composition, stopCondition = 0.999)
```

Like the Python `ParsePeptideSequence`, this returns the residue composition without a free peptide's terminal H2O -- add `c(H = 2, O = 1)` to it for the neutral molecule. A modification whose delta removes more atoms of an element than the bare sequence provides gives a negative count (`"G[UNIMOD:11]"` nets `H -1`, `S -1`); `IsoSpecify` rejects such a vector rather than quietly dropping those elements.

See `Examples/R/` for radiolabelling, full-spectrum extraction, and Unimod-annotated sequences.

## Advanced features

- **Alternative algorithms** — for cases where `IsoTotalProb` isn't the right shape:
  - `IsoThreshold` — all isotopologues with probability above a fixed threshold.
  - `IsoStochastic` — simulate a measured spectrum by sampling integer ion counts.
  - `IsoBinned` — histogram-style envelope at a chosen bin width.
- **FASTA/peptide-sequence support** — build an `Iso` directly from an amino-acid sequence, with `[UNIMOD:<id>]` modification tags recognized inline; optionally include the N/C-terminal water.
- **Custom isotopic tables** — override natural abundances per element (e.g. for radio- or stable-isotope labelling). See `Examples/*/radiolabelling.*`.
- **Fixed-envelope arithmetic** — addition, normalization, convolution, Wasserstein distance.
- **Nominal-mass mode** — compute distributions over nucleon counts instead of real masses.
- **Streaming generators (for performance)** — the tabulated functions above materialize the full envelope into arrays. For large molecules where that uses too much memory, or for pipelines that consume one configuration at a time (e.g. binning on the fly), use the generator classes instead: `IsoThresholdGenerator`, `IsoLayeredGenerator`, `IsoOrderedGenerator`, `IsoStochasticGenerator`. They yield one isotopologue per call to `advanceToNextConfiguration()` without ever storing the full set. `IsoOrderedGenerator` additionally streams in strict order of decreasing probability.

Algorithmic details are in the papers cited below; the Supporting Information of each is the better starting point for implementation specifics.

## Repository layout

```
src/IsoSpec++/    C++ core library (the algorithms live here)
src/IsoSpecPy/    Python binding, loaded via cffi
src/IsoSpecR/     R binding, via Rcpp
Examples/         Working examples in each language
tests/            C, C++, and Python test suites
skbuild/          scikit-build-core entry point used when installing the Python wheel
```

## Citation

If IsoSpec is useful in your research, please cite:

- Łącki, M. K.; Valkenborg, D.; Startek, M. P. *IsoSpec2: Ultrafast Fine Structure Calculator.* Analytical Chemistry **2020**, *92* (14), 9472–9475. <https://doi.org/10.1021/acs.analchem.0c00959>
- Łącki, M. K.; Startek, M.; Valkenborg, D.; Gambin, A. *IsoSpec: Hyperfast Fine Structure Calculator.* Analytical Chemistry **2017**, *89* (6), 3272–3277. <https://doi.org/10.1021/acs.analchem.6b01459>

The Supporting Information of each paper has the algorithmic detail.

## License

2-clause BSD (see `LICENCE`). If you need different licensing terms, contact the authors. The authors appreciate (but do not require) being told when IsoSpec is used in other software.
