# Composition expansion: one isotope-resolution, reached from every language

## What this is

`composition.h`/`.cpp` holds the library's vocabulary for "a molecule as element counts" and the
one implementation of turning that into the isotope masses and probabilities an `Iso` is built
from:

- **`ElementComposition`** — parallel `(element_table_first_index, count)` arrays. Counts signed.
- **`ExpandedComposition`** — the flat per-isotope form the generic
  `Iso(dimNumber, isotopeNumbers, atomCounts, masses, probs)` constructor takes, which is also
  exactly what the C ABI's `setupIso` wants.
- **`expand_composition_into`** — the resolution itself.
- **`reject_negative_counts`** — the check every path from a composition to a molecule owes.
- **`build_iso_from_composition`** — the two of those, then the constructor.

Four callers reach `expand_composition_into`, and no other code resolves a symbol to its
isotopes:

| Caller | Route |
|---|---|
| Formula strings (C++) | `parse_formula` (`isoSpec++.cpp`) |
| Peptide sequences (C++) | `build_iso_from_composition` (`composition.cpp`) |
| C ABI / Python | `expandCompositionC` (`cwrapper.cpp`) → `IsoSpecPy.IsoParamsFromDict` |
| R | not yet — see "What is still duplicated" |

## The gap it filled

Before this, the C ABI offered exactly two ways to obtain an `Iso`, at opposite extremes:

- `setupIso` — the caller supplies *every* isotope mass and probability, i.e. must carry its
  own copy of the periodic table.
- `isoFromFasta` / `isoFromFastaWithMods` — accepts nothing but a peptide sequence.

A binding holding a *composition* (symbol → count) has neither. IsoSpecPy holds one on every
single `Iso(...)` call — from `ParseFormula`, from `ParsePeptideSequence`, or both merged — and
so it carried its own resolution loop in `IsoParamsFromDict`, walking `PeriodicTbl`'s
`symbol_to_masses` / `symbol_to_massNo` / `symbol_to_probs` dicts.

`expandCompositionC` is the missing middle rung. `IsoParamsFromDict` now forwards to it and
`PeriodicTbl` is no longer consulted for isotope resolution at all — only `ParseFormula` still
reads it, for its unknown-symbol check on a formula *string*.

## Why it needed its own header

`ElementComposition` and the expansion originally lived in `fasta_mods.h`. That worked for the
peptide path and the C ABI, but it meant the *formula* path could only share the resolution by
making the core include a peptide-modification header — backwards, so it didn't, and
`parse_formula` kept an identical loop of its own.

`composition.h` is the fix: three unrelated consumers (the formula parser, the sequence parsers,
the C ABI) include one small header that depends on nothing but `isoSpec++.h`. `fasta_mods.h`
includes it and keeps only the sequence-parsing declarations, so every existing includer still
compiles unchanged. `parse_formula` then collapsed onto it, because
`parse_formula_tokens`' two output vectors *are* an `ElementComposition`'s two members:

```cpp
ElementComposition composition;
parse_formula_tokens(formula, composition.element_first_index, composition.count);
reject_negative_counts(composition);
ExpandedComposition expanded;
expand_composition_into(composition, expanded, use_nominal_masses);
```

`parse_formula`'s mass/probability parameters are out-parameters of a public function that has
always *appended* to whatever the caller passed in, so the expansion is `insert`ed rather than
moved or assigned. `tests/C++/test_composition.cpp` pins that; it costs one copy of a few dozen
doubles on a path walked once per `Iso`.

## Why Python forwards here and not to `isoFromFastaWithMods`

`isoFromFastaWithMods` would collapse rather more Python than this does, and is still called by
nobody. It cannot be used, for three reasons, all of them about `Iso.__init__`'s public surface
rather than about sequences:

1. It takes a sequence and nothing else. `Iso(peptide_sequence=..., formula="H2O1")` — the
   documented way to get a neutral peptide's masses, since the sequence parsers deliberately
   omit the terminal water — composes a sequence *with* a formula. No C entry point does that.
2. `Iso.__init__` divides every mass by `charge` before `setupIso`. An entry point that
   resolves masses internally leaves nowhere to do it.
3. `self.atomCounts`, `self.isotopeMasses` and `self.isotopeNumbers` are public and read by
   `Advanced.py`, `approximations.py` and `tests/Python/test_all_configs_output.py`. If Python
   stops building them, they have to come back out of the `Iso` — four more C getters.

Expanding a composition sidesteps all three: Python keeps `setupIso`, keeps the charge
division, keeps the attributes, and loses only the duplicated lookup.

## Negative counts, and where they are rejected

`expand_composition_into` and `expandCompositionC` copy counts through **unchanged, negatives
included**. This follows the doctrine `ElementComposition`'s own comment states: a composition is
not a molecule, a Unimod delta can legitimately net a negative count (`G[UNIMOD:11]` → H−1 S−1),
and a second modification could bring it back up.

The non-negativity check therefore lives at each point a composition *becomes* a molecule.
There are four such points and two implementations:

- `reject_negative_counts` (`composition.cpp`) — shared by `build_iso_from_composition` and
  `parse_formula`, which used to have separate checks with separate messages. It names the
  offending elements, following R's example.
- `Iso.__init__` (Python) and `IsoSpecify` (R) keep their own, in their own idiom, because each
  has a better error to give than a C++ exception surfacing through a binding.

A negative atom count reaching `Iso`'s constructor is undefined behaviour, not merely a wrong
answer — `Marginal::computeModeConf()` sizes its buffer from the count and
`writeInitialConfiguration` then walks off the end (ASan catches it as a heap-buffer-overflow).

## What is still duplicated

One copy remains. **R's `isotopicData$IsoSpec`** (`src/IsoSpecR/data/isotopicData.rda`, 288 rows)
is an independent copy of the isotope table that `IsoSpecify` resolves symbols against. Verified
against `element_tables.cpp`: identical symbol order, max mass difference 0, max abundance
difference 1.1e-16 — in sync, with nothing enforcing it. The 4-row gap is the `E`/`Me`/`Pn`
pseudo-elements, which R's table omits.

Python solved this differently and better: `PeriodicTbl.py` *derives* its dicts from the C++
tables over cffi at import time, so it cannot drift. R could do the same — `Rinterface.cpp`
binds C++ directly, so it would call `expand_composition_into` rather than `expandCompositionC`
— but `isotopicData` is also documented, exported and user-overridable (`IsoSpecify(isotopes=)`),
so that is a user-facing change, not a refactor.

Separately, `IsoSpecPy.ParseFormula` is still a second *formula-grammar* parser (a regex)
alongside `parse_formula_tokens`. Different problem, untouched here.

## Build wiring

`composition.cpp` is a new translation unit, so it is registered in `unity-build.cpp` and
`src/IsoSpec++/Makefile`'s `SRCFILES`. The R package picks it up automatically (it compiles
every `.cpp` in `src/`, and its copies of the core are git-tracked **symlinks** into
`src/IsoSpec++/`, not vendored duplicates); the tests' `LIB_SRCS` is a wildcard; the wheel and
CMake builds go through `unity-build.cpp`; `pyproject.toml`'s sdist globs `*.cpp`/`*.h`.

## Tests

- `tests/C++/test_composition.cpp` — the expansion's contract (element order preserved, nominal
  masses, negatives passed through, reused out-parameter cleared), `reject_negative_counts`
  naming only the guilty elements, `build_iso_from_composition` against the equivalent formula,
  `parse_formula` agreeing with the composition path isotope for isotope (not merely in the
  masses it produces) and still appending to its out-parameters, and the C ABI round trip
  including `parseFastaWithModsC` → `expandCompositionC` → `setupIso`.
- `tests/Python/test_composition.py` — agreement with the old `PeriodicTbl` oracle for **every**
  element the library knows, the returned per-element tuple shape, and
  `test_resolution_no_longer_reads_periodic_tbl`, which empties `PeriodicTbl`'s dicts under
  `monkeypatch` and asserts everything still works. That last one is the guard: if a
  Python-side copy of the lookup ever returns, it is the only test that will notice.
