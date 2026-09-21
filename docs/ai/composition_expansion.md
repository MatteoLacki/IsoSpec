# Composition expansion: one isotope-resolution, reached from every language

## What this is

`expand_composition_into` (`src/IsoSpec++/fasta_mods.h`/`.cpp`) turns an
`ElementComposition` — parallel `(element_table_first_index, count)` arrays — into an
`ExpandedComposition`, the flat per-isotope form the generic
`Iso(dimNumber, isotopeNumbers, atomCounts, masses, probs)` constructor takes. `cwrapper.h`'s
`expandCompositionC` exposes it over the C ABI, keyed by element *symbol* rather than table
index, so a binding never needs a table index.

It is the single implementation of "resolve an element symbol to its isotopes' masses and
probabilities". Three callers reach it:

| Caller | Route |
|---|---|
| C++ | `build_iso_from_composition` (same file) |
| C ABI / Python | `expandCompositionC` → `IsoSpecPy.IsoParamsFromDict` |
| R | not yet — see "What is still duplicated" |

## The gap it fills

Before this, the C ABI offered exactly two ways to obtain an `Iso`, at opposite extremes:

- `setupIso` — the caller supplies *every* isotope mass and probability, i.e. must carry its
  own copy of the periodic table.
- `isoFromFasta` / `isoFromFastaWithMods` — accepts nothing but a peptide sequence.

A binding holding a *composition* (symbol → count) has neither. IsoSpecPy holds one on every
single `Iso(...)` call — from `ParseFormula`, from `ParsePeptideSequence`, or both merged — and
so it carried its own resolution loop in `IsoParamsFromDict`, walking `PeriodicTbl`'s
`symbol_to_masses` / `symbol_to_massNo` / `symbol_to_probs` dicts. That was a second
implementation of a rule `parse_formula` and `build_iso_from_composition` already owned, kept
correct only by both sides reading the same underlying tables.

`expandCompositionC` is the missing middle rung. `IsoParamsFromDict` now forwards to it and
`PeriodicTbl` is no longer consulted for isotope resolution at all — only
`ParseFormula` still reads it, for its unknown-symbol check on a formula *string*.

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
included**. This follows the doctrine `ElementComposition`'s own comment already states: a
composition is not a molecule, a Unimod delta can legitimately net a negative count
(`G[UNIMOD:11]` → H−1 S−1), and a second modification could bring it back up.

The non-negativity check therefore lives at each point a composition *becomes* a molecule, and
there are four of them, all still in place:

- `build_iso_from_composition` (C++) — moved above the expansion, not into it.
- `parse_formula` (C++) — its own copy, for the formula-string path.
- `Iso.__init__` (Python) — before `IsoParamsFromDict` is called.
- `IsoSpecify` (R) — `stop()` on `any(molecule < 0)`.

A negative atom count reaching `Iso`'s constructor is undefined behaviour, not merely a wrong
answer — `Marginal::computeModeConf()` sizes its buffer from the count and
`writeInitialConfiguration` then walks off the end (ASan catches it as a heap-buffer-overflow).
`parse_formula`'s comment has the long version.

## What is still duplicated

Two copies of this resolution remain, deliberately out of scope:

- **`parse_formula` (`isoSpec++.cpp`)** has the same loop over `element_isotope_count` and the
  same mass/probability pushes. It was left alone: it works on raw `int**` out-parameters with
  `array_copy` ownership, and routing the core formula path through `fasta_mods.h` would make
  the core depend on a peptide-modification header. The honest fix is to move
  `ElementComposition` / `ExpandedComposition` / `expand_composition_into` /
  `build_iso_from_composition` into their own `composition.h` that both include — a rename-shaped
  change touching `cwrapper.cpp`, R's `Rinterface.cpp` and `docs/ai/unimod.md`, worth doing on
  its own rather than smuggled in here.
- **R's `isotopicData$IsoSpec`** (`src/IsoSpecR/data/isotopicData.rda`, 288 rows) is an
  independent copy of the isotope table that `IsoSpecify` resolves symbols against. Verified
  against `element_tables.cpp` as of this change: identical symbol order, max mass difference
  0, max abundance difference 1.1e-16 — in sync, with nothing enforcing it. The 4-row gap is the
  `E`/`Me`/`Pn` pseudo-elements, which R's table omits. Python solved this differently and
  better: `PeriodicTbl.py` *derives* its dicts from the C++ tables over cffi at import time,
  so it cannot drift. R could call `expandCompositionC` the same way, but `Rinterface` binds
  C++ directly and would call `expand_composition_into` itself.

## Tests

- `tests/C++/test_composition.cpp` — the expansion's contract (element order preserved, nominal
  masses, negatives passed through, reused out-parameter cleared), `build_iso_from_composition`
  against the equivalent formula, and the C ABI round trip including
  `parseFastaWithModsC` → `expandCompositionC` → `setupIso`.
- `tests/Python/test_composition.py` — agreement with the old `PeriodicTbl` oracle for **every**
  element the library knows, the returned per-element tuple shape, and
  `test_resolution_no_longer_reads_periodic_tbl`, which empties `PeriodicTbl`'s dicts under
  `monkeypatch` and asserts everything still works. That last one is the guard: if a
  Python-side copy of the lookup ever returns, it is the only test that will notice.
