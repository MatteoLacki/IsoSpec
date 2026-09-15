> Historical implementation plan. Current runtime CSV storage and packaging: see `docs/ai/unimod.md`.

# Unimod-aware peptide sequence parsing for IsoSpec

## Context

`git/isospec` (`MatteoLacki/IsoSpec`, cloned into this monorepo this session) computes
isotope distributions from elemental composition, but its only "peptide sequence" entry
point (`Iso::FromFASTA` / `parse_fasta`, `src/IsoSpec++/fasta.h`/`.cpp`) walks the string
one byte at a time through a fixed 256-entry CHNOSSe-only lookup table with no concept of
a bracket, digit, or colon. This repo's own SAGE fork (`git/sage`) already writes peptide
strings using `[UNIMOD:<id>]` modification notation (`crates/sage/src/peptide.rs:404-419`:
`[UNIMOD:<id>]-SEQUENCE` for N-term, `X[UNIMOD:<id>]` inline for internal, `-[UNIMOD:<id>]`
suffix for C-term) — and feeding one of those strings into IsoSpec today doesn't error, it
silently corrupts, because every letter in the word "UNIMOD" happens to be a valid 1-letter
amino acid code (confirmed empirically this session: `ParseFASTA('MPEPTC[UNIMOD:4]DEK')`
returns `C76 H123 N19 O27 S3 Se1` instead of the correct base composition plus
Carbamidomethyl's delta).

The goal: teach IsoSpec to recognize this exact notation, in C++, exposed to Python, so a
peptide string produced anywhere in this pipeline (e.g. `dump_peptides` output) can be fed
straight into IsoSpec without reformatting.

Decisions already made this session (do not re-litigate):
- **Scope trim**: ship only Unimod entries IsoSpec can represent correctly as a pure
  element-count delta — excludes 132/1560 isotope-labeled entries (`2H`, `13C`, …; IsoSpec
  has no pinned-isotope pseudo-elements) and 449/1560 glycan/derivatization "brick" entries
  (`Hex`, `HexNAc`, `dHex`, `NeuAc`, `Sulf`, `Ac`, …; each needs its own monosaccharide→formula
  expansion table, not built here). **979/1560 (~63%) ship**, verified to include every
  practically-important common PTM (Carbamidomethyl, Oxidation, Phospho, Acetyl, Methyl,
  GG/ubiquitin, Deamidated) and to require 27 distinct elements, not just CHNOSSe — so the
  new resolver must be general, not a fixed small array.
- **DB override**: no process-global flag/env var — pass an explicit `unimod_db_path` into
  the constructor/function at call time. Default (no path given) uses a table **embedded at
  compile time** (zero runtime I/O); an explicit path is parsed at runtime and cached per
  path (this pipeline processes millions of peptides, so re-parsing a ~979-row CSV on every
  call would be wasteful — cache keyed by path instead).
- **Naming**: `fasta=` has never accepted an actual FASTA file (no `>` header, no
  multi-record support) — only ever a bare sequence string. Keep `fasta=`/`ParseFASTA` as a
  compatibility alias; add a properly-named `peptide_sequence=`/`ParsePeptideSequence` as the
  recommended primary spelling; document plainly, including the user's requested aside that
  co-author Michał Startek "has completely no idea how things are supposed to be called."
- **Default behavior change**: `ParseFASTA`/`fasta=` now *always* scans for `[UNIMOD:<id>]`
  brackets (fixes the silent-corruption bug by default) — a literal `[`/digit was never valid
  input before either way, so nothing legitimate regresses. Flag prominently in docs/changelog.

## Useful discovery: half the low-level plumbing already exists

`parse_formula` (`src/IsoSpec++/isoSpec++.cpp:424-499`) already contains a working, general
"element symbol string → its isotopes in `element_tables.h`'s flat arrays" resolver (lines
461-475: linear scan of `elem_table_symbol[292]`; lines 480-493: walk forward while
`elem_table_ID[]` stays constant to collect that element's isotope run). This is exactly the
resolver the Unimod composition-delta path needs — extract it into a small reusable helper
rather than write a second copy. (Superseded below: rather than caching this scan in a map,
it gets replaced outright by a direct-indexed `[26][27]` table — see "Design".) Change
provably safe to make since `test_formula.cpp`'s existing suite is the regression gate.

Error handling: IsoSpec already has one clean convention worth reusing exactly, not
reinventing — core code throws `std::invalid_argument` (`isoSpec++.cpp:435`,
`marginalTrek++.cpp:131`), and every C ABI wrapper in `cwrapper.cpp` funnels through a
`c_guard()` helper (`cwrapper.cpp:45-67`, confirmed by direct read) that turns any escaped
exception into `nullptr`/`NaN`/`0`/no-op depending on return type. New wrappers should use
this, not a new error-signaling scheme.

## Design

**Unimod DB build script** — new `git/isospec/scripts/build_unimod_table.py`, same shape as
`git/sage/scripts/build_unimod_table.py` (fetch `unimod.obo`, regex over `[Term]` blocks),
extended to also capture `xref: delta_composition "..."` and apply the isotope/brick filter
above. Outputs two files in one run, always together (so they can't drift apart):
- `git/isospec/data/unimod.csv` — `id,name,mono_mass,composition` (composition = the Unimod
  `delta_composition` string **converted into IsoSpec's own native formula-string grammar**,
  e.g. `"H(3) C(2) N O"` → `"H3C2N1O1"`, `"H(-2) O(-1)"` → `"H-2O-1"` — so table loading can
  call the existing, minimally-extended `parse_formula` directly on this column, rather than
  a bespoke tokenizer for Unimod's own space/parenthesis notation).
- `git/isospec/src/IsoSpec++/unimod_table_data.h` — the same CSV embedded as a
  `R"UNIMODCSV(...)UNIMODCSV"` raw string literal, checked into git like `element_tables.cpp`
  already is (real periodic-table data also lives as a checked-in, generated-once C source
  file, not re-derived at build time — same precedent, no new CMake step needed: confirmed no
  `configure_file`/`xxd -i`/`#embed` pattern exists anywhere in this build today, and none is
  needed since `unity-build.cpp` already `#include`s sibling `.cpp` files directly and
  `skbuild`/sdist/`MANIFEST.in` already glob `*.h`).
  The script's own filter needs IsoSpec's known-element set — read it from the *installed*
  `IsoSpecPy.PeriodicTbl.symbol_to_masses.keys()` rather than hardcoding, so it can't drift
  from `element_tables.h`.
- A test asserts the embedded header's content matches `data/unimod.csv` read fresh, to catch
  the two ever going out of sync.

**C++ core** (new files, `fasta.h`/`.cpp` stay untouched — their fixed `int atomCounts[6]`
contract doesn't fit arbitrary-element deltas). Revised per feedback: no `std::map`/
`unordered_map`-based structures anywhere in this path — direct array indexing throughout,
and the *existing* formula-composition-reading code (`parse_formula`) is reused directly to
read each mod's atomic composition, not reimplemented as a separate bespoke tokenizer:
- `parse_formula` (`isoSpec++.cpp:424-499`) gets one small, additive extension: accept an
  optional leading `-` immediately before an element's digit run, so it can parse signed
  counts (`"H-2O-1"`). Existing formula strings never contain `-`, so this changes nothing for
  any current caller — purely additive, `test_formula.cpp` is the regression gate.
- The build script (below) converts each shipped entry's Unimod `delta_composition`
  (`"H(3) C(2) N O"`) into this exact native grammar with explicit counts and no separators
  (`"H3C2N1O1"`; `"H(-2) O(-1)"` → `"H-2O-1"`) at generation time — so loading a table entry's
  composition is a literal call to the same (now signed-capable) `parse_formula`, not a new
  parser.
- `element_lookup.h`/`.cpp`: **no map, no linear scan either** — a true direct-indexed lookup
  table, `int symbol_index_table[26][27]` (first letter A-Z × optional second letter a-z or
  "none"), built once at static-init by walking `elem_table_symbol`'s 292 rows and writing each
  element's first-isotope-row index into its slot. `find_element_table_first_index(symbol)`
  becomes two array indexes and a bounds/sentinel check, no scanning, no hashing.
  `element_isotope_count(first_index)` unchanged (walks `elem_table_ID` forward while constant,
  as `parse_formula` already does). Verified this is safe against IsoSpec's real data: all 87
  element symbols in `element_tables.cpp` are unique with zero collisions (checked directly,
  every symbol appears in exactly one contiguous isotope run) — a direct `[first][second]`
  table has no ambiguity to resolve.
- **Tokenization-boundary safety** (the real "two symbols" hazard, not a table-data collision):
  IsoSpec's native concatenated grammar (`"H3C2N1O1"`) relies on every element having an
  *explicit* trailing count digit — the scanner greedily consumes all consecutive letters as one
  element name, bounded on both sides by digits, so e.g. `"Co"` (Cobalt) can never be misread as
  `"C"` next to a stray `"o"`, and a 2-letter symbol like `Cl`/`Ca`/`Cr`/`Cu`/`Co` sitting next
  to plain `C` in the same composition is unambiguous *as long as no element's count is ever
  omitted* (Unimod's own notation frequently omits count-1, e.g. bare `"N"` — the build script
  must always emit it explicitly, e.g. `"N1"`, never bare `"N"`). This is a hard invariant on
  the generator, not a parser change; add a test with a real shipped composition that mixes a
  2-letter and an adjacent 1-letter symbol to prove the boundary holds.
- `unimod.h`/`.cpp`: `UnimodTable` is a **dense `std::vector<UnimodEntry>` indexed directly by
  Unimod id** (sized `max_id+1` among shipped entries), not a hashmap — matches the requested
  "make sure it is >=0 and < total number of entries in the table" bounds check exactly.
  `UnimodEntry{bool present; std::string name; double mono_mass; /* parse_formula's resolved
  isotopeNumbers/atomCounts/masses/probs for this entry's composition, computed once at table
  build/load time */}`. `embedded_unimod_table()` (lazy magic-static over
  `unimod_table_data.h`); `unimod_table_for_path(path)` — **not a map either**: a single cached
  `(last_path, table)` pair (process realistically uses 0 or 1 distinct override paths, so a
  full path-keyed map buys nothing), compared by string equality; single-threaded design per
  this repo's `CLAUDE.md`, so no locking needed — document that on the function.
- `fasta_mods.h`/`.cpp`: `parse_fasta_with_mods(sequence, mods)` — on seeing `[`, a direct
  `strncmp(p, "UNIMOD:", 7)` check (not a generic bracket-content grammar dispatch), then scans
  digits up to `]`, parses the id, bounds-checks `id < mods.size() && mods[id].present` (throws
  `std::invalid_argument` otherwise — covers unknown and deliberately-excluded ids alike, same
  error). Base-residue counts come from the existing `aa_symbol_to_elem_counts` table; each
  matched mod's pre-resolved composition is merged into a running **fixed-size accumulator
  array indexed by `element_lookup`'s element-table index** (small, dense — on the order of the
  ~30 distinct elements this feature ever touches, not a string-keyed map), converted to the
  final `isotopeNumbers`/`atomCounts`/`masses`/`probs` arrays for the generic
  `Iso(dimNumber, ...)` constructor `parse_formula` already uses.
  `Iso::FromFASTAWithMods(sequence, unimod_db_path=nullptr, ...)` ties it together, declared in
  `isoSpec++.h` next to `FromFASTA`.
- `unity-build.cpp` gets new `#include`s for the new `.cpp` files.

**C ABI + cffi**: `parse_fasta_c`/`isoFromFasta` stay byte-for-byte unchanged (no breaking
change to existing callers). New: `isoFromFastaWithMods`, and a handle+accessor group
(`parseFastaWithModsC` / `compositionSizeC` / `compositionSymbolsC` / `compositionCountsC` /
`deleteCompositionC`) mirroring the existing `FixedEnvelope` opaque-handle idiom rather than
the separate `...WithDeleter` idiom (that one exists for a different purpose — releasing
array ownership to the caller — not needed here since Python reads the result once into a
dict and discards it). New `isoFFI.py` cdef entries added by hand, mirroring `cwrapper.h`
exactly (existing convention — no cdef auto-generation exists today).

**Python surface** (`IsoSpecPy.py`): new `ParsePeptideSequence(sequence, unimod_db_path=None)`
replaces `ParseFASTA`'s internals (general arbitrary-symbol `OrderedDict`, superseding the old
fixed `["C","H","N","O","S"] + maybe "Se"` list entirely — no need to special-case elements
anymore, the new C ABI call already returns exactly what's present). `ParseFASTA` becomes a
one-line compatibility alias. `Iso.__init__` gains `peptide_sequence=""` and
`unimod_db_path=None` kwargs; `fasta`/`peptide_sequence` given non-default and unequal raises
`ValueError` (they're the same argument under two names); whichever is set flows through the
exact same `formula=`-summation code path that already exists (`IsoSpecPy.py:145-166`) — no
new summation logic needed, `Iso(formula=..., peptide_sequence=...)` falls out for free.

**Documentation**: new `git/isospec/docs/ai/unimod.md`, matching `git/sage/docs/ai/unimod.md`'s
shape and this monorepo's own AI-agent-doc convention — bracket placement/notation, how to
regenerate the DB, the `unimod_db_path` override + caching, the two explicit exclusions
(isotope-labeled, glycan-brick) with counts, and other known non-goals (no `[+mass]` support,
no real multi-record FASTA-file input). Plus the naming note in the docstrings of
`ParseFASTA`/`ParsePeptideSequence` themselves: `fasta=` never parsed an actual FASTA file,
kept only for backward compatibility, `peptide_sequence=` is the honest name — including the
requested remark that Michał Startek (co-author, `LICENCE`/`pyproject.toml`) has no idea how
things are supposed to be called. Add a "Documentation index" table to `git/isospec/CLAUDE.md`
pointing at the new doc file, without touching the pre-existing "Persona" section (its removal
is already in progress, uncommitted, on branch `experimental/unimod` from earlier this
session — separate, out of scope here).

**R**: checked `src/IsoSpecR/src/Rinterface.cpp` directly — R has **no** existing FASTA/sequence
entry point today, only a raw-atom-count-vector + isotope-DataFrame interface (`Rinterface`,
line 46). So R does not get this "for free"; it needs a genuine new
`// [[Rcpp::export]]` function (e.g. `RParsePeptideSequence(std::string sequence, std::string
unimod_db_path = "")`, returning an R named list/vector of symbol→count) calling
`parse_fasta_with_mods` directly — same shape as the new C ABI additions for Python, just
bound straight to the C++ classes the way `Rinterface.cpp` already does. Separately, this repo
contains no `Makevars`/copy-script showing *how* `src/IsoSpec++/*.cpp` reaches
`src/IsoSpecR/src/` before `R CMD build` — that out-of-repo step (owned elsewhere) needs to
also carry over the new files; document as a TODO in the new doc file rather than guessing at
a mechanism not visible in this repo.

## Critical files

- `git/isospec/scripts/build_unimod_table.py` (new)
- `git/isospec/data/unimod.csv` + `git/isospec/src/IsoSpec++/unimod_table_data.h` (new, generated)
- `git/isospec/src/IsoSpec++/element_lookup.h`/`.cpp` (new, extracted from `isoSpec++.cpp:461-493`)
- `git/isospec/src/IsoSpec++/unimod.h`/`.cpp` (new)
- `git/isospec/src/IsoSpec++/fasta_mods.h`/`.cpp` (new)
- `git/isospec/src/IsoSpec++/isoSpec++.cpp`/`.h` (refactor `parse_formula`'s resolver out; add `FromFASTAWithMods` declaration)
- `git/isospec/src/IsoSpec++/cwrapper.h`/`.cpp` (new C ABI entries)
- `git/isospec/src/IsoSpec++/unity-build.cpp` (new includes)
- `git/isospec/src/IsoSpecPy/isoFFI.py` (new cdef entries)
- `git/isospec/src/IsoSpecPy/IsoSpecPy.py` (`ParsePeptideSequence`, `Iso.__init__` kwargs)
- `git/isospec/src/IsoSpecR/src/Rinterface.cpp` (new `RParsePeptideSequence` export)
- `git/isospec/docs/ai/unimod.md` (new), `git/isospec/CLAUDE.md` (doc index entry)
- Tests: `git/isospec/tests/C++/test_unimod.cpp` (new, registered in `tests/C++/Makefile`'s `TEST_SRCS`), `git/isospec/tests/Python/test_iso_api.py` or a new `test_unimod.py`

Explicitly not touched: `fasta.h`/`.cpp` (unchanged), `parse_fasta_c`/`isoFromFasta` C ABI (unchanged signatures/behavior).

## Explicit non-goals

No glycan/brick support, no isotope-labeled mods, no `[+mass]` numeric notation, no real
multi-record FASTA-file parsing. All four excluded deliberately this round, not deferred
implicitly.

## Verification

1. `cd git/isospec/tests/C++ && make -j"$(nproc)" test` — new `test_unimod.cpp` cases (plain
   sequence unaffected, internal/N-term/C-term mod applies correct delta, unknown/excluded id
   throws, malformed bracket throws, path override works, path cache actually caches, embedded
   table matches `data/unimod.csv`) plus unchanged `test_fasta.cpp`/`test_formula.cpp` passing
   (regression gate for the `parse_formula` refactor).
2. `cd git/isospec/tests/C && make test` — new `cwrapper.h` entries still compile as valid C11.
3. `cd git/isospec && ./reinstall.sh` — rebuild/reinstall the Python wheel with the new sources.
4. `python -m pytest git/isospec/tests/Python` — new alias/summation/error/override tests pass.
5. Manual re-check of this session's own finding, expected to flip:
   ```python
   from IsoSpecPy.IsoSpecPy import ParseFASTA, ParsePeptideSequence
   print(ParseFASTA('MPEPTC[UNIMOD:4]DEK'))   # was: C76 H123 N19 O27 S3 Se1 (wrong)
   print(ParseFASTA('MPEPTCDEK'))             # baseline, unaffected
   print(ParsePeptideSequence('MPEPTC[UNIMOD:4]DEK'))  # same result under the new name
   ```
   Bracketed vs. unbracketed composition should now differ by exactly Carbamidomethyl's delta
   (`H3 C2 N1 O1`, ~57.021 Da), not by six phantom residues.
