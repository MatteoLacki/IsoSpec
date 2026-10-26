# The R package (`IsoSpecR`)

Rcpp straight to the C++ classes (not the C ABI, which Python uses), from
`src/IsoSpecR/src/Rinterface.cpp`; the R functions are in
`src/IsoSpecR/R/`. `src/IsoSpecR/src` holds symlinks into `src/IsoSpec++`,
which `R CMD build` turns into copies.

## Shape of the API

- A molecule is a named integer vector (`c(C = 6, H = 12, O = 6)`) plus an
  isotope table (`isotopicData$IsoSpec` by default, rebuilt from the C++
  element table on load), as `IsoSpecify` has always taken it. Every
  molecule-level function (`IsoSpecify`, `IsoBinned`, `IsoStochastic`,
  `IsoProperties`) goes through the same `.check_molecule` (negative counts
  are an error) and the same C++ `iso_from_R`.
- An envelope is what `IsoSpecify` returns: a matrix with columns `mass` and
  `prob` (and isotope counts with `showCounts`). The envelope functions take
  such a matrix or a data frame, hand its two columns to C++
  (`FixedEnvelope::FromArrays`, which copies and checks them), and return a
  new matrix. R holds no C++ object between calls: no external pointers, no
  finalisers, nothing that a saved workspace could leave dangling.
- Left to plain R, deliberately: shifting or scaling masses, scaling or
  normalising probabilities, sorting, weighted means and variances. Each is
  one line on the matrix's columns, and a wrapper would only be more API.
- No generators: their state lives between calls, which the matrix-in,
  matrix-out design above avoids, and R code would step through them one
  isotopologue per interpreted iteration. The materialised functions cover
  what R users do with them.
- Seeds are numbers, `NULL` for the calling thread's randomly seeded
  generator; whole numbers up to 2^53 (exactly representable in a double),
  turned into `seeded_rng(uint64)`, so a seed gives the same stream as in
  C++, Python and Rust. Molecule counts (`ionicCurrent`) are checked the same
  way before the cast to `size_t`.

## Routine registration

`src/IsoSpecR_init.c` registers the `.Call` routines by hand
(`useDynLib(IsoSpecR)` without `.registration = TRUE`, so
`Rcpp::compileAttributes` does not generate it). A new `[[Rcpp::export]]`
needs a declaration and a `CallEntries` line there, with its argument count,
or R cannot find it.

## Quickselect pivots

`FromTotalProb`'s trim and `quickselect()` once took pivots from
`random_gen` in every build but R's, which took the middle element instead
(CRAN packages may not draw from their own random number generator). Since
the pivots come from `PivotSequence`, a fixed arithmetic sequence with no
RNG behind it, R uses them too: before, R's total-probability envelopes came
out in a different order from every other binding's (so a seeded
`envelopeResample` of the "same" envelope gave different counts), and a
middle pivot can go quadratic on adversarial input.

## Testing

testthat, `src/IsoSpecR/tests/testthat/`; `test-envelopes.R` checks the
newer functions against IsoSpecPy's values (binned probabilities, seeded
stochastic and resample counts, Wasserstein distances). The release workflow
runs `R CMD check` and `R CMD check --as-cran` on the built tarball and fails
on any WARNING or ERROR.
