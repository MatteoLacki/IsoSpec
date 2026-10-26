# IsoSpecR 2.6.0

The package now follows the version numbers of the IsoSpec library it
contains (2.3.3 was the previous CRAN release).

## New

* `IsoBinned()`: a molecule's isotopic distribution on a mass grid, computed
  directly rather than by binning every isotopologue.
* `IsoStochastic()`: a simulated measurement of a given number of ions, with
  a `seed` for reproducible results (the same counts as IsoSpecPy gives for
  that seed).
* `IsoProperties()`: the monoisotopic, lightest, heaviest and most probable
  isotopologues, and the mean and variance of the mass.
* Operations on envelopes (the `mass`/`prob` matrices `IsoSpecify()`
  returns): `envelopeConvolve()`, `envelopeLinearCombination()`,
  `envelopeBin()`, `envelopeResample()`, `wassersteinDistance()`,
  `orientedWassersteinDistance()`, `abyssalWassersteinDistance()` and
  `wassersteinMatch()`.

## Changed

* Total-probability envelopes come out in the same order as from the other
  IsoSpec interfaces. Among equally probable isotopologues at the cut-off,
  which are kept may differ from earlier versions.
* Invalid input (NaN or out-of-range probabilities, isotope abundances not
  summing to 1, non-finite values) is an error rather than a wrong result.
* Faster total-probability and binned computations for large molecules.

## Fixed

* The package failed to load when only the base package was attached
  (as `R CMD check` loads it).
