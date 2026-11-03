# Changelog of the `isospec` crate

The crate's own versions; the C++ library it bundles has its own
(`isospec::LIBRARY_VERSION`), carried as semver build metadata after the `+`
(`0.1.0+2.6.0`: crate 0.1.0, IsoSpec 2.6.0) and listed with each release.

## Unreleased

- Envelopes from an `Iso` given by value: `Envelope::from_threshold`,
  `from_total_prob` and `binned` take `&Iso` or `Iso` (the new sealed trait
  `IsoArg`), and `Iso::into_envelope()` returns an `OwnedEnvelopeBuilder`. An
  owned `Iso` is moved into the computation instead of copied, which saves
  0.5-0.7 us per call for a 30-residue peptide (4.6 to 3.9 us binned at 1 Da).

## 0.1.0+2.6.0 (2026-10-26)

First release, bundling IsoSpec 2.6.0.

- Molecules (`Iso`) from formulas, peptide sequences with Unimod
  modifications, FASTA sequences, compositions or raw isotope data.
- Envelopes: threshold, total probability, stochastic and binned; arithmetic,
  statistics and Wasserstein distances.
- Generators: threshold, layered, ordered and stochastic, as iterators or
  step by step with configurations.
- The element table, compositions and Unimod tables.
