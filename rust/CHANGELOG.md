# Changelog of the `isospec` crate

The crate's own versions; the C++ library it bundles has its own
(`isospec::LIBRARY_VERSION`), listed with each release.

## 0.1.0 (unreleased)

First release, bundling IsoSpec 2.6.0.

- Molecules (`Iso`) from formulas, peptide sequences with Unimod
  modifications, FASTA sequences, compositions or raw isotope data.
- Envelopes: threshold, total probability, stochastic and binned; arithmetic,
  statistics and Wasserstein distances.
- Generators: threshold, layered, ordered and stochastic, as iterators or
  step by step with configurations.
- The element table, compositions and Unimod tables.
