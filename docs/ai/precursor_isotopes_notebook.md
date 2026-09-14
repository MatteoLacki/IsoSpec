# Pipeline precursor isotope notebook

Entry point: `notebooks/precursor_isotopes_99.ipynb` (JupyterLab, pipeline common kernel).
Checked notebook git blob: `57317ad16011fa1076b63847e82e841bdd27891a`.
Check freshness from the IsoSpec root with
`git hash-object notebooks/precursor_isotopes_99.ipynb`; refresh this summary and hash
when notebook behavior changes. Only this notebook belongs to this summary node.

## Input and calculations

- The actual checkout is `git/isospec` (lowercase).
- Input is the pipeline `dumped_peptides` artifact's `peptide` column, read through
  `pandas_ops.io.read_df`. Default points at
  `results/f9468_dump_peptides_recreate/dumped_peptides/peptides.parquet`.
- Default processes the first 1,000 of 6,300,884 input rows. Set
  `MAX_PEPTIDE_ROWS = None` for the full input. Column reading precedes slicing.
- Strip UNIMOD and signed numeric mass-shift brackets, then terminal hyphens.
  Reject unsupported notation and ambiguous residues; support canonical residues plus U/O.
- Deduplicate bare sequences in first-occurrence order; retain every selected source
  row's positional mapping, including target/decoy and modification variants.
- `IsoTotalProb(0.99, fasta=sequence, formula="H2O", get_minimal_pset=True)`:
  neutral free peptide, optimal fine-structure coverage, natural isotope abundances.
  IsoSpec's sequence parser returns residue composition, so H2O must be added explicitly.
- Round absolute exact masses via `floor(mass + 0.5)` and sum probabilities per
  occupied integer bin **after** selecting the optimal set. No renormalization.
  This is nearest 1 Da rounding, with half-integers rounded upward.

## Output

Each execution uses a fresh `notebooks/output/precursor_isotopes_99_*` directory.
`notebooks/.gitignore` excludes outputs and Jupyter checkpoints.

| Artifact | Columns / contents |
|---|---|
| `isotopes.mmappet` | sequence_id:uint64, mass_da:int64, probability:float64 |
| `source_rows.mmappet` | source_row:uint64, sequence_id:uint64 |
| `sequences.tsv` | sequence_id, bare sequence |
| `metadata.json` | input path/stat, counts, coverage, rounding, IsoSpec version |

mmappet tables are numeric; the sequence dictionary uses TSV to retain readable text.
Sparse rows sort by sequence ID then mass; append batches contain 1,000 sequences.
Metadata is written after export completes. A failed execution may leave partial files.
Input identity uses path/stat metadata, not a content checksum.

## Verification cells

Checks strip internal/terminal UNIMOD and numeric tags, reject unsupported strings,
exercise rounding boundaries and probability summation, and compare PEPTIDE against
C34H53N7O15. A minimal-set check verifies removing its least likely retained configuration
would drop coverage below 99%. Every exported sequence checks coverage and conservation.
Reopening mmappet verifies source mappings, schema, occupied-bin uniqueness, counts,
and probability sums against pre-write values. Full-input performance is not established.
