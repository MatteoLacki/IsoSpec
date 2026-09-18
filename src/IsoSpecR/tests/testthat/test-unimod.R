context("RParsePeptideSequence: [UNIMOD:<id>] modification parsing")

test_that("a plain sequence with no brackets parses as before", {
  res <- RParsePeptideSequence("PEPTIDE")
  expect_true(all(c("C", "H", "N", "O") %in% names(res)))
  expect_true(all(res > 0))
})

test_that("an internal UNIMOD mod applies its composition delta", {
  plain <- RParsePeptideSequence("PEPTCDEK")
  modded <- RParsePeptideSequence("PEPTC[UNIMOD:4]DEK")  # Carbamidomethyl: H3 C2 N1 O1

  expect_equal(unname(modded["H"]) - unname(plain["H"]), 3L)
  expect_equal(unname(modded["C"]) - unname(plain["C"]), 2L)
  expect_equal(unname(modded["N"]) - unname(plain["N"]), 1L)
  expect_equal(unname(modded["O"]) - unname(plain["O"]), 1L)
})

test_that("an N-terminal UNIMOD mod applies before the first residue", {
  plain <- RParsePeptideSequence("PEPTIDE")
  modded <- RParsePeptideSequence("[UNIMOD:1]-PEPTIDE")  # Acetyl: H2 C2 O1

  expect_equal(unname(modded["H"]) - unname(plain["H"]), 2L)
  expect_equal(unname(modded["C"]) - unname(plain["C"]), 2L)
  expect_equal(unname(modded["O"]) - unname(plain["O"]), 1L)
})

test_that("a C-terminal UNIMOD mod applies after the last residue", {
  plain <- RParsePeptideSequence("PEPTIDE")
  modded <- RParsePeptideSequence("PEPTIDE-[UNIMOD:35]")  # Oxidation: O1

  expect_equal(unname(modded["O"]) - unname(plain["O"]), 1L)
})

test_that("an unknown or excluded UNIMOD id errors", {
  expect_error(RParsePeptideSequence("PEPTC[UNIMOD:999999]DEK"))
})

test_that("a malformed modification bracket errors", {
  expect_error(RParsePeptideSequence("PEPTC[UNIMOD:4DEK"))     # missing ']'
  expect_error(RParsePeptideSequence("PEPTC[FOO:4]DEK"))       # wrong prefix
  expect_error(RParsePeptideSequence("PEPTC[UNIMOD:]DEK"))     # missing numeric id
})

test_that("an empty unimod_db_path uses the embedded default table", {
  res <- RParsePeptideSequence("PEPTC[UNIMOD:4]DEK", unimod_db_path = "")
  expect_true("C" %in% names(res))
})

# A modification's delta can remove more atoms of an element than the bare
# sequence provides -- Met->Hsl (UNIMOD:11) is H-4 C-1 S-1, and a lone
# glycine cannot pay for it. RParsePeptideSequence reports the negative net
# (the arithmetic is correct and a further modification could bring it back
# up); IsoSpecify must refuse it rather than silently dropping those elements
# from every peak, which is what it did before.

test_that("a net-negative composition is reported by the parser", {
  composition <- RParsePeptideSequence("G[UNIMOD:11]")
  expect_true(any(composition < 0))
})

test_that("IsoSpecify refuses a negative atom count instead of ignoring it", {
  composition <- RParsePeptideSequence("G[UNIMOD:11]")
  expect_error(IsoSpecR::IsoSpecify(composition, 0.99), "[Nn]egative")
})

test_that("an id past 32 bits is unknown, not silently truncated to a real one", {
  expect_error(RParsePeptideSequence("PEPTC[UNIMOD:4294967300]DEK"))
})
