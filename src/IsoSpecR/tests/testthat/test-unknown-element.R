context("IsoSpecify: unknown/unresolvable element symbols")

# github.com/MatteoLacki/IsoSpec/issues/49: a molecule element with zero
# matching rows in the isotopes table used to go through silently as a
# zero-isotope dimension -- wrong output (that atom's mass just missing, no
# error), and in at least one reported case, memory corruption severe enough
# to crash the R session on repeated calls. Should error clearly instead.

test_that("an element absent from the isotopes table errors, not returns mass 0", {
  expect_error(IsoSpecR::IsoSpecify(c(Tc = 1), 0.95), "Tc")
})

test_that("a real formula with one unresolvable element errors, not silently drops it", {
  formula_with_typo <- c(C = 18, H = 25, Ac = 1, N = 1, O = 5, S = 1)
  expect_error(IsoSpecR::IsoSpecify(formula_with_typo, 0.95), "Ac")
})

test_that("a normal, fully-resolvable formula is unaffected", {
  res <- IsoSpecR::IsoSpecify(c(H = 2, O = 1), 0.999)
  expect_true(nrow(res) > 0)
  expect_true(all(res[, "mass"] > 0))
})
