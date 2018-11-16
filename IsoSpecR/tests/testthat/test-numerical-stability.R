context("Numerical stability")

test_that("IsoSpecR calculates correctly full isotopic distributions.", {
  expect_equal(nrow(IsoSpecify(c(C=100, H=202), 1)), 101*203)
  expect_equal(nrow(IsoSpecify(c(C=100, H=202, S=5), 1)), 101*203*choose(4+5-1, 5))
  expect_equal(nrow(IsoSpecify(c(O=100, N=10, S=6), 1)), choose(100+3-1,100)*11*choose(4+6-1, 6))
})
