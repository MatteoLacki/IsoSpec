context("IsoSpecify: D (deuterium) as its own symbol")

# github.com/MatteoLacki/IsoSpec/issues/51

test_that("D is present in the default isotopes table as a pure, monoisotopic entry", {
  d_rows <- isotopicData$IsoSpec[isotopicData$IsoSpec$element == "D", ]
  expect_equal(nrow(d_rows), 1)
  expect_equal(d_rows$abundance, 1.0)
  # Full core-table precision now (was a 6-decimal-rounded static value
  # before the storage-unification fix), so this needs a tolerance, not an
  # exact match against the old rounded literal.
  expect_equal(d_rows$mass, 2.01410177819, tolerance = 1e-6)
})

test_that("heavy water (D2O) is heavier than normal water (H2O) by ~2 neutron masses", {
  heavy <- IsoSpecR::IsoSpecify(c(D = 2, O = 1), 0.999)
  normal <- IsoSpecR::IsoSpecify(c(H = 2, O = 1), 0.999)
  d_mono <- min(heavy[, "mass"])
  h_mono <- min(normal[, "mass"])
  expect_equal(d_mono - h_mono, 2 * (2.014102 - 1.007825), tolerance = 1e-5)
})
