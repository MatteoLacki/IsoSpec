context("IsoSpecify: D (deuterium) as its own symbol")

# github.com/MatteoLacki/IsoSpec/issues/51

test_that("D is present in the default isotopes table as a pure, monoisotopic entry", {
  data(isotopicData, package = "IsoSpecR")
  d_rows <- isotopicData$IsoSpec[isotopicData$IsoSpec$element == "D", ]
  expect_equal(nrow(d_rows), 1)
  expect_equal(d_rows$abundance, 1.0)
  expect_equal(d_rows$mass, 2.014102)
})

test_that("heavy water (D2O) is heavier than normal water (H2O) by ~2 neutron masses", {
  heavy <- IsoSpecR::IsoSpecify(c(D = 2, O = 1), 0.999)
  normal <- IsoSpecR::IsoSpecify(c(H = 2, O = 1), 0.999)
  d_mono <- min(heavy[, "mass"])
  h_mono <- min(normal[, "mass"])
  expect_equal(d_mono - h_mono, 2 * (2.014102 - 1.007825), tolerance = 1e-5)
})
