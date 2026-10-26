context("Binned and stochastic distributions, molecule properties, envelope operations")

# Expected values from IsoSpecPy 2.5.0 (the same C++ library), where noted.

sorted <- function(env) env[order(env[, "mass"]), , drop = FALSE]

test_that("IsoBinned matches IsoSpecPy's IsoBinned", {
  # IsoBinned(1.0, formula="C100H202", target_total_prob=0.999)
  bins <- sorted(IsoBinned(c(C = 100, H = 202), binWidth = 1, stopCondition = 0.999))
  expect_equal(nrow(bins), 7)
  expect_equal(unname(bins[1:3, "mass"]), c(1404, 1405, 1406))
  expect_equal(unname(bins[1:2, "prob"]), c(0.33020491691200565, 0.367830792929928), tolerance = 1e-12)
  expect_equal(sum(bins[, "prob"]), 0.9990856895148157, tolerance = 1e-12)
  # Bins of half a dalton, centred on half-daltons, and charge 2.
  half <- IsoBinned(c(C = 100, H = 202), binWidth = 0.5, binMiddle = 0.5, stopCondition = 0.999, charge = 2)
  expect_true(all(abs(half[, "mass"] * 2 - round(half[, "mass"] * 2)) < 1e-9))
  expect_equal(sum(half[, "prob"]), sum(bins[, "prob"]), tolerance = 1e-9)
  expect_error(IsoBinned(c(C = 1), binWidth = 0))
})

test_that("IsoStochastic with a seed gives IsoSpecPy's counts", {
  # IsoStochastic(1000, formula="C100", seed=42)
  spectrum <- sorted(IsoStochastic(c(C = 100), ionicCurrent = 1000, seed = 42))
  expect_equal(unname(spectrum[, "prob"]), c(329, 362, 202, 84, 17, 3, 3))
  expect_equal(unname(spectrum[1, "mass"]), 1200)
  again <- sorted(IsoStochastic(c(C = 100), ionicCurrent = 1000, seed = 42))
  expect_identical(again, spectrum)
  unseeded <- IsoStochastic(c(C = 100), ionicCurrent = 1000)
  expect_lte(sum(unseeded[, "prob"]), 1000)
  counted <- IsoStochastic(c(C = 10, H = 22), ionicCurrent = 100, seed = 1, showCounts = TRUE)
  expect_equal(colnames(counted), c("mass", "prob", "C12", "C13", "H1", "H2"))
  expect_true(all(rowSums(counted[, c("C12", "C13")]) == 10))
  expect_error(IsoStochastic(c(C = 1), ionicCurrent = 10, seed = -1))
  expect_error(IsoStochastic(c(C = 1), ionicCurrent = 10, seed = 1.5))
  expect_error(IsoStochastic(c(C = 1), ionicCurrent = 10.5))
  expect_error(IsoStochastic(c(C = 1), ionicCurrent = 10, precision = 2))
})

test_that("IsoProperties gives the distinguished peaks and moments", {
  glucose <- IsoProperties(c(C = 6, H = 12, O = 6))
  expect_equal(glucose$monoisotopicMass, 180.06338810844, tolerance = 1e-12)
  expect_equal(glucose$averageMass, 180.15669542527945, tolerance = 1e-12)
  expect_equal(glucose$variance, 0.11749298880254992, tolerance = 1e-10)
  expect_equal(glucose$modeMass, glucose$monoisotopicMass)
  expect_equal(glucose$lightestMass, glucose$monoisotopicMass)
  expect_equal(unname(glucose$heaviestCounts), c(0, 6, 0, 12, 0, 0, 6))
  expect_equal(names(glucose$monoisotopicCounts), c("C12", "C13", "H1", "H2", "O16", "O17", "O18"))
  expect_equal(exp(glucose$monoisotopicLogProb), 0.922, tolerance = 1e-3)
  expect_lt(glucose$heaviestLogProb, glucose$lightestLogProb)
  insulin <- IsoProperties(c(C = 254, H = 377, N = 65, O = 75, S = 6))
  expect_gt(insulin$modeMass, insulin$monoisotopicMass + 2)
  expect_gt(insulin$modeLogProb, insulin$monoisotopicLogProb)
  expect_error(IsoProperties(c(C = 6, H = -1)), "Negative atom count")
})

test_that("resample with a seed gives IsoSpecPy's counts", {
  # x = IsoTotalProb(0.999, formula="C6H12O6"); x.resample(1000, seed=7)
  glucose <- IsoSpecify(c(C = 6, H = 12, O = 6), 0.999)
  measured <- sorted(envelopeResample(glucose, ionicCurrent = 1000, seed = 7))
  measured <- measured[measured[, "prob"] > 0, , drop = FALSE]
  expect_equal(unname(measured[, "prob"]), c(925, 59, 1, 1, 8, 2, 4))
  expect_equal(unname(measured[1, "mass"]), 180.06338810844, tolerance = 1e-12)
  expect_error(envelopeResample(glucose, ionicCurrent = 10, betaBias = -1))
})

test_that("Wasserstein distances match IsoSpecPy's", {
  normalised <- function(env) { env[, "prob"] <- env[, "prob"] / sum(env[, "prob"]); env }
  a <- normalised(IsoSpecify(c(C = 6, H = 12, O = 6), 0.999))
  b <- normalised(IsoSpecify(c(C = 6, H = 12, O = 5, N = 1), 0.999))
  expect_equal(wassersteinDistance(a, b), 1.9930987487082414, tolerance = 1e-10)
  expect_equal(orientedWassersteinDistance(a, b), -1.9930987487082414, tolerance = 1e-10)
  expect_equal(abyssalWassersteinDistance(a, b, abyssDepth = 100), 1.9930987487082414, tolerance = 1e-10)
  m <- wassersteinMatch(a, b, flowDistance = 0.5)
  expect_equal(names(m), c("unmatchedThis", "unmatchedOther", "flow"))
  expect_equal(unname(m["flow"]), 0.011739440031703186, tolerance = 1e-10)
  expect_equal(unname(m["unmatchedThis"]), 0.988260559968297, tolerance = 1e-10)
  doubled <- a
  doubled[, "prob"] <- 2 * doubled[, "prob"]
  expect_error(wassersteinDistance(doubled, b))
  # Data frames work as well as matrices.
  expect_equal(wassersteinDistance(as.data.frame(a), as.data.frame(b)), 1.9930987487082414, tolerance = 1e-10)
})

test_that("convolution, linear combination and binning", {
  part <- cbind(mass = c(1, 2), prob = c(0.5, 0.5))
  whole <- envelopeConvolve(part, part)
  expect_equal(nrow(whole), 4)  # peaks of equal mass are not merged
  expect_equal(sum(whole[, "mass"] * whole[, "prob"]), 3)
  # Water twice is H4O2.
  water <- IsoSpecify(c(H = 2, O = 1), 0.99999)
  h4o2 <- IsoProperties(c(H = 4, O = 2))
  expect_equal(min(envelopeConvolve(water, water)[, "mass"]), h4o2$monoisotopicMass, tolerance = 1e-9)

  mix <- envelopeLinearCombination(list(part, cbind(mass = 100, prob = 1)), c(2, 0.5))
  expect_equal(unname(mix[, "prob"]), c(1, 1, 0.5))
  expect_equal(nrow(envelopeLinearCombination(list(part, part))), 4)
  expect_error(envelopeLinearCombination(list(part), c(1, 2)))

  binned <- sorted(envelopeBin(cbind(mass = c(100.2, 100.4, 101.1), prob = c(0.25, 0.25, 0.5)), width = 1))
  expect_equal(unname(binned), cbind(c(100, 101), c(0.5, 0.5)))
  # Binning computed envelopes agrees with IsoBinned.
  fine <- IsoSpecify(c(C = 100, H = 202), 0.999)
  expect_equal(sum(envelopeBin(fine)[, "prob"]), sum(fine[, "prob"]), tolerance = 1e-12)
})

test_that("bad envelopes are rejected", {
  expect_error(envelopeBin(cbind(m = 1, p = 1)), "columns")
  expect_error(envelopeBin(cbind(mass = NaN, prob = 1)))
  expect_error(envelopeBin(cbind(mass = 1, prob = -1)))
  expect_error(envelopeConvolve(cbind(mass = 1e300, prob = 1), cbind(mass = 1e300, prob = 1)))
})
