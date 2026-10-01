test_that("ordered isotope counts agree with their masses and composition", {
  isotopes <- IsoSpecR::RIsotopicTable()
  for (molecule in list(c(C = 2), c(C = 2, H = 2, O = 1))) {
    result <- IsoSpecR::IsoSpecify(molecule, 1, isotopes = isotopes,
                                  algo = 1, showCounts = TRUE)
    counts <- result[, -(1:2), drop = FALSE]
    rows <- match(colnames(counts), isotopes$isotope)
    expect_false(anyNA(rows))
    expect_equal(as.vector(counts %*% isotopes$mass[rows]),
                 as.vector(result[, "mass"]), tolerance = 1e-10)
    for (symbol in names(molecule)) {
      element_counts <- counts[, isotopes$element[rows] == symbol, drop = FALSE]
      expect_equal(as.vector(rowSums(element_counts)),
                   rep(unname(molecule[symbol]), nrow(result)))
    }
    expect_true(all(diff(result[, "prob"]) <= 0))
  }
})
