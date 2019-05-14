# context("Numerical stability")
# 
# test_that("IsoSpecR calculates correctly full isotopic distributions.", {
#   expect_equal(nrow(IsoSpecify(c(C=100, H=202), 1)), 101*203)
#   expect_equal(nrow(IsoSpecify(c(C=100, H=202, S=5), 1)), 101*203*choose(4+5-1, 5))
#   expect_equal(nrow(IsoSpecify(c(O=100, N=10, S=6), 1)), choose(100+3-1,100)*choose(10+2-1, 10)*choose(4+6-1, 6))
# })
# 
# test_that("IsoSpecR provides results similar to envipat.", {
#   load("envipat.Rd")
#   isospec = lapply(mols, IsoSpecify, stopCondition=1)
#   expect_equal(sapply(isospec, nrow), sapply(envipat, nrow))
#   isospec = lapply(mols, IsoSpecify, stopCondition=1)
# 
#   for(i in 1:3){
#     envi_1_probs = envipat[[i]][,"abundance"] / sum(envipat[[i]][,"abundance"]) 
#     envi_1_probs = sort(envi_1_probs)
#     isos_1_probs = sort(exp(isospec[[i]][,2]))
#     expect_that(max(abs(isos_1_probs - envi_1_probs)) < 1.0e-10, is_true())
#   }
# })
