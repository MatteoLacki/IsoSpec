#
#   Copyright (C) 2015-2018 Mateusz Łącki and Michał Startek.
#
#   This file is part of IsoSpec.
#
#   IsoSpec is free software: you can redistribute it and/or modify
#   it under the terms of the Simplified ("2-clause") BSD licence.
#
#   IsoSpec is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the Simplified BSD Licence
#   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
#


#' @useDynLib  IsoSpecR
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("IsoSpecR", libpath)
}

#' Calculate the isotopic fine structure peaks.
#'
#' \code{IsoSpecify} is a wrapper around \code{Rinterface} that calls the C++ implementation of the IsoSpec algorithm. Given a molecular formula, it will calculate the smallest set of infinitely resolved peaks (isotopologues) that jointly is \code{p} probable, where \code{p} is provided by the user.
#'
#' @param molecule  A named integer vector, e.g. \code{c(C=2,H=6,O=1)}, containing the chemical formula of the substance of interest.
#' @param stopCondition A numeric value between 0 and 1.
#' @param showCounts Logical. If \code{TRUE}, then we output matrix contains additionally counts of isotopes for each isotopologue.
#' @param trim Logical. If \code{FALSE}, then we output matrix contains additionally isotopologues that otherwise would get trimmed in order to find the smalles possible p-set. Therefore, switching to \code{FALSE} results in a slightly larger set then the optimal p-set.
#' @param algo An integer: 0 - use standard IsoSpec algoritm,
#' where \code{stopCondition} specifies the probability of the optimal p-set,
#' 1 - use a version of algorithm that uses priority queue. Slower than 0, but does not require sorting.
#' 2 - use a threshold version of the algorithm, where \code{stopCondition} specifies the height of the pruned peaks.
#' 3 - for the threshold version of IsoSpec with \code{stopCondition} being
#' the percentage of the highest peak below which isotopologues get pruned.
#' @param isotopes  A named list of isotopic information required for IsoSpec. The names must be valid element symbols, see \code{isotopicData} for examples. Each enlisted object should be a \code{data.frame} containing columns \code{element} (specifying the symbol of the element), \code{mass} (specifying the mass of the isotope), \code{abundance} (specyfying the assumed frequency of finding that isotope).
#' @param step      The percent of the the percentile of isotopologues in the current isolayer, specyfying the cutoff for the next isolayer. It has been optimised and better not change the default value.
#' @param charge    The charge state of the molecule. All masses will be divided by this to obtain m/z values.
#' @return A numeric matrix containing the masses, the logarithms of probability, and, optionally, counts of isotopologues. Attention: this matrix does not have to be sorted. Sorting it would also compromise the linear complexity of our algorithm.
#' @export
#' @examples
#' library(IsoSpecR)
#' res <- IsoSpecify( molecule = c(C=10,H=22,O=1), stopCondition = .9999 )
#' print(res)
IsoSpecify <- function(
        molecule,
        stopCondition,
        isotopes= NULL,
        showCounts = FALSE,
        trim    = TRUE,
        algo    = 0,
        step    = .25,
        charge  = 1.0
){
    if(is.null(isotopes)){
        isotopes <- isotopicData$IsoSpec
    }

    Rinterface(
        molecule        = molecule[ molecule > 0 ],
        isotopes        = isotopes,
        stopCondition   = stopCondition,
        algo            = as.integer(algo),
        tabSize         = 1000,
        hashSize        = 1000,
        step            = step,
        showCounts      = showCounts,
        trim            = trim,
        charge          = charge
    )
}
