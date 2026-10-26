#
#   Copyright (C) 2015-2026 Mateusz Łącki and Michał Startek.
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

# Binned and stochastic distributions, molecule properties, and operations on
# envelopes. An envelope here is what IsoSpecify() returns: a matrix (or a
# data frame) with columns "mass" and "prob", one row per peak. Shifting
# masses, scaling or normalising probabilities, sorting and weighted means are
# ordinary R on those columns, so there are no functions for them.

.envelope_columns <- function(x, argument) {
    if(!all(c("mass", "prob") %in% colnames(x)))
        stop(argument, " must be a matrix or data frame with columns \"mass\" and \"prob\".")
    list(mass = as.numeric(x[, "mass"]), prob = as.numeric(x[, "prob"]))
}

.molecule_counts <- function(molecule) {
    .check_molecule(molecule)
    molecule[molecule > 0]
}

#' Isotopic distribution on a mass grid.
#'
#' Computes the isotopic distribution of a molecule summed into bins of width
#' \code{binWidth}, one bin centred on \code{binMiddle} (so the default bins
#' are centred on whole daltons). The isotopologues are enumerated, most
#' probable first, until \code{stopCondition} of the probability is covered,
#' and summed into their bins as they come: no envelope of individual
#' isotopologues is ever stored, which makes this far cheaper than binning the
#' output of \code{\link{IsoSpecify}} for large molecules.
#'
#' @param molecule A named integer vector, e.g. \code{c(C=2,H=6,O=1)}, the chemical formula.
#' @param binWidth The width of a bin, in daltons (of m/z, with \code{charge}).
#' @param binMiddle The centre of one of the bins; the others follow at multiples of \code{binWidth}.
#' @param stopCondition The probability to cover, between 0 and 1.
#' @param isotopes Isotope data, as in \code{\link{IsoSpecify}}; the default is \code{isotopicData$IsoSpec}.
#' @param charge The charge state: masses are divided by it.
#' @return A numeric matrix with columns \code{mass} (the bin centres) and \code{prob}, one row per nonempty bin, in no particular order.
#' @export
#' @examples
#' library(IsoSpecR)
#' bins <- IsoBinned(c(C=6, H=12, O=6), binWidth = 1)
#' bins[order(bins[, "mass"]), ]
IsoBinned <- function(molecule, binWidth = 1.0, binMiddle = 0.0, stopCondition = 0.9999,
                      isotopes = NULL, charge = 1.0) {
    RBinned(.molecule_counts(molecule), .default_isotopes(isotopes),
            stopCondition, binWidth, binMiddle, charge)
}

#' Simulated measurement of a molecule's isotopic distribution.
#'
#' Samples \code{ionicCurrent} molecules from the isotopic distribution and
#' reports how many fell on each isotopologue: the spectrum an instrument
#' counting that many ions would see. The isotopologues are not enumerated
#' first, so this works for molecules whose distributions are far too large to
#' list.
#'
#' @param molecule A named integer vector, e.g. \code{c(C=2,H=6,O=1)}, the chemical formula.
#' @param ionicCurrent The number of molecules to sample.
#' @param precision The share of the distribution sampled from, in (0, 1]; molecules that would fall in the remaining tail are dropped.
#' @param betaBias Tuning of the sampling algorithm (nonnegative); it changes the speed, not the distribution of the result.
#' @param seed \code{NULL} for a random result, or a whole number for a reproducible one: the same seed gives the same counts on every platform, and in IsoSpecPy and the C++ and Rust libraries.
#' @param isotopes Isotope data, as in \code{\link{IsoSpecify}}.
#' @param showCounts Logical: also return each isotopologue's isotope counts, as \code{\link{IsoSpecify}} does.
#' @param charge The charge state: masses are divided by it.
#' @return A numeric matrix with columns \code{mass} and \code{prob}, the latter the number of molecules at that isotopologue (and the isotope counts, with \code{showCounts}).
#' @export
#' @examples
#' library(IsoSpecR)
#' spectrum <- IsoStochastic(c(C=100, H=202), ionicCurrent = 1000, seed = 42)
#' sum(spectrum[, "prob"])
IsoStochastic <- function(molecule, ionicCurrent, precision = 0.9999, betaBias = 5.0, seed = NULL,
                          isotopes = NULL, showCounts = FALSE, charge = 1.0) {
    RStochastic(.molecule_counts(molecule), .default_isotopes(isotopes),
                ionicCurrent, precision, betaBias, seed, showCounts, charge)
}

#' Properties of a molecule's isotopic distribution.
#'
#' The masses and probabilities of the distinguished isotopologues -- the
#' monoisotopic (every atom its element's most abundant isotope), the lightest,
#' the heaviest and the most probable (the mode) -- and the mean and variance
#' of the mass, computed without enumerating the distribution.
#'
#' @param molecule A named integer vector, e.g. \code{c(C=2,H=6,O=1)}, the chemical formula.
#' @param isotopes Isotope data, as in \code{\link{IsoSpecify}}.
#' @param charge The charge state: masses are divided by it.
#' @return A named list: \code{monoisotopicMass}, \code{lightestMass}, \code{heaviestMass} and \code{modeMass}, with the natural logarithm of each one's probability (\code{monoisotopicLogProb}, ...) and, but for the mode, its isotope counts (\code{monoisotopicCounts}, ..., named integer vectors); \code{averageMass} and \code{variance}.
#' @export
#' @examples
#' library(IsoSpecR)
#' glucose <- IsoProperties(c(C=6, H=12, O=6))
#' glucose$monoisotopicMass
#' glucose$averageMass
IsoProperties <- function(molecule, isotopes = NULL, charge = 1.0) {
    RIsoProperties(.molecule_counts(molecule), .default_isotopes(isotopes), charge)
}

#' Convolution of two isotopic distributions.
#'
#' The distribution of the sum of the masses: that of a molecule made of the
#' two parts. Peaks of equal mass are not merged.
#'
#' @param x,y Envelopes: matrices or data frames with columns \code{mass} and \code{prob}, as returned by \code{\link{IsoSpecify}}.
#' @return An envelope: a numeric matrix with columns \code{mass} and \code{prob}.
#' @export
#' @examples
#' library(IsoSpecR)
#' water <- IsoSpecify(c(H=2, O=1), 0.999)
#' two <- envelopeConvolve(water, water)
#' sum(two[, "prob"])
envelopeConvolve <- function(x, y) {
    a <- .envelope_columns(x, "x")
    b <- .envelope_columns(y, "y")
    REnvelopeConvolve(a$mass, a$prob, b$mass, b$prob)
}

#' Weighted sum of isotopic distributions.
#'
#' All the envelopes' peaks together, each envelope's probabilities multiplied
#' by its intensity: a mixture.
#'
#' @param envelopes A list of envelopes (matrices or data frames with columns \code{mass} and \code{prob}).
#' @param intensities The weight of each envelope.
#' @return An envelope: a numeric matrix with columns \code{mass} and \code{prob}.
#' @export
#' @examples
#' library(IsoSpecR)
#' light <- IsoSpecify(c(C=6, H=12, O=6), 0.999)
#' heavy <- IsoSpecify(c(C=6, H=12, O=6), 0.999, charge = 2)
#' mixture <- envelopeLinearCombination(list(light, heavy), c(0.7, 0.3))
envelopeLinearCombination <- function(envelopes, intensities = rep(1, length(envelopes))) {
    columns <- lapply(seq_along(envelopes), function(i) .envelope_columns(envelopes[[i]], "each envelope"))
    REnvelopeLinearCombination(lapply(columns, `[[`, "mass"), lapply(columns, `[[`, "prob"),
                               as.numeric(intensities))
}

#' Bin an isotopic distribution.
#'
#' Sums the peaks into bins of width \code{width}, one bin centred on
#' \code{middle}; each bin is reported at its centre.
#'
#' @param x An envelope: a matrix or data frame with columns \code{mass} and \code{prob}.
#' @param width The width of a bin.
#' @param middle The centre of one of the bins.
#' @return An envelope: a numeric matrix with columns \code{mass} and \code{prob}, one row per nonempty bin.
#' @seealso \code{\link{IsoBinned}}, which bins a molecule's distribution without computing it peak by peak first.
#' @export
#' @examples
#' library(IsoSpecR)
#' fine <- IsoSpecify(c(C=6, H=12, O=6), 0.9999)
#' envelopeBin(fine, width = 1)
envelopeBin <- function(x, width = 1.0, middle = 0.0) {
    a <- .envelope_columns(x, "x")
    REnvelopeBin(a$mass, a$prob, width, middle)
}

#' Simulated measurement of an isotopic distribution.
#'
#' Replaces the probabilities by the counts of \code{ionicCurrent} molecules
#' sampled from them.
#'
#' @param x An envelope: a matrix or data frame with columns \code{mass} and \code{prob}.
#' @param ionicCurrent The number of molecules to sample.
#' @param betaBias Tuning of the sampling algorithm (nonnegative); it changes the speed, not the distribution of the result.
#' @param seed \code{NULL} for a random result, or a whole number for a reproducible one.
#' @return An envelope: a numeric matrix with columns \code{mass} and \code{prob}, the latter molecule counts.
#' @seealso \code{\link{IsoStochastic}}, which samples a molecule's distribution without computing it first.
#' @export
#' @examples
#' library(IsoSpecR)
#' envelope <- IsoSpecify(c(C=100), 0.999)
#' measured <- envelopeResample(envelope, ionicCurrent = 1000, seed = 1)
#' sum(measured[, "prob"])
envelopeResample <- function(x, ionicCurrent, betaBias = 1.0, seed = NULL) {
    a <- .envelope_columns(x, "x")
    REnvelopeResample(a$mass, a$prob, ionicCurrent, betaBias, seed)
}

#' Wasserstein distances between isotopic distributions.
#'
#' The Wasserstein (earth mover's) distance: how much probability must be
#' moved how far to turn one distribution into the other. \code{wassersteinDistance}
#' needs both normalised (probabilities summing to 1). \code{orientedWassersteinDistance}
#' is signed: positive when \code{y} lies to the heavier side of \code{x}.
#' \code{abyssalWassersteinDistance} also lets probability be moved to or from
#' an abyss, at cost \code{abyssDepth} per unit, so the distributions need not
#' have equal totals; \code{y}'s probabilities are multiplied by
#' \code{otherScale} first.
#'
#' @param x,y Envelopes: matrices or data frames with columns \code{mass} and \code{prob}.
#' @param abyssDepth The cost of moving a unit of probability into or out of the abyss.
#' @param otherScale The factor \code{y}'s probabilities are multiplied by.
#' @return The distance, a number.
#' @export
#' @examples
#' library(IsoSpecR)
#' x <- cbind(mass = 100, prob = 1)
#' y <- cbind(mass = 101, prob = 1)
#' wassersteinDistance(x, y)
#' orientedWassersteinDistance(y, x)
#' abyssalWassersteinDistance(x, cbind(mass = 200, prob = 1), abyssDepth = 1)
wassersteinDistance <- function(x, y) {
    a <- .envelope_columns(x, "x")
    b <- .envelope_columns(y, "y")
    RWasserstein(a$mass, a$prob, b$mass, b$prob, 0L, 0.0, 1.0)
}

#' @rdname wassersteinDistance
#' @export
orientedWassersteinDistance <- function(x, y) {
    a <- .envelope_columns(x, "x")
    b <- .envelope_columns(y, "y")
    RWasserstein(a$mass, a$prob, b$mass, b$prob, 1L, 0.0, 1.0)
}

#' @rdname wassersteinDistance
#' @export
abyssalWassersteinDistance <- function(x, y, abyssDepth, otherScale = 1.0) {
    a <- .envelope_columns(x, "x")
    b <- .envelope_columns(y, "y")
    RWasserstein(a$mass, a$prob, b$mass, b$prob, 2L, abyssDepth, otherScale)
}

#' Match the probability of two isotopic distributions.
#'
#' Moves probability between peaks of \code{x} and of \code{y} at most
#' \code{flowDistance} apart, \code{y}'s probabilities multiplied by
#' \code{otherScale} first, and reports how much was matched and how much of
#' each was left over.
#'
#' @param x,y Envelopes: matrices or data frames with columns \code{mass} and \code{prob}.
#' @param flowDistance The largest mass difference over which probability may be matched.
#' @param otherScale The factor \code{y}'s probabilities are multiplied by.
#' @return A named numeric vector: \code{unmatchedThis} (of \code{x}), \code{unmatchedOther} (of \code{y}, scaled) and \code{flow} (matched).
#' @export
#' @examples
#' library(IsoSpecR)
#' x <- cbind(mass = 100, prob = 1)
#' y <- cbind(mass = 100.001, prob = 0.5)
#' wassersteinMatch(x, y, flowDistance = 0.01)
wassersteinMatch <- function(x, y, flowDistance, otherScale = 1.0) {
    a <- .envelope_columns(x, "x")
    b <- .envelope_columns(y, "y")
    RWassersteinMatch(a$mass, a$prob, b$mass, b$prob, flowDistance, otherScale)
}
