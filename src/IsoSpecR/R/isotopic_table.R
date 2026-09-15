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

#' Read the core element/isotope table directly from the compiled library.
#'
#' Unlike \code{isotopicData$IsoSpec} (a static, hand-maintained copy), this
#' is built fresh from \code{element_tables.cpp} every time it's called --
#' the same authoritative source \code{IsoSpecPy.PeriodicTbl} reads on the
#' Python side. Excludes the synthetic charge-state pseudo-elements E/Me/Pn.
#'
#' @return A data.frame with columns \code{element}, \code{isotope},
#'   \code{mass}, \code{abundance}.
#' @export
#' @examples
#' library(IsoSpecR)
#' head(RIsotopicTable())
RIsotopicTable <- function() {
    .Call('_IsoSpecR_RIsotopicTable', PACKAGE = 'IsoSpecR')
}

# isotopicData's IsoSpec/IsoSpecShort/IsoSpecShortZero sub-tables used to be
# static, hand-maintained copies of the core element_tables.cpp data -- kept
# in TWO separate files (data/isotopicData.rda, the public data()-accessible
# copy, and R/sysdata.rda, the internal one IsoSpecify() actually reads at
# runtime), which had already silently drifted from each other before the
# #51 (deuterium) fix even started. Rebuilding them here from RIsotopicTable()
# on every package load makes the object IsoSpecify() actually uses
# impossible to go stale again, the same guarantee IsoSpecPy.PeriodicTbl
# already has on the Python side.
#
# enviPat/enviPatShort are untouched -- they're vendored verbatim from the
# external enviPat package (cross-validation reference data, not derived
# from our own core table), and stay exactly as sysdata.rda already has them
# (lazy-loaded automatically, same as before this change).
#
# data/isotopicData.rda (the old data()-accessible static copy) is gone --
# deleted, not just superseded -- so there is now exactly one isotopicData,
# this one, exported directly from the package namespace (see
# data_description.R's `@rawNamespace export(isotopicData)`; a bare string
# there is roxygen's doc-only convention, it creates no binding itself).
# `data(isotopicData)` no longer works and nothing should call it; plain
# `isotopicData` after library(IsoSpecR) is the only access path now, and it
# can no longer silently drift from what IsoSpecify() uses -- they're the
# same object. R/sysdata.rda still supplies the *initial* binding this
# function reads below (that's where enviPat/enviPatShort still live,
# lazy-loaded same as always) -- only its IsoSpec* fields are now always
# stale-on-disk-but-irrelevant, since this hook overwrites them in the
# namespace on every load regardless of what's on disk there.
.onLoad <- function(libname, pkgname) {
    ns <- asNamespace(pkgname)
    current <- get("isotopicData", envir = ns)

    fresh <- RIsotopicTable()
    fresh$ratioC <- NA_real_

    # "D" (deuterium) is included alongside the CHNOSSe set itself -- the
    # #51 fix added it to IsoSpecShort/IsoSpecShortZero too, not just the
    # full IsoSpec table, since it shares H's chemistry; preserved here.
    chnose <- c("H", "C", "N", "O", "S", "Se", "D")
    fresh_short <- fresh[fresh$element %in% chnose, ]
    # IsoSpecShortZero additionally carries one synthetic zero-abundance
    # placeholder row (a real radioactive sulfur isotope) -- preserved as a
    # template for a user adding a trace/synthetic isotope, matching what
    # was already in both static files before this change.
    fresh_short_zero <- rbind(
        fresh_short,
        data.frame(element = "S", isotope = "S35", mass = 35.0,
                   abundance = 0, ratioC = NA_real_, stringsAsFactors = FALSE)
    )

    current$IsoSpec <- fresh
    current$IsoSpecShort <- fresh_short
    current$IsoSpecShortZero <- fresh_short_zero

    assignInNamespace("isotopicData", current, ns = ns)
}
