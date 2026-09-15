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

#' Parse a peptide sequence with optional [UNIMOD:<id>] modification tags.
#'
#' Parses \code{sequence} for \code{[UNIMOD:<id>]} bracket notation in the
#' same placement this monorepo's SAGE search-engine fork writes: an
#' N-terminal modification as \code{[UNIMOD:<id>]-SEQUENCE}, an internal one
#' immediately after the modified residue as \code{X[UNIMOD:<id>]}, a
#' C-terminal one as \code{SEQUENCE-[UNIMOD:<id>]}. Unmodified letters are
#' read as plain amino acids, exactly as before this notation was supported.
#'
#' @param sequence A peptide sequence string, e.g. \code{"PEPTC[UNIMOD:4]DEK"}.
#' @param unimod_db_path Optional path to an override Unimod composition-delta
#'   CSV (same \code{id,name,mono_mass,composition} shape as the table
#'   shipped in the package). Empty string (the default) uses the packaged
#'   table.
#' @return A named integer vector of element symbol -> atom count, e.g.
#'   \code{c(C=44, H=69, N=11, O=17, S=2)}.
#' @export
#' @examples
#' library(IsoSpecR)
#' RParsePeptideSequence("PEPTC[UNIMOD:4]DEK")
RParsePeptideSequence <- function(sequence, unimod_db_path = "") {
    .Call('_IsoSpecR_RParsePeptideSequence', PACKAGE = 'IsoSpecR', sequence, unimod_db_path)
}
