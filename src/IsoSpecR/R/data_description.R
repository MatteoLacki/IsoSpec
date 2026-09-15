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


#' Data on isotope masses, abundances and other.
#'
#' A list of data frames or table data frames (dplyr like), containing different information on isotopes.
#'
#' \code{isotopicData} is a live object refreshed every time the package is
#' loaded (see \code{.onLoad} in \code{isotopic_table.R}), not a
#' \code{data()}-loaded dataset -- just use it directly after
#' \code{library(IsoSpecR)}, no \code{data(isotopicData)} call needed (that
#' used to work here up to and including a stale copy that could silently
#' drift from what \code{IsoSpecify} actually used; removed together with
#' that possibility). \code{IsoSpec}, \code{IsoSpecShort} and
#' \code{IsoSpecShortZero} are (re)built on every load directly from the
#' compiled core element/isotope table (see \code{\link{RIsotopicTable}});
#' \code{enviPat} and \code{enviPatShort} are vendored, unchanged, from the
#' external \code{enviPat} package below.
#' @format A list of 5 data.frames, each containing:
#' \describe{
#' 	\item{element}{The symbol of an element from Mendeleev's periodic table.}
#' 	\item{isotope}{String composed of the nucleon number and the symbol of element.}
#' 	\item{mass}{Isotope's Mass in Daltons.}
#' 	\item{abundance}{The abundance of the isotopes. In case of enviPat data abundances do not sum to one. In case of all other, they do.}
#' 	\item{ratioC}{As in enviPat reference manual: "Maximum number of atoms of an element for one C-atom in a molecule, based on 99.99 \% of case molecules". Always \code{NA} outside the \code{enviPat}/\code{enviPatShort} entries -- not used anywhere in this package, kept only for column-shape consistency.}
#' }
#' @source R Package enviPat and Commission on Isotopic Abundances and Atomic Weights, CIAAW, \url{https://www.ciaaw.org/index.htm}, for \code{enviPat}/\code{enviPatShort}; this package's own compiled element/isotope table (see \code{\link{RIsotopicTable}}) for \code{IsoSpec}/\code{IsoSpecShort}/\code{IsoSpecShortZero}.
#' @usage isotopicData
#' @rawNamespace export(isotopicData)
"isotopicData"
