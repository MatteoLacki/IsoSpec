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
#' @format A list of 6 tbl_df's or data frames, each constaining:
#' \describe{
#' 	\item{element}{The symbol of an element from Mendeleev's periodic table.}	
#' 	\item{isotope}{String composed of the nucleon number and the symbol of element.}
#' 	\item{mass}{Isotope's Mass in Daltons.}
#' 	\item{abundance}{The abundance of the isotopes. In case of enviPat data abundances do not sum to one. In case of all other, they do.}
#' 	\item{ratioC}{As in enviPat reference manual: "Maximum number of atoms of an element for one C-atom in a molecule, based on 99.99 \% of case molecules".}
#' }
#' @source R Package enviPat and Commission on Isotopic Abundances and Atomic Weights, CIAAW, \url{https://www.ciaaw.org/index.htm}.
"isotopicData"
