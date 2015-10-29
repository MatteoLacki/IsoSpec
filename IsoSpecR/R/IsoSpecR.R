#
#    This file is part of IsoSpec.
#
#    IsoSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License
#    version 3, as published by the Free Software Foundation.
#
#    IsoSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with IsoSpec.  If not, see <http://www.gnu.org/licenses/>.
#

#' @useDynLib  IsoSpecR
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("IsoSpecR", libpath)
}

# molecule <- c(C=2000,H=3000,N=300,O=300,S=20)
# cutOff <- .9999
# data(isotopicData)

#' Calculate the isotopic fine structure peaks.
#'
#' \code{IsoSpecify} is a wrapper around \code{Rinterface} for the C++ implementation of the IsoSpec algorith.
#' 
#' @param molecule  A named integer vector, e.g. \code{c(C=2,H=6,O=1)}, containing the chemical formula of the substance of interest.
#' @param stopCondition A numeric value between 0 and 1.
#' @param algo      An integer: 0 - for the standard IsoStar algoritm, where \code{stopCondition} specifies the probability of the optimal p-set, 3 - for the threshold version of IsoStar with \code{stopCondition} being the percentage of the highest peak below which isotopologues get pruned. 
#' @param isotopes  A named list of isotopic information required for IsoStar, e.g. \code{isotopicData$IsoSpecShortList}. The names must be valid element symbols. Each enlisted object should be a \code{data.frame} containing columns \code{element} (specifying the symbol of the element), \code{mass} (specifying the mass of the isotope), \code{abundance} (specyfying the assumed frequency of finding that isotope).
#' @param step      The percent of the the percentile of isotopologues in the current isolayer, specyfying the cutoff for the next isolayer. It has been optimised and better not change the default value.
#' @param tabSize   A technical parameter: the initial size of the \code{C++} dynamic table containing the results. Better not change the default value.
#' @return      
#' @export
IsoSpecify <- function(
        molecule,
        stopCondition,
        isotopes= isotopicData$IsoSpec,
        algo    = 0,
        step    = .3,
        tabSize = 1000
){
        # molecule <- c(C=10,H=20, O=2)
        molecule <- molecule[molecule>0]

        if( !all( names(molecule) %in% isotopes$element ) ) 
            stop(
                paste0('Elements: ',
                    paste0(names(molecule)[!(names(molecule) %in% isotopes$element)],collapse=' ',''), 
                    ' are not in the default/provided isotope data.frame. Check their name or insert an isotope data.frame containing this/these tags.'
                )
            )

        massAbundance <- as.matrix( isotopicData$IsoSpec[ isotopicData$IsoSpec$element %in% moleculeNames, c('mass','abundance') ] )
        rownames(massAbundance) <- NULL
        
        molecule <- as.integer(molecule)
        # isotopes <- isotopes[moleculeNames]
        # massAbundance   <- sapply(
        #         c('mass','abundance'),
        #         function(x) unlist(
        #                 sapply( isotopes, '[', x,USE.NAMES=FALSE ),
        #                 use.names=FALSE
        #         ),
        #         USE.NAMES=FALSE
        # )

        dims            <- as.integer(sapply(isotopes, nrow))
        elementsNo      <- as.integer(length(dims))


        Rinterface(
                dims,
                molecule,
                as.double(massAbundance[,'mass']),
                as.double(massAbundance[,'abundance']),
                stopCondition,
		        algo,
                tabSize,
                1000,
                step
        )
}
