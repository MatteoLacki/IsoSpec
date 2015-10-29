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
#' \code{IsoSpecify} is a wrapper around \code{Rinterface} for the C++ implementation of the IsoSpec algorithm.
#' 
#' @param molecule  A named integer vector, e.g. \code{c(C=2,H=6,O=1)}, containing the chemical formula of the substance of interest.
#' @param stopCondition A numeric value between 0 and 1.
#' @param algo      An integer: 0 - use standard IsoStar algoritm, 
#' where \code{stopCondition} specifies the probability of the optimal p-set, 
#' 1 - use a version of algorithm that uses priority queue. Slower than 0, but does not require sorting.
#' 2 - use a threshold version of the algorithm, where \code{stopCondition} specifies the height of the pruned peaks. 
#' 3 - for the threshold version of IsoStar with \code{stopCondition} being 
#' the percentage of the highest peak below which isotopologues get pruned. 
#' @param isotopes  A named list of isotopic information required for IsoStar, e.g. \code{isotopicData$IsoSpecShortList}. The names must be valid element symbols. Each enlisted object should be a \code{data.frame} containing columns \code{element} (specifying the symbol of the element), \code{mass} (specifying the mass of the isotope), \code{abundance} (specyfying the assumed frequency of finding that isotope).
#' @param step      The percent of the the percentile of isotopologues in the current isolayer, specyfying the cutoff for the next isolayer. It has been optimised and better not change the default value.
#' @param tabSize   A technical parameter: the initial size of the \code{C++} dynamic table containing the results. Better not change the default value.
#' @return A list constaining the masses, logarithms of probability, and the tags of isotopes making up the molecule.     
#' @export 
#' @examples
#' res <- IsoSpecify( molecule = c(C=10,H=22,O=1), stopCondition = .9999 )
#' res2df(res)
IsoSpecify <- function(
        molecule,
        stopCondition,
        isotopes= isotopicData$IsoSpec,
        algo    = 0,
        step    = .25,
        tabSize = 1000
){
    molecule <- molecule[molecule>0]

    if( !all( names(molecule) %in% isotopes$element ) ) 
        stop(
            paste0('Elements: ',
                paste0(names(molecule)[!(names(molecule) %in% isotopes$element)],collapse=' ',''), 
                ' are not in the default/provided isotope data.frame. Check their name or insert an isotope data.frame containing this/these tags.'
            )
        )

    isotopesTmp <- isotopes[ 
        isotopes$element %in% names(molecule), 
        c('element','isotope','mass','abundance') 
    ] 

    molecule<- as.integer(molecule)
    dims    <- as.integer(table(isotopesTmp[,'element']))

    res <- Rinterface( isotopeNumbers = dims, atomCounts = molecule, isotopeMasses = isotopesTmp[,'mass'],
        isotopeProbabilities = isotopesTmp[,'abundance'], stopCondition = stopCondition, algo = as.integer(algo),
        tabSize = tabSize, hashSize = 1000, step = step )
    
    c(res, isotopes = list(isotopesTmp$isotope))
}

#' Represents the results of IsoSpec in a user friendly way.
#'
#' ATTENTION: this function sorts the results with descending probability. 
#' The sorting is a \code{nlog(n)} procedure which asymptotically dominates the calculations of IsoSpec and should be applied only for a good reason.
#' 
#' @param IsoSpecResults The list of results from calling \code{IsoSpecify}.
#' @return A \code{data.frame} containing the masses, logarithms of probability, and probability of the isotopologues.     
#' @export 
#' @examples
#' res <- IsoSpecify( molecule = c(C=10,H=22,O=1), stopCondition = .9999 )
#' res2df(res)
res2df <- function( IsoSpecResults ){
    confs <- as.data.frame( matrix(
            IsoSpecResults$configurations, 
            nrow = length(IsoSpecResults$mass), 
            byrow= TRUE ) )
    colnames(confs) <- IsoSpecResults$isotope
    res <- cbind( mass = IsoSpecResults$mass, logProb = IsoSpecResults$logProb, prob = exp(IsoSpecResults$logProb), cumProb = cumsum(exp(IsoSpecResults$logProb)), confs )
    res[order(res$logProb,decreasing=TRUE),]
}
