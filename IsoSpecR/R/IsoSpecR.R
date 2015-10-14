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
IsoSpecify <- function(
        molecule,
        stopCondition,
        tabSize = 1000,
        hashSize= 1000,
        isotopes= isotopicData$IsoSpecShortList,
        step    = .3,
	algo    = 0
){
        molecule        <- molecule[molecule>0]
        moleculeNames   <- names(molecule)

        if( !all( names(molecule) %in% names(isotopes) ) ) error('The elements you inserted have a non standard symbol. E.G. you should have inserted c(C=100,H=200) not c(Carbon=12,Hydrogen=423,Apple=12)')

        molecule <- as.integer(molecule)
        isotopes <- isotopes[moleculeNames]

        massAbundance   <- sapply(
                c('mass','abundance'),
                function(x) unlist(
                        sapply( isotopes, '[', x,USE.NAMES=FALSE ),
                        use.names=FALSE
                ),
                USE.NAMES=FALSE
        )

        dims            <- as.integer(sapply(isotopes, nrow))
        elementsNo      <- as.integer(length(dims))


        stupidRinterface(
                dims,
                molecule,
                as.double(massAbundance[,1]),
                as.double(massAbundance[,2]),
                stopCondition,
		algo,
                tabSize,
                hashSize,
                step
        )
}
