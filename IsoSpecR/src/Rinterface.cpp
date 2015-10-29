#include <Rcpp.h>
#include "cwrapper.h"
#include "misc.hpp"
#include "isoSpec++.hpp"
#include <vector>
#include <iostream>

using namespace Rcpp;

//' Call the C++ routine for calculating isotopic fine structure peaks.
//' 
//' @param isotopeNumbers 	Interger, numbers of isotopes of elements composing the studied chemical compound. 
//' @param atomCounts 		Interger, numbers of atoms of elements composing the studied chemical compound. 
//' @param isotopeMasses  	Numeric, masses of isotopes of elements composing the studied chemical compound. 
//' @param isotopeProbabilities Numeric, frequencies of isotopes of elements composing the studied chemical compound. 
//' @param cutOff 			Numeric, either the probability of the optimal p-set, or the threshold value of probability below which isotopologues are trimmed.
//' @param tabSize 			Integer, technical parameter.
//' @param hashSize 		Integer, initial size of the hash-table used.
//' @param step  			Numeric, the percentage of the percentile from the previous layer below defining the threshold for the next layer.
// [[Rcpp::export]]
List Rinterface(
	IntegerVector 	isotopeNumbers,
	IntegerVector 	atomCounts,
	NumericVector 	isotopeMasses,
	NumericVector 	isotopeProbabilities,
	double 			stopCondition,
	int				algo = 0,
	int 			tabSize 	= 1000,
	int				hashSize 	= 1000,
	double 			step 		= .25
){

	int dimNumber = isotopeNumbers.size();
	std::vector<int> stdIsotopeNumbers 	= Rcpp::as<std::vector<int> >(isotopeNumbers);
	std::vector<int> stdAtomCounts 		= Rcpp::as<std::vector<int> >(atomCounts);
	std::vector<double> stdIsotopeMasses 	= Rcpp::as<std::vector<double> >(isotopeMasses);
	std::vector<double> stdIsotopeProbabilities = Rcpp::as<std::vector<double> >(isotopeProbabilities);

	IsoSpec* iso = reinterpret_cast<IsoSpec*>(setupIso(
		dimNumber,
		stdIsotopeNumbers.data(),
		stdAtomCounts.data(),
		stdIsotopeMasses.data(),
		stdIsotopeProbabilities.data(),
		stopCondition,
		algo,
		tabSize,
		hashSize,
		step		
	));

	int confsNo = iso->getNoVisitedConfs();

	int isotopesNo = 0;
	for( auto it=isotopeNumbers.begin(); it != isotopeNumbers.end(); it++ ){
		isotopesNo += *it;
	}

	NumericVector res_mass(confsNo);
	NumericVector res_logProb(confsNo);
	IntegerVector res_isoCounts(confsNo*isotopesNo);

	int i = 0;
    int j = 0;

    for(auto it = iso->newaccepted.cbegin(); it != iso->newaccepted.cend(); it++)
    {
        int* curr_conf  = getConf(*it);
        res_mass[i] 	= combinedSum( curr_conf, iso->masses, dimNumber );
        res_logProb[i]	= getLProb(*it);

        for(int isotopeNumber=0; isotopeNumber<dimNumber; isotopeNumber++)
        {
            int currentConfIndex = curr_conf[isotopeNumber];
            int locIsoNo = isotopeNumbers[isotopeNumber];

            for( int k=0; k<locIsoNo; k++){
            	res_isoCounts[j+k] = (*(iso->marginalConfs[isotopeNumber]))[currentConfIndex][k];
            }
            j += locIsoNo;
        }
        i++;
    }

	destroyIso(iso);	

	List ret;

	ret["mass"] 	= res_mass; 
	ret["logProb"] 	= res_logProb;
	ret["configurations"] = res_isoCounts;

	return(ret);
} 
