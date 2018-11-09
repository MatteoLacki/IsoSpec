/*
 *   Copyright (C) 2015-2018 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */


#include <Rcpp.h>
#include "cwrapper.h"
#include "misc.h"
#include "isoSpec++.h"
#include <vector>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Rinterface(
    const IntegerVector&    molecule,
    const DataFrame&        isotopes,
    double                  stopCondition,
    int                     algo = 0,
    int                     tabSize = 1000,
    int                     hashSize = 1000,
    double                  step = .3,
    bool                    showCounts = false,
    bool                    trim = true
){

    unsigned int dimNumber = molecule.size();
    std::vector<int>     stdIsotopeNumbers;
    std::vector<double> stdIsotopeMasses;
    std::vector<double> stdIsotopeProbabilities;

    const CharacterVector& element = isotopes["element"];
    const CharacterVector& isotope = isotopes["isotope"];
    const NumericVector& mass        = isotopes["mass"];
    const NumericVector& abundance = isotopes["abundance"];
    const CharacterVector& molecule_names = molecule.names();

    CharacterVector stdIsotopeTags = CharacterVector::create("mass", "logProb");

    for (unsigned int i=0; i<dimNumber; i++)
    {
        unsigned int counter = 0;
        for (int j=0; j<element.size(); j++)
            if( element[j] == molecule_names[i] )
            {
                counter++;
                stdIsotopeMasses.push_back( mass[j] );
                stdIsotopeProbabilities.push_back( abundance[j] );
                if( showCounts )
                    stdIsotopeTags.push_back( isotope[j] );
            }
        stdIsotopeNumbers.push_back(counter);
    }

    IsoSpec* iso = reinterpret_cast<IsoSpec*>(setupIso(
        dimNumber,
        stdIsotopeNumbers.data(),
        Rcpp::as<std::vector<int> >( molecule).data(),
        stdIsotopeMasses.data(),
        stdIsotopeProbabilities.data(),
        stopCondition,
        algo,
        tabSize,
        hashSize,
        step,
        trim
    ));

    int confsNo = iso->getNoVisitedConfs();
    int columnsNo = stdIsotopeTags.size(); // standard

    int isotopesNo = 0;
    for( std::vector<int>::iterator it=stdIsotopeNumbers.begin(); it != stdIsotopeNumbers.end(); it++ ){
        isotopesNo += *it;
    }

    NumericMatrix res(confsNo, columnsNo);

    int i = 0;
    int j = 0;

    for( std::vector<void*>::iterator it = iso->newaccepted.begin(); it != iso->newaccepted.end(); it++)
    {
        int* curr_conf  = getConf(*it);
        res(i,0) = combinedSum( curr_conf, iso->masses, dimNumber );
        res(i,1) = getLProb(*it);

        j = 0;
        for(unsigned int isotopeNumber=0; isotopeNumber < dimNumber; isotopeNumber++)
        {
            int currentConfIndex = curr_conf[isotopeNumber];
            int locIsoNo = stdIsotopeNumbers[isotopeNumber];

            if( showCounts )
                for( int k=0; k<locIsoNo; k++)
                    res(i,2+j+k) = (*(iso->marginalConfs[isotopeNumber]))[currentConfIndex][k];
            j += locIsoNo;
        }
        i++;
    }

    colnames(res) = stdIsotopeTags; //This is RCPP sugar. It sucks.
    destroyIso(iso);

    return(res);
}
