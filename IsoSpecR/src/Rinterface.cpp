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
#include "tabulator.h"
#include <vector>
#include <iostream>

using namespace Rcpp;
using namespace IsoSpec;

Tabulator* mkIsoG(Iso& iso, int algo, double stopCondition, bool trim, bool get_confs)
{
    switch(algo)
    {
        case ISOSPEC_ALGO_LAYERED_ESTIMATE: // Not used anymore, just fall through to the next case
        case ISOSPEC_ALGO_LAYERED: return new LayeredTabulator(std::move(iso), stopCondition, trim, true, true, false, get_confs);
        case ISOSPEC_ALGO_ORDERED: return new LayeredTabulator(std::move(iso), stopCondition, true, true, true, false, get_confs); // Using the ordered algo in R is *completely* pointless today
                                                                                                                             // The only point of ordered algo is when one is using the generator
                                                                                                                             // interface, which we are not exposing in R
                                                                                                                             // We'll just do layered, trim and sort it afterwards - it's equivalent 
                                                                                                                             // and much faster
        case ISOSPEC_ALGO_THRESHOLD_ABSOLUTE: return new ThresholdTabulator(std::move(iso), stopCondition, false, true, true, false, get_confs);
        case ISOSPEC_ALGO_THRESHOLD_RELATIVE: return new ThresholdTabulator(std::move(iso), stopCondition, true, true, true, false, get_confs);
    }
    throw std::logic_error("Invalid algorithm selected");
}

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
    const NumericVector& mass      = isotopes["mass"];
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

    std::vector<double*> IM;
    std::vector<double*> IP;
    size_t tot = 0;
    for(size_t ii = 0; ii<stdIsotopeNumbers.size(); ii++)
    {
        IM.push_back(&(stdIsotopeMasses.data()[tot]));
        IP.push_back(&(stdIsotopeProbabilities.data()[tot]));
        tot += stdIsotopeNumbers[ii];
    }
    Iso iso(dimNumber, stdIsotopeNumbers.data(), Rcpp::as<std::vector<int> >( molecule).data(), IM.data(), IP.data());
    Tabulator* TAB = mkIsoG(iso, algo, stopCondition, trim, showCounts);

    unsigned int columnsNo = stdIsotopeTags.size(); // standard

    unsigned int isotopesNo = iso.getAllDim();

    double* logProbs = TAB->lprobs();
    double* masses = TAB->masses();
    int* confs = TAB->confs();
    std::vector<size_t> ordering;
    size_t data_offset = 0;
    size_t confs_no = TAB->confs_no();


    const bool needs_sorting = (ISOSPEC_ALGO_ORDERED == algo);

    if(needs_sorting)
    {
        // We need to sort the confs for backward compatibility
        ordering.reserve(confs_no);
        for(size_t i = 0; i < confs_no; i++)
            ordering.push_back(i);
        std::sort(ordering.begin(), ordering.end(), [&logProbs](size_t idx1, size_t idx2) -> bool { return logProbs[idx1] > logProbs[idx2]; });
    }

    NumericMatrix res(confs_no, columnsNo);

    size_t idx, confs_idx;

    for(size_t i = 0; i < confs_no; i++)
    {
        idx = needs_sorting ? ordering[i] : i;
        res(i,0) = masses[idx];
        res(i,1) = logProbs[idx];

        if(showCounts)
        {
            confs_idx = idx*isotopesNo;
            for(size_t j = 0; j < isotopesNo; j++)
                res(i,2+j) = confs[confs_idx+j];
        }
    }

    delete TAB;
    TAB = nullptr;

    colnames(res) = stdIsotopeTags; //This is RCPP sugar. It sucks.

    return(res);
}
