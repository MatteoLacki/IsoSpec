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
using namespace IsoSpec;

IsoGenerator* mkIsoG(Iso& iso, int algo, double stopCondition, int tabSize, int hashSize, int step, int trim)
{
    switch(algo)
    {
        case ISOSPEC_ALGO_LAYERED_ESTIMATE: // Not used anymore, just fall through to the next case
        case ISOSPEC_ALGO_LAYERED: return new IsoLayeredGenerator(std::move(iso), stopCondition, step, tabSize, hashSize, trim);
        case ISOSPEC_ALGO_ORDERED: return new IsoLayeredGenerator(std::move(iso), stopCondition, step, tabSize, hashSize, true); // Using the ordered algo in R is *completely* pointless today
                                                                                                                             // The only point of ordered algo is when one is using the generator
                                                                                                                             // interface, which we are not expising in R
                                                                                                                             // We'll just do layered, trim and sort it afterwards - it's equivalent 
                                                                                                                             // and much faster
        case ISOSPEC_ALGO_THRESHOLD_ABSOLUTE: return new IsoThresholdGenerator(std::move(iso), stopCondition, true, tabSize, hashSize, true);
        case ISOSPEC_ALGO_THRESHOLD_RELATIVE: return new IsoThresholdGenerator(std::move(iso), stopCondition, true, tabSize, hashSize, true);
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
    IsoGenerator* IG = mkIsoG(iso, algo, stopCondition, tabSize, hashSize, step, trim);

    unsigned int columnsNo = stdIsotopeTags.size(); // standard

    unsigned int isotopesNo = iso.getAllDim();

    // Code doing useless copying around of memory follows, as NumericMatrix apparently can't resize dynamically like std::vector does, so we can't directly
    // write into it, as we don't know how many configurations we're going to get upfront and we can't preallocate size.
    std::vector<double> logProbs;
    std::vector<double> masses;
    std::vector<int> confs;
    std::vector<size_t> ordering;
    size_t data_offset = 0;

    while(IG->advanceToNextConfiguration())
    {
        logProbs.push_back(IG->lprob());
        masses.push_back(IG->mass());
        if(showCounts)
        {
            confs.resize(confs.size()+isotopesNo);
            IG->get_conf_signature(&(confs.data()[data_offset]));
            data_offset += isotopesNo;
        }
    }

    delete IG;
    IG = nullptr;

    const bool needs_sorting = (ISOSPEC_ALGO_ORDERED == algo);

    if(needs_sorting)
    {
        // We need to sort the confs for backward compatibility
        ordering.reserve(logProbs.size());
        for(size_t i = 0; i < logProbs.size(); i++)
            ordering.push_back(i);
        std::sort(ordering.begin(), ordering.end(), [&logProbs](size_t idx1, size_t idx2) -> bool { return logProbs[idx1] > logProbs[idx2]; });
    }

    NumericMatrix res(logProbs.size(), columnsNo);


    size_t idx, confs_idx;

    for(size_t i = 0; i < logProbs.size(); i++)
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

    colnames(res) = stdIsotopeTags; //This is RCPP sugar. It sucks.

    return(res);
}
