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
#include "fixedEnvelopes.h"
#include <vector>
#include <iostream>

using namespace Rcpp;
using namespace IsoSpec;

FixedEnvelope mkIsoG(Iso& iso, int algo, double stopCondition, bool trim, bool get_confs)
{
    switch(algo)
    {
        case ISOSPEC_ALGO_LAYERED_ESTIMATE: // Not used anymore, just fall through to the next case
        case ISOSPEC_ALGO_LAYERED: return FixedEnvelope::FromTotalProb(std::move(iso), stopCondition, trim, get_confs);
        case ISOSPEC_ALGO_ORDERED: return FixedEnvelope::FromTotalProb(std::move(iso), stopCondition, true, get_confs); // Using the ordered algo in R is *completely* pointless today
                                                                                                                  // The only point of ordered algo is when one is using the generator
                                                                                                                  // interface, which we are not exposing in R
                                                                                                                  // We'll just do layered, trim and sort it afterwards - it's equivalent
                                                                                                                  // and much faster
        case ISOSPEC_ALGO_THRESHOLD_ABSOLUTE:
        case ISOSPEC_ALGO_THRESHOLD_RELATIVE: throw std::logic_error("");
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
    bool                    trim = true,
    double                  charge = 1.0
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

    CharacterVector stdIsotopeTags = CharacterVector::create("mass", "prob");

    for (unsigned int i=0; i<dimNumber; i++)
    {
        unsigned int counter = 0;
        for (int j=0; j<element.size(); j++)
            if( element[j] == molecule_names[i] )
            {
                counter++;
                stdIsotopeMasses.push_back( mass[j] / charge );
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

    unsigned int columnsNo = stdIsotopeTags.size(); // standard

    unsigned int isotopesNo = iso.getAllDim();

    if(algo == ISOSPEC_ALGO_THRESHOLD_ABSOLUTE || algo == ISOSPEC_ALGO_THRESHOLD_RELATIVE)
    {
        IsoThresholdGenerator ITG(std::move(iso), stopCondition, (algo == ISOSPEC_ALGO_THRESHOLD_ABSOLUTE));

        size_t no_confs = ITG.count_confs();
        size_t ii = 0;

        NumericMatrix res(no_confs, columnsNo);

        if(showCounts)
        {
            int* conf_sig = new int[isotopesNo];
            while(ITG.advanceToNextConfiguration())
            {
                res(ii,0) = ITG.mass();
                res(ii,1) = ITG.prob();
                ITG.get_conf_signature(conf_sig);
                for(size_t jj = 0; jj < isotopesNo; jj++)
                    res(ii, 2+jj) = conf_sig[jj];
                ii++;
            }
            delete[] conf_sig;
        }
        else
            while(ITG.advanceToNextConfiguration())
            {
                res(ii,0) = ITG.mass();
                res(ii,1) = ITG.prob();
                ii++;
            }

	colnames(res) = stdIsotopeTags; //This is RCPP sugar. It sucks.
        return res;
    }

    // The remaining (layered) algos
    FixedEnvelope TAB = mkIsoG(iso, algo, stopCondition, trim, showCounts);

    const double* probs = TAB.probs();
    const double* masses = TAB.masses();
    const int* confs = TAB.confs();
    std::vector<size_t> ordering;
    size_t confs_no = TAB.confs_no();


    const bool needs_sorting = (ISOSPEC_ALGO_ORDERED == algo);

    const unsigned int isotopesNoplus2 = isotopesNo + 2;

    NumericMatrix res(confs_no, columnsNo);

    if(needs_sorting)
    {
        // We need to sort the confs for backward compatibility
        ordering.reserve(confs_no);
        for(size_t i = 0; i < confs_no; i++)
            ordering.push_back(i);
        std::sort(ordering.begin(), ordering.end(), [&probs](size_t idx1, size_t idx2) -> bool { return probs[idx1] > probs[idx2]; });

        unsigned int idx;
        for(size_t i = 0; i < confs_no; i++)
        {
            idx = ordering[i];
            res(i,0) = masses[idx];
            res(i,1) = probs[idx];
        }

        if(showCounts)
        {
            unsigned int confs_idx;
            for(size_t i = 0; i < confs_no; i++)
            {
                confs_idx = ordering[i];
                for(size_t j = 2; j < isotopesNoplus2; j++)
                {
                    res(i,j) = confs[confs_idx];
                    confs_idx++;
                }
            }
        }
    }
    else
    {
        for(size_t i = 0; i < confs_no; i++)
        {
            res(i,0) = masses[i];
            res(i,1) = probs[i];
        }

        if(showCounts)
        {
            size_t confs_idx = 0;
            for(size_t i = 0; i < confs_no; i++)
            {
                for(size_t j = 2; j < isotopesNoplus2; j++)
                {
                    res(i,j) = confs[confs_idx];
                    confs_idx++;
                }
            }
        }
    }

    colnames(res) = stdIsotopeTags; //This is RCPP sugar. It sucks.

    return(res);
}
