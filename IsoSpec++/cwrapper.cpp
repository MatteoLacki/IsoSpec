/*
 *   Copyright (C) 2015-2016 Mateusz Łącki and Michał Startek.
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


#include <tuple>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "cwrapper.h"
#include "misc.h"
#include "marginalTrek++.h"
#include "isoSpec++.h"


extern "C"
{

void* setupMarginal(
    const double* masses,   // masses size = logProbs size = isotopeNo
    const double* probs,
    int isotopeNo,                  // No of isotope configurations.
    int atomCnt,
    int tabSize,
    int hashSize
)
{
    /*
        *        std::tuple<double*,double*,int*,int> res_tmp =
        *                        getMarginal(
        *                                masses,
        *                                probs,
        *                                isotopeNo,
        *                                atomCnt,
        *                                cutOff,
        *                                tabSize,
        *                                hashSize
        *                        );
        *
        *        std::tuple<double*,double*,int*,int,int>* res =
        *                        new std::tuple<double*,double*,int*,int,int>(
        *                                                std::get<0>(res_tmp),
        *                                                std::get<1>(res_tmp),
        *                                                std::get<2>(res_tmp),
        *                                                std::get<3>(res_tmp),
        *                                                isotopeNo
        *                                            );
        *
        *        return reinterpret_cast<void*>(res);
        */
    MarginalTrek* MT = new MarginalTrek(
        masses,
        probs,
        isotopeNo,
        atomCnt,
        tabSize,
        hashSize
    );
    return reinterpret_cast<void*>(MT);
}

int probeConfigurationIdx(void* MT, int idx)
{
    return reinterpret_cast<MarginalTrek*>(MT)->probeConfigurationIdx(idx);
}

int getConfMT(void* MT, int idx, double* mass, double* logProb, int* configuration)
{
    MarginalTrek* IMT = reinterpret_cast<MarginalTrek*>(MT);
    if(!IMT->probeConfigurationIdx(idx))
        return 0;
    *mass = IMT->conf_masses()[idx];
    *logProb = IMT->conf_probs()[idx];
    memcpy(configuration, IMT->confs()[idx], IMT->_isotopeNo*sizeof(int));
    return 1;
}

int processMTUntilCutoff(void* MT, double cutoff)
{
    return reinterpret_cast<MarginalTrek*>(MT)->processUntilCutoff(cutoff);
}

int getConfNo(void* MT)
{
    return reinterpret_cast<MarginalTrek*>(MT)->conf_probs().size();
}

void getConfs(int howmany, void* MT, double* masses, double* logprobs, int* configurations)
{
    MarginalTrek* IMT = reinterpret_cast<MarginalTrek*>(MT);
    int cnt = std::min(howmany, (int)IMT->conf_probs().size());
    memcpy(masses, IMT->conf_probs().data(), cnt*sizeof(double));
    memcpy(logprobs, IMT->conf_probs().data(), cnt*sizeof(double));
    int dim = IMT->get_isotopeNo();
    for(int i=0; i<cnt; i++)
        memcpy(&configurations[i*dim], IMT->confs()[i], dim*sizeof(int));
}


void destroyConf(void* marginals)
{
    if (marginals != NULL)
    {
        delete reinterpret_cast<MarginalTrek*>(marginals);
    }
}


// =================================================================================

void* setupIsoLayered( int      _dimNumber,
                        const int*      _isotopeNumbers,
                        const int*      _atomCounts,
                        const double*   _isotopeMasses,
                        const double*   _isotopeProbabilities,
                        const double    _cutOff,
                        int             tabSize,
                        int             hashSize,
                        double          step,
                        bool            estimate,
                        bool            trim
)
{
    const double** IM = new const double*[_dimNumber];
    const double** IP = new const double*[_dimNumber];
    int idx = 0;
    for(int i=0; i<_dimNumber; i++)
    {
        IM[i] = &_isotopeMasses[idx];
        IP[i] = &_isotopeProbabilities[idx];
        idx += _isotopeNumbers[i];
    }


    IsoSpec* iso = new IsoSpecLayered(
        _dimNumber,
        _isotopeNumbers,
        _atomCounts,
        IM,
        IP,
        _cutOff,
        tabSize,
        hashSize,
        step,
	estimate,
	trim
    );

    try {
        iso->processConfigurationsUntilCutoff();
    }
    catch (std::bad_alloc& ba) {
        delete iso;
        iso = NULL;
    }

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}

void* setupIsoOrdered( int             _dimNumber,
                        const int*      _isotopeNumbers,
                        const int*      _atomCounts,
                        const double*   _isotopeMasses,
                        const double*   _isotopeProbabilities,
                        const double    _cutOff,
                        int             tabSize,
                        int             hashSize
)
{
    const double** IM = new const double*[_dimNumber];
    const double** IP = new const double*[_dimNumber];
    int idx = 0;
    for(int i=0; i<_dimNumber; i++)
    {
        IM[i] = &_isotopeMasses[idx];
        IP[i] = &_isotopeProbabilities[idx];
        idx += _isotopeNumbers[i];
    }


    IsoSpec* iso = new IsoSpecOrdered(
        _dimNumber,
        _isotopeNumbers,
        _atomCounts,
        IM,
        IP,
        _cutOff,
        tabSize,
        hashSize
    );

    try {
        iso->processConfigurationsUntilCutoff();
    }
    catch (std::bad_alloc& ba) {
        delete iso;
        iso = NULL;
    }

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}

void* setupIsoThreshold(    int      _dimNumber,
                            const int*      _isotopeNumbers,
                            const int*      _atomCounts,
                            const double*   _isotopeMasses,
                            const double*   _isotopeProbabilities,
                            const double    _threshold,
                            int             _absolute,
                            int             tabSize,
                            int             hashSize
)
{
    const double** IM = new const double*[_dimNumber];
    const double** IP = new const double*[_dimNumber];
    int idx = 0;
    for(int i=0; i<_dimNumber; i++)
    {
        IM[i] = &_isotopeMasses[idx];
        IP[i] = &_isotopeProbabilities[idx];
        idx += _isotopeNumbers[i];
    }


    IsoSpecThreshold* iso = new IsoSpecThreshold(
        _dimNumber,
        _isotopeNumbers,
        _atomCounts,
        IM,
        IP,
        _threshold,
        _absolute,
        tabSize,
        hashSize
    );

    try {
        iso->processConfigurationsAboveThreshold();
    }
    catch (std::bad_alloc& ba) {
        delete iso;
	iso = NULL;
    }

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}

void* setupIso( int             _dimNumber,
                const int*      _isotopeNumbers,
                const int*      _atomCounts,
                const double*   _isotopeMasses,
                const double*   _isotopeProbabilities,
                const double    _StopCondition,
                int             algo,
                int             tabSize,
                int             hashSize,
                double          step,
                bool            trim
)
{
    switch(algo)
    {
        case ALGO_LAYERED:
            return setupIsoLayered(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses,
                                    _isotopeProbabilities, _StopCondition, tabSize, hashSize, step, false, trim);
            break;
	case ALGO_LAYERED_ESTIMATE:
	    return setupIsoLayered(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses,
                                    _isotopeProbabilities, _StopCondition, tabSize, hashSize, step, true, trim);

        case ALGO_ORDERED:
            return setupIsoOrdered(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses,
                                    _isotopeProbabilities, _StopCondition, tabSize, hashSize);
            break;
        case ALGO_THRESHOLD_ABSOLUTE:
            return setupIsoThreshold(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses,
                                        _isotopeProbabilities, _StopCondition, true, tabSize, hashSize);
            break;
        case ALGO_THRESHOLD_RELATIVE:
            return setupIsoThreshold(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses,
                                        _isotopeProbabilities, _StopCondition, false, tabSize, hashSize);
            break;
    }
    return NULL;
}


int getIsotopesNo(void* iso)
{
    return reinterpret_cast<IsoSpec*>(iso)->getNoIsotopesTotal();
}

void* IsoFromFormula(const char* formula, double cutoff, int tabsize, int hashsize)
{
    IsoSpec* iso;
    try{
        iso = IsoSpec::IsoFromFormula<IsoSpecLayered>(
            formula, cutoff,
                tabsize, hashsize);
    }
    catch (const std::invalid_argument& e)
    {
        return NULL;
    }

    return reinterpret_cast<void*>(iso);
}

int getIsoConfNo(void* iso)
{
    return reinterpret_cast<IsoSpec*>(iso)->getNoVisitedConfs();
}

void getIsoConfs(void* iso, double* res_mass, double* res_logProb, int* res_isoCounts)
{
    reinterpret_cast<IsoSpec*>(iso)->getProduct(res_mass, res_logProb, res_isoCounts);
}

void destroyIso(void* iso)
{
    if (iso != NULL)
    {
        delete reinterpret_cast<IsoSpec*>(iso);
    }
}


}
