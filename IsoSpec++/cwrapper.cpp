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



// =================================================================================

void* setupIsoLayered( int      _dimNumber,
                        const int*      _isotopeNumbers,
                        const int*      _atomCounts,
                        const double*   _isotopeMasses,
                        const double*   _isotopeProbabilities,
                        const double    _cutOff,
                        int             tabSize,
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
                                    _isotopeProbabilities, _StopCondition, tabSize, step, false, trim);
            break;
	case ALGO_LAYERED_ESTIMATE:
	    return setupIsoLayered(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses,
                                    _isotopeProbabilities, _StopCondition, tabSize, step, true, trim);

        case ALGO_ORDERED:
            return setupIsoOrdered(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses,
                                    _isotopeProbabilities, _StopCondition, tabSize, hashSize);
            break;
        case ALGO_THRESHOLD_ABSOLUTE:
            return NULL;

            break;
        case ALGO_THRESHOLD_RELATIVE:
            return NULL;

            break;
    }
    return NULL;
}


int getIsotopesNo(void* iso)
{
    return reinterpret_cast<IsoSpec*>(iso)->getNoIsotopesTotal();
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
