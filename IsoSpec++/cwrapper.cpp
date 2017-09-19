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
#include "tabulator.h"


extern "C"
{

void* setupIso( int      _dimNumber,
                const int*      _isotopeNumbers,
                const int*      _atomCounts,
                const double*   _isotopeMasses,
                const double*   _isotopeProbabilities)
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

    Iso* iso = new Iso(
        _dimNumber,
        _isotopeNumbers,
        _atomCounts,
        IM,
        IP
    );

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);

}


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

// ATTENTION! BELOW THIS LINE MATTEO WAS CODING AND IT IS BETTER NOT TO COMPILE THAT
#define C_CODE(generatorType, dataType, method)\
dataType method##generatorType(void* generator){ return reinterpret_cast<generatorType*>(generator)->method(); }

#define DELETE(generatorType) void delete##generatorType(void* generator){ delete reinterpret_cast<generatorType*>(generator); }

// #define C_GENERATOR_TO_ARRAY_CODE(generatorType)\
// MassSpectrum set_tables##generatorType(void* generator, int init_size) \
// { \
//     DoubleArray Masses, LogProbs; \
//     DoubleArrayConstructor(&Masses, init_size); \
//     DoubleArrayConstructor(&LogProbs, init_size); \
//     \
//     while(advanceToNextConfiguration##generatorType(generator)) \
//     { \
//         DoubleArrayPush(&Masses, mass##generatorType(generator)); \
//         DoubleArrayPush(&LogProbs, lprob##generatorType(generator)); \
//     } \
//     \
//     MassSpectrum S = { Masses.array, LogProbs.array, Masses.used }; \
//     DoubleArrayDestructor(&Masses, true); \
//     DoubleArrayDestructor(&LogProbs, true); \
//     return S; \
// }\

#define C_CODES(generatorType)\
C_CODE(generatorType, double, mass) \
C_CODE(generatorType, double, lprob) \
C_CODE(generatorType, const int*, get_conf_signature) \
C_CODE(generatorType, bool, advanceToNextConfiguration) \
DELETE(generatorType)



//______________________________________________________THRESHOLD GENERATOR
void* setupIsoThresholdGenerator(int dimNumber,
                                 const int* isotopeNumbers,
                                 const int* atomCounts,
                                 const double* isotopeMasses,
                                 const double* isotopeProbabilities,
                                 const double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize)
{
    const double** IM = new const double*[dimNumber];
    const double** IP = new const double*[dimNumber];
    int idx = 0;
    for(int i=0; i<dimNumber; i++)
    {
        IM[i] = &isotopeMasses[idx];
        IP[i] = &isotopeProbabilities[idx];
        idx += isotopeNumbers[i];
    }
    //TODO in place (maybe pass a numpy matrix??)

    IsoThresholdGenerator* iso = new IsoThresholdGenerator(
        Iso(dimNumber, isotopeNumbers, atomCounts, IM, IP),
        threshold,
        _absolute,
        _tabSize,
        _hashSize);

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}
C_CODES(IsoThresholdGenerator)

//______________________________________________________ORDERED GENERATOR
void* setupIsoOrderedGenerator(int dimNumber,
                               const int* isotopeNumbers,
                               const int* atomCounts,
                               const double* isotopeMasses,
                               const double* isotopeProbabilities,
                               int _tabSize,
                               int _hashSize)
{
    const double** IM = new const double*[dimNumber];
    const double** IP = new const double*[dimNumber];
    int idx = 0;
    for(int i=0; i<dimNumber; i++)
    {
        IM[i] = &isotopeMasses[idx];
        IP[i] = &isotopeProbabilities[idx];
        idx += isotopeNumbers[i];
    }
    //TODO in place (maybe pass a numpy matrix??)

    IsoOrderedGenerator* iso = new IsoOrderedGenerator(
        Iso(dimNumber, isotopeNumbers, atomCounts, IM, IP),
        _tabSize,
        _hashSize);

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}
C_CODES(IsoOrderedGenerator)

//______________________________________________________ Threshold Tabulator 1.0

void* setupTabulator(void* generator,
                     bool  get_masses,
                     bool  get_probs,
                     bool  get_lprobs,
                     bool  get_confs)
{
    Tabulator* tabulator = new Tabulator(reinterpret_cast<IsoThresholdGenerator*>(generator),
                                         get_masses,
                                         get_probs,
                                         get_lprobs,
                                         get_confs);

    return reinterpret_cast<void*>(tabulator);
}

DELETE(Tabulator)
C_CODE(Tabulator, double*, masses)
C_CODE(Tabulator, double*, lprobs)
C_CODE(Tabulator, double*, probs)
C_CODE(Tabulator, int*,    confs)
C_CODE(Tabulator, int,     confs_no)

}  //extern "C" ends here
