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


// ATTENTION! BELOW THIS LINE MATTEO WAS CODING AND IT IS BETTER NOT TO COMPILE THAT
#define C_CODE(generatorType, dataType, method)\
dataType method##generatorType(void* generator){ return reinterpret_cast<generatorType*>(generator)->method(); }

#define C_CODE_GET_CONF_SIGNATURE(generatorType)\
void get_conf_signature##generatorType(void* generator, int* space)\
{ reinterpret_cast<generatorType*>(generator)->get_conf_signature(space); }


#define DELETE(generatorType) void delete##generatorType(void* generator){ delete reinterpret_cast<generatorType*>(generator); }

#define C_CODES(generatorType)\
C_CODE(generatorType, double, mass) \
C_CODE(generatorType, double, lprob) \
C_CODE_GET_CONF_SIGNATURE(generatorType) \
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


//______________________________________________________LAYERED GENERATOR
void* setupIsoLayeredGenerator(int dimNumber,
                               const int* isotopeNumbers,
                                 const int* atomCounts,
                                 const double* isotopeMasses,
                                 const double* isotopeProbabilities,
                                 double _delta,
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

    IsoThresholdGenerator* iso = new IsoLayeredGenerator(
        Iso(dimNumber, isotopeNumbers, atomCounts, IM, IP),
        _delta,
        _tabSize,
        _hashSize);

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}
C_CODES(IsoLayeredGenerator)


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
