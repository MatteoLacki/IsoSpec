/*
 *   Copyright (C) 2015-2019 Mateusz Łącki and Michał Startek.
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
#include "fixedEnvelopes.h"

using namespace IsoSpec;


extern "C"
{
void * setupIso(int             dimNumber,
                const int*      isotopeNumbers,
                const int*      atomCounts,
                const double*   isotopeMasses,
                const double*   isotopeProbabilities)
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

    Iso* iso = new Iso(dimNumber, isotopeNumbers, atomCounts, IM, IP);

    delete[] IM;
    delete[] IP;

    return reinterpret_cast<void*>(iso);
}

void deleteIso(void* iso)
{
    delete reinterpret_cast<Iso*>(iso);
}

double getLightestPeakMassIso(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->getLightestPeakMass();
}

double getHeaviestPeakMassIso(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->getHeaviestPeakMass();
}

double getMonoisotopicPeakMassIso(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->getMonoisotopicPeakMass();
}

double getModeLProbIso(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->getModeLProb();
}

double getModeMassIso(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->getModeMass();
}

double getTheoreticalAverageMassIso(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->getTheoreticalAverageMass();
}



#define ISOSPEC_C_FN_CODE(generatorType, dataType, method)\
dataType method##generatorType(void* generator){ return reinterpret_cast<generatorType*>(generator)->method(); }

#define ISOSPEC_C_FN_CODE_GET_CONF_SIGNATURE(generatorType)\
void get_conf_signature##generatorType(void* generator, int* space)\
{ reinterpret_cast<generatorType*>(generator)->get_conf_signature(space); }


#define ISOSPEC_C_FN_DELETE(generatorType) void delete##generatorType(void* generator){ delete reinterpret_cast<generatorType*>(generator); }

#define ISOSPEC_C_FN_CODES(generatorType)\
ISOSPEC_C_FN_CODE(generatorType, double, mass) \
ISOSPEC_C_FN_CODE(generatorType, double, lprob) \
ISOSPEC_C_FN_CODE(generatorType, double, prob) \
ISOSPEC_C_FN_CODE_GET_CONF_SIGNATURE(generatorType) \
ISOSPEC_C_FN_CODE(generatorType, bool, advanceToNextConfiguration) \
ISOSPEC_C_FN_DELETE(generatorType)



//______________________________________________________THRESHOLD GENERATOR
void* setupIsoThresholdGenerator(void* iso,
                                 double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize)
{
    IsoThresholdGenerator* iso_tmp = new IsoThresholdGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        threshold,
        _absolute,
        _tabSize,
        _hashSize);

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoThresholdGenerator)


//______________________________________________________LAYERED GENERATOR
void* setupIsoLayeredGenerator(void* iso,
                     int _tabSize,
                     int _hashSize
                )
{
    IsoLayeredGenerator* iso_tmp = new IsoLayeredGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        _tabSize,
        _hashSize
    );

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoLayeredGenerator)


//______________________________________________________ORDERED GENERATOR
void* setupIsoOrderedGenerator(void* iso,
                               int _tabSize,
                               int _hashSize)
{
    IsoOrderedGenerator* iso_tmp = new IsoOrderedGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        _tabSize,
        _hashSize);

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoOrderedGenerator)

//______________________________________________________ Threshold FixedEnvelope

void* setupThresholdFixedEnvelope(void* iso,
                     double threshold,
                     bool absolute,
                     bool  get_confs,
                     bool  get_lprobs,
                     bool  get_masses,
                     bool  get_probs)
{
    ThresholdFixedEnvelope* tabulator = new ThresholdFixedEnvelope(*reinterpret_cast<const Iso*>(iso),
                                         threshold,
                                         absolute,
                                         get_confs,
                                         get_lprobs,
                                         get_masses,
                                         get_probs);

    return reinterpret_cast<void*>(tabulator);
}

void deleteThresholdFixedEnvelope(void* t)
{
    delete reinterpret_cast<ThresholdFixedEnvelope*>(t);
}

const double* massesThresholdFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<ThresholdFixedEnvelope*>(tabulator)->release_masses();
}

const double* lprobsThresholdFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<ThresholdFixedEnvelope*>(tabulator)->release_lprobs();
}

const double* probsThresholdFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<ThresholdFixedEnvelope*>(tabulator)->release_probs();
}

const int*    confsThresholdFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<ThresholdFixedEnvelope*>(tabulator)->release_confs();
}

int confs_noThresholdFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<ThresholdFixedEnvelope*>(tabulator)->confs_no();
}


//______________________________________________________ Layered FixedEnvelope

void* setupTotalProbFixedEnvelope(void* iso,
                     double target_coverage,
                     bool optimize,
                     bool  get_confs,
                     bool  get_lprobs,
                     bool  get_masses,
                     bool  get_probs)
{
    TotalProbFixedEnvelope* tabulator = new TotalProbFixedEnvelope(*reinterpret_cast<const Iso*>(iso),
                                         target_coverage,
                                         optimize,
                                         get_confs,
                                         get_lprobs,
                                         get_masses,
                                         get_probs);

    return reinterpret_cast<void*>(tabulator);
}

void deleteTotalProbFixedEnvelope(void* t)
{
    delete reinterpret_cast<TotalProbFixedEnvelope*>(t);
}

const double* massesTotalProbFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<TotalProbFixedEnvelope*>(tabulator)->release_masses();
}

const double* lprobsTotalProbFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<TotalProbFixedEnvelope*>(tabulator)->release_lprobs();
}

const double* probsTotalProbFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<TotalProbFixedEnvelope*>(tabulator)->release_probs();
}

const int*    confsTotalProbFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<TotalProbFixedEnvelope*>(tabulator)->release_confs();
}

int confs_noTotalProbFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<TotalProbFixedEnvelope*>(tabulator)->confs_no();
}

void freeReleasedArray(void* array)
{
    free(array);
}
}  //extern "C" ends here
