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
#include "tabulator.h"

using namespace IsoSpec;

#define TABULATOR_LAYERED 1
#define TABULATOR_THRESHOLD 2

//dbl_param, bool_param, get_lprobs, get_masses, get_probs, get_confs)
#define MK_TABULATOR(get_lprobs, get_masses, get_probs, get_confs) \
{ \
if (tabulator_type == TABULATOR_LAYERED) \
    v_tabulator = new LayeredTabulator<get_lprobs, get_masses, get_probs, get_confs>(std::move(iso), dbl_param, bool_param); \
if (tabulator_type == TABULATOR_THRESHOLD) \
    v_tabulator = new ThresholdTabulator<get_lprobs, get_masses, get_probs, get_confs>(std::move(iso), dbl_param, bool_param); \
}


/* A trivial helper class for the various tabulators. 
class TabulatorShell : public TabulatorParentCls
{
public:
    TabulatorShell(Iso&& iso, int tabulator_type, double dbl_param, bool bool_param, bool get_lprobs, bool get_masses, bool get_probs, bool get_confs) :
    TabulatorParentCls()
    {
        TabulatorParentCls* v_tabulator = nullptr;

        if(get_lprobs)
        {
            if(get_masses)
            {
                if(get_probs)
                {
                    if(get_confs)
                        MK_TABULATOR(true, true, true, true)
                    else
                        MK_TABULATOR(true, true, true, false)
                }
                else
                {
                    if(get_confs)
                        MK_TABULATOR(true, true, false, true)
                    else
                        MK_TABULATOR(true, true, false, false)
                }
            }
            else
            {
                if(get_probs)
                {
                    if(get_confs)
                        MK_TABULATOR(true, false, true, true)
                    else
                        MK_TABULATOR(true, false, true, false)
                }
                else
                {
                    if(get_confs)
                        MK_TABULATOR(true, false, false, true)
                    else
                        MK_TABULATOR(true, false, false, false)
                }
            }
        }
        else
        {
            if(get_masses)
            {
                if(get_probs)
                {
                    if(get_confs)
                        MK_TABULATOR(false, true, true, true)
                    else
                        MK_TABULATOR(false, true, true, false)
                }
                else
                {
                    if(get_confs)
                        MK_TABULATOR(false, true, false, true)
                    else
                        MK_TABULATOR(false, true, false, false)
                }
            }
            else
            {
                if(get_probs)
                {
                    if(get_confs)
                        MK_TABULATOR(false, false, true, true)
                    else
                        MK_TABULATOR(false, false, true, false)
                }
                else
                {
                    if(get_confs)
                        MK_TABULATOR(false, false, false, true)
                    else
                        MK_TABULATOR(false, false, false, false)
                }
            }
        }


        _lprobs = get_lprobs ? v_tabulator->lprobs(true) : nullptr;
        _masses = get_masses ? v_tabulator->masses(true) : nullptr;
        _probs  = get_probs  ? v_tabulator->probs(true)  : nullptr;
        _confs  = get_confs  ? v_tabulator->confs(true)  : nullptr;

        _confs_no = v_tabulator->confs_no();
        allDim = v_tabulator->getAllDim();

        delete v_tabulator;
    }

    virtual ~TabulatorShell() {}
};
*/

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

//______________________________________________________ Threshold Tabulator

void* setupThresholdTabulator(void* iso,
                     double threshold,
                     bool absolute,
                     bool  get_masses,
                     bool  get_probs,
                     bool  get_lprobs,
                     bool  get_confs)
{
    ThresholdTabulator* tabulator = new ThresholdTabulator(std::move(*reinterpret_cast<Iso*>(iso)),
                                         threshold,
                                         absolute,
                                         get_lprobs,
                                         get_masses,
                                         get_probs,
                                         get_confs);

    return reinterpret_cast<void*>(tabulator);
}

void deleteThresholdTabulator(void* t)
{
    delete reinterpret_cast<ThresholdTabulator*>(t);
}

const double* massesThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<ThresholdTabulator*>(tabulator)->masses(true);
}

const double* lprobsThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<ThresholdTabulator*>(tabulator)->lprobs(true);
}

const double* probsThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<ThresholdTabulator*>(tabulator)->probs(true);
}

const int*    confsThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<ThresholdTabulator*>(tabulator)->confs(true);
}

int confs_noThresholdTabulator(void* tabulator)
{
    return reinterpret_cast<ThresholdTabulator*>(tabulator)->confs_no();
}


//______________________________________________________ Layered Tabulator

void* setupLayeredTabulator(void* iso,
                     bool  get_masses,
                     bool  get_probs,
                     bool  get_lprobs,
                     bool  get_confs,
                     double target_coverage,
                     bool optimize)
{
    LayeredTabulator* tabulator = new LayeredTabulator(std::move(*reinterpret_cast<Iso*>(iso)),
                                         target_coverage,
                                         optimize,
                                         get_lprobs,
                                         get_masses,
                                         get_probs,
                                         get_confs);

    return reinterpret_cast<void*>(tabulator);
}

void deleteLayeredTabulator(void* t)
{
    delete reinterpret_cast<LayeredTabulator*>(t);
}

const double* massesLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<LayeredTabulator*>(tabulator)->masses(true);
}

const double* lprobsLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<LayeredTabulator*>(tabulator)->lprobs(true);
}

const double* probsLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<LayeredTabulator*>(tabulator)->probs(true);
}

const int*    confsLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<LayeredTabulator*>(tabulator)->confs(true);
}

int confs_noLayeredTabulator(void* tabulator)
{
    return reinterpret_cast<LayeredTabulator*>(tabulator)->confs_no();
}

void freeReleasedArray(void* array)
{
    free(array);
}
}  //extern "C" ends here
