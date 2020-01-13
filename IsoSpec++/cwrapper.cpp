/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
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

double getIsoVariance(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->variance();
}

double getIsoStddev(void* iso)
{
    return reinterpret_cast<Iso*>(iso)->stddev();
}


double* getMarginalLogSizeEstimates(void* iso, double target_total_prob)
{
    Iso* i = reinterpret_cast<Iso*>(iso);
    double* ret = reinterpret_cast<double*>(malloc(sizeof(double)*i->getDimNumber()));
    i->saveMarginalLogSizeEstimates(ret, target_total_prob);
    return ret;
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
                                 int _hashSize,
                                 bool reorder_marginals)
{
    IsoThresholdGenerator* iso_tmp = new IsoThresholdGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        threshold,
        _absolute,
        _tabSize,
        _hashSize,
        reorder_marginals);

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoThresholdGenerator)


//______________________________________________________LAYERED GENERATOR
void* setupIsoLayeredGenerator(void* iso,
                     int _tabSize,
                     int _hashSize,
                     bool reorder_marginals,
                     double t_prob_hint
                )
{
    IsoLayeredGenerator* iso_tmp = new IsoLayeredGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        _tabSize,
        _hashSize,
        reorder_marginals,
        t_prob_hint
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

//______________________________________________________STOCHASTIC GENERATOR
void* setupIsoStochasticGenerator(void* iso,
                                  size_t no_molecules,
                                  double precision,
                                  double beta_bias)
{
    IsoStochasticGenerator* iso_tmp = new IsoStochasticGenerator(
        std::move(*reinterpret_cast<Iso*>(iso)),
        no_molecules,
        precision,
        beta_bias);

    return reinterpret_cast<void*>(iso_tmp);
}
ISOSPEC_C_FN_CODES(IsoStochasticGenerator)

//______________________________________________________ Threshold FixedEnvelope

void* setupThresholdFixedEnvelope(void* iso,
                     double threshold,
                     bool absolute,
                     bool  get_confs,
                     bool  get_lprobs,
                     bool  get_masses,
                     bool  get_probs)
{
    ThresholdFixedEnvelope* tabulator = new ThresholdFixedEnvelope(Iso(*reinterpret_cast<const Iso*>(iso), true),
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
    TotalProbFixedEnvelope* tabulator = new TotalProbFixedEnvelope(Iso(*reinterpret_cast<const Iso*>(iso), true),
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

//______________________________________________________ Generic FixedEnvelope

void* setupFixedEnvelope(double* masses, double* probs, size_t size, bool mass_sorted, bool prob_sorted, double total_prob)
{
    FixedEnvelope* ret = new FixedEnvelope(masses, probs, size, mass_sorted, prob_sorted, total_prob);
    return reinterpret_cast<void*>(ret);
}

void deleteFixedEnvelope(void* t, bool release_everything)
{
    FixedEnvelope* tt = reinterpret_cast<FixedEnvelope*>(t);
    if(release_everything)
    {
        tt->release_lprobs();
        tt->release_masses();
        tt->release_probs();
        tt->release_confs();
    }
    delete tt;
}

const double* massesFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->release_masses();
}

const double* lprobsFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->release_lprobs();
}

const double* probsFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->release_probs();
}

const int* confsFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->release_confs();
}

int confs_noFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->confs_no();
}

double wassersteinDistance(void* tabulator1, void* tabulator2)
{
    try
    {
        return reinterpret_cast<FixedEnvelope*>(tabulator1)->WassersteinDistance(*reinterpret_cast<FixedEnvelope*>(tabulator2));
    }
    catch(std::logic_error&)
    {
        return NAN;
    }
}

void* addEnvelopes(void* tabulator1, void* tabulator2)
{
    // Hopefully the compiler will do the copy elision...
    return reinterpret_cast<void*>(new FixedEnvelope(*reinterpret_cast<FixedEnvelope*>(tabulator1)+*reinterpret_cast<FixedEnvelope*>(tabulator2)));
}

void* convolveEnvelopes(void* tabulator1, void* tabulator2)
{
    // Hopefully the compiler will do the copy elision...
    return reinterpret_cast<void*>(new FixedEnvelope(*reinterpret_cast<FixedEnvelope*>(tabulator1)**reinterpret_cast<FixedEnvelope*>(tabulator2)));
}

double getTotalProbOfEnvelope(void* envelope)
{
    return reinterpret_cast<FixedEnvelope*>(envelope)->get_total_prob();
}

void scaleEnvelope(void* envelope, double factor)
{
    reinterpret_cast<FixedEnvelope*>(envelope)->scale(factor);
}

void normalizeEnvelope(void* envelope)
{
    reinterpret_cast<FixedEnvelope*>(envelope)->normalize();
}

void* binnedEnvelope(void* envelope, double width, double middle)
{
    // Again, counting on copy elision...
    return reinterpret_cast<void*>(new FixedEnvelope(reinterpret_cast<FixedEnvelope*>(envelope)->bin(width, middle)));
}

void* linearCombination(void* const * const envelopes, const double* intensities, size_t count)
{
    // Same...
    return reinterpret_cast<void*>(new FixedEnvelope(FixedEnvelope::LinearCombination(reinterpret_cast<const FixedEnvelope* const *>(envelopes), intensities, count)));
}

void sortEnvelopeByMass(void* envelope)
{
    reinterpret_cast<FixedEnvelope*>(envelope)->sort_by_mass();
}

void sortEnvelopeByProb(void* envelope)
{
    reinterpret_cast<FixedEnvelope*>(envelope)->sort_by_prob();
}




void freeReleasedArray(void* array)
{
    free(array);
}

}  //extern "C" ends here
