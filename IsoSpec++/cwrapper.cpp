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


#include <cstring>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include "cwrapper.h"
#include "misc.h"
#include "marginalTrek++.h"
#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include "fasta.h"

using namespace IsoSpec;  // NOLINT(build/namespaces) - all of this really should be in a namespace IsoSpec, but C doesn't have them...


extern "C"
{
void * setupIso(int             dimNumber,
                const int*      isotopeNumbers,
                const int*      atomCounts,
                const double*   isotopeMasses,
                const double*   isotopeProbabilities)
{
    Iso* iso = new Iso(dimNumber, isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities);

    return reinterpret_cast<void*>(iso);
}

void * isoFromFasta(const char* fasta, bool use_nominal_masses, bool add_water)
{
    Iso* iso = new Iso(Iso::FromFASTA(fasta, use_nominal_masses, add_water));

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
    if(ret != nullptr)
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



// ______________________________________________________THRESHOLD GENERATOR
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


// ______________________________________________________LAYERED GENERATOR
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


// ______________________________________________________ORDERED GENERATOR
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

// ______________________________________________________STOCHASTIC GENERATOR
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

// ______________________________________________________ FixedEnvelopes

void* setupThresholdFixedEnvelope(void* iso,
                     double threshold,
                     bool absolute,
                     bool  get_confs)
{
    FixedEnvelope* ret = new FixedEnvelope(  // Use copy elision to allocate on heap with named constructor
            FixedEnvelope::FromThreshold(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                         threshold,
                                         absolute,
                                         get_confs));

    return reinterpret_cast<void*>(ret);
}

void* setupTotalProbFixedEnvelope(void* iso,
                     double target_coverage,
                     bool optimize,
                     bool get_confs)
{
    FixedEnvelope* ret = new FixedEnvelope(  // Use copy elision to allocate on heap with named constructor
            FixedEnvelope::FromTotalProb(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                         target_coverage,
                                         optimize,
                                         get_confs));

    return reinterpret_cast<void*>(ret);
}

void* setupStochasticFixedEnvelope(void* iso,
                    size_t no_molecules,
                    double precision,
                    double beta_bias,
                    bool get_confs)
{
    FixedEnvelope* ret = new FixedEnvelope(  // Use copy elision to allocate on heap with named constructor
            FixedEnvelope::FromStochastic(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                          no_molecules,
                                          precision,
                                          beta_bias,
                                          get_confs));

    return reinterpret_cast<void*>(ret);
}


void* setupBinnedFixedEnvelope(void* iso,
                    double target_total_prob,
                    double bin_width,
                    double bin_middle)
{
    FixedEnvelope* ret = new FixedEnvelope(  // Use copy elision to allocate on heap with named constructor
            FixedEnvelope::Binned(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                  target_total_prob,
                                  bin_width,
                                  bin_middle));

    return reinterpret_cast<void*>(ret);
}

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
        tt->release_masses();
        tt->release_probs();
        tt->release_confs();
    }
    delete tt;
}

void* copyFixedEnvelope(void* other)
{
    FixedEnvelope* ret = new FixedEnvelope(*reinterpret_cast<FixedEnvelope*>(other));
    return reinterpret_cast<void*>(ret);
}

const double* massesFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->release_masses();
}

const double* probsFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->release_probs();
}

const int* confsFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->release_confs();
}

size_t confs_noFixedEnvelope(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->confs_no();
}

double empiricAverageMass(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->empiric_average_mass();
}

double empiricVariance(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->empiric_variance();
}

double empiricStddev(void* tabulator)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator)->empiric_stddev();
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

double orientedWassersteinDistance(void* tabulator1, void* tabulator2)
{
    try
    {
        return reinterpret_cast<FixedEnvelope*>(tabulator1)->OrientedWassersteinDistance(*reinterpret_cast<FixedEnvelope*>(tabulator2));
    }
    catch(std::logic_error&)
    {
        return NAN;
    }
}

double abyssalWassersteinDistance(void* tabulator1, void* tabulator2, double abyss_depth, double other_scale)
{
    return reinterpret_cast<FixedEnvelope*>(tabulator1)->AbyssalWassersteinDistance(*reinterpret_cast<FixedEnvelope*>(tabulator2), abyss_depth, other_scale);
}

#if 0
double abyssalWassersteinDistanceGrad(void* const* envelopes, const double* scales, double* ret_gradient, size_t N, double abyss_depth_exp, double abyss_depth_the)
{
    return AbyssalWassersteinDistanceGrad(reinterpret_cast<FixedEnvelope* const*>(envelopes), scales, ret_gradient, N, abyss_depth_exp, abyss_depth_the);
}
#endif

struct ws_match_res wassersteinMatch(void* tabulator1, void* tabulator2, double flow_dist, double other_scale)
{
    struct ws_match_res res;
    auto tuple = reinterpret_cast<FixedEnvelope*>(tabulator1)->WassersteinMatch(*reinterpret_cast<FixedEnvelope*>(tabulator2), flow_dist, other_scale);
    res.res1 = std::get<0>(tuple);
    res.res2 = std::get<1>(tuple);
    res.flow = std::get<2>(tuple);
    return res;
}

void* addEnvelopes(void* tabulator1, void* tabulator2)
{
    //  Hopefully the compiler will do the copy elision...
    return reinterpret_cast<void*>(new FixedEnvelope((*reinterpret_cast<FixedEnvelope*>(tabulator1))+(*reinterpret_cast<FixedEnvelope*>(tabulator2))));
}

void* convolveEnvelopes(void* tabulator1, void* tabulator2)
{
    //  Hopefully the compiler will do the copy elision...
    return reinterpret_cast<void*>(new FixedEnvelope((*reinterpret_cast<FixedEnvelope*>(tabulator1))*(*reinterpret_cast<FixedEnvelope*>(tabulator2))));
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

void shiftMassEnvelope(void* envelope, double d_mass)
{
    reinterpret_cast<FixedEnvelope*>(envelope)->shift_mass(d_mass);
}

void resampleEnvelope(void* envelope, size_t ionic_current, double beta_bias)
{
    reinterpret_cast<FixedEnvelope*>(envelope)->resample(ionic_current, beta_bias);
}


void* binnedEnvelope(void* envelope, double width, double middle)
{
    //  Again, counting on copy elision...
    return reinterpret_cast<void*>(new FixedEnvelope(reinterpret_cast<FixedEnvelope*>(envelope)->bin(width, middle)));
}

void* linearCombination(void* const * const envelopes, const double* intensities, size_t count)
{
    //  Same...
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

void array_add(double* array, size_t N, double what)
{
    for(size_t ii = 0; ii < N; ii++)
        array[ii] += what;
}

void array_mul(double* array, size_t N, double what)
{
    for(size_t ii = 0; ii < N; ii++)
        array[ii] *= what;
}

void array_fma(double* array, size_t N, double mul, double add)
{
    for(size_t ii = 0; ii < N; ii++)
#if defined(FP_FAST_FMA)
        array[ii] = std::fma(array[ii], mul, add);
#else
        array[ii] += (array[ii] * mul) + add;
#endif
}

void parse_fasta_c(const char* fasta, int atomCounts[6])
{
    // Same thing, only this time with C linkage
    parse_fasta(fasta, atomCounts);
}
}  //  extern "C" ends here
