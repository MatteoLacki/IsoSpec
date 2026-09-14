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
#include <type_traits>
#include <limits>
#include <cmath>
#include <memory>
#include <vector>
#include "cwrapper.h"
#include "misc.h"
#include "marginalTrek++.h"
#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include "fasta.h"
#include "fasta_mods.h"
#include "element_tables.h"

using namespace IsoSpec;  // NOLINT(build/namespaces) - all of this really should be in a namespace IsoSpec, but C doesn't have them...


// Letting a C++ exception unwind across an extern "C" function into a C (or, in
// our case, Python/cffi) caller is undefined behaviour. The library internals
// throw std::invalid_argument / std::length_error / std::bad_alloc on bad input
// or OOM, so every wrapper below funnels its body through c_guard(), which turns
// any escaped exception into an out-of-band error value:
//   * pointer returns  -> nullptr   (callers must null-check the returned handle)
//   * floating returns -> NaN       (matches the pre-existing Wasserstein behaviour)
//   * integral/bool    -> zero / false
//   * void             -> no-op
// See the matching note in cwrapper.h.
namespace
{
template<typename T> struct c_fallback         { static T      value() noexcept { return T(); } };
template<>           struct c_fallback<double>  { static double value() noexcept { return std::numeric_limits<double>::quiet_NaN(); } };
template<>           struct c_fallback<float>   { static float  value() noexcept { return std::numeric_limits<float>::quiet_NaN(); } };

template<typename F>
static inline auto c_guard(F&& f) noexcept -> decltype(f())
{
    using R = decltype(f());
    try
    {
        return f();
    }
    catch(...)
    {
        if constexpr (std::is_void_v<R>)
            return;
        else
            return c_fallback<R>::value();
    }
}
}  // anonymous namespace


extern "C"
{
void * setupIso(int             dimNumber,
                const int*      isotopeNumbers,
                const int*      atomCounts,
                const double*   isotopeMasses,
                const double*   isotopeProbabilities)
{
    return c_guard([&]() -> void*
    {
        return new Iso(dimNumber, isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities);
    });
}

void * isoFromFasta(const char* fasta, bool use_nominal_masses, bool add_water)
{
    return c_guard([&]() -> void*
    {
        return new Iso(Iso::FromFASTA(fasta, use_nominal_masses, add_water));
    });
}

void deleteIso(void* iso)
{
    delete reinterpret_cast<Iso*>(iso);
}

double getLightestPeakMassIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getLightestPeakMass(); });
}

double getLightestPeakLProbIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getLightestPeakLProb(); });
}

void getLightestPeakSignature(void* iso, int* space)
{
    return c_guard([&]{ reinterpret_cast<Iso*>(iso)->getLightestPeakSignature(space); });
}

double getHeaviestPeakMassIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getHeaviestPeakMass(); });
}

double getHeaviestPeakLProbIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getHeaviestPeakLProb(); });
}

void getHeaviestPeakSignature(void* iso, int* space)
{
    return c_guard([&]{ reinterpret_cast<Iso*>(iso)->getHeaviestPeakSignature(space); });
}

double getMonoisotopicPeakMassIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getMonoisotopicPeakMass(); });
}

double getMonoisotopicPeakLProbIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getMonoisotopicPeakLProb(); });
}

void getMonoisotopicPeakSignature(void* iso, int* space)
{
    return c_guard([&]{ reinterpret_cast<Iso*>(iso)->getMonoisotopicPeakSignature(space); });
}

double getModeLProbIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getModeLProb(); });
}

double getModeMassIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getModeMass(); });
}

double getTheoreticalAverageMassIso(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->getTheoreticalAverageMass(); });
}

double getIsoVariance(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->variance(); });
}

double getIsoStddev(void* iso)
{
    return c_guard([&]{ return reinterpret_cast<Iso*>(iso)->stddev(); });
}


double* getMarginalLogSizeEstimates(void* iso, double target_total_prob)
{
    return c_guard([&]() -> double*
    {
        Iso* i = reinterpret_cast<Iso*>(iso);
        double* ret = reinterpret_cast<double*>(malloc(sizeof(double)*i->getDimNumber()));
        if(ret == nullptr)
            throw std::bad_alloc();
        i->saveMarginalLogSizeEstimates(ret, target_total_prob);
        return ret;
    });
}



#define ISOSPEC_C_FN_CODE(generatorType, dataType, method)\
dataType method##generatorType(void* generator){ return c_guard([&]{ return reinterpret_cast<generatorType*>(generator)->method(); }); }

#define ISOSPEC_C_FN_CODE_GET_CONF_SIGNATURE(generatorType)\
void get_conf_signature##generatorType(void* generator, int* space)\
{ return c_guard([&]{ reinterpret_cast<generatorType*>(generator)->get_conf_signature(space); }); }


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
    return c_guard([&]() -> void*
    {
        return new IsoThresholdGenerator(
            std::move(*reinterpret_cast<Iso*>(iso)),
            threshold,
            _absolute,
            _tabSize,
            _hashSize,
            reorder_marginals);
    });
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
    return c_guard([&]() -> void*
    {
        return new IsoLayeredGenerator(
            std::move(*reinterpret_cast<Iso*>(iso)),
            _tabSize,
            _hashSize,
            reorder_marginals,
            t_prob_hint);
    });
}
ISOSPEC_C_FN_CODES(IsoLayeredGenerator)


// ______________________________________________________ORDERED GENERATOR
void* setupIsoOrderedGenerator(void* iso,
                               int _tabSize,
                               int _hashSize)
{
    return c_guard([&]() -> void*
    {
        return new IsoOrderedGenerator(
            std::move(*reinterpret_cast<Iso*>(iso)),
            _tabSize,
            _hashSize);
    });
}
ISOSPEC_C_FN_CODES(IsoOrderedGenerator)

// ______________________________________________________STOCHASTIC GENERATOR
void* setupIsoStochasticGenerator(void* iso,
                                  size_t no_molecules,
                                  double precision,
                                  double beta_bias)
{
    return c_guard([&]() -> void*
    {
        return new IsoStochasticGenerator(
            std::move(*reinterpret_cast<Iso*>(iso)),
            no_molecules,
            precision,
            beta_bias);
    });
}
ISOSPEC_C_FN_CODES(IsoStochasticGenerator)

// ______________________________________________________ FixedEnvelopes

void* setupThresholdFixedEnvelope(void* iso,
                     double threshold,
                     bool absolute,
                     bool  get_confs)
{
    return c_guard([&]() -> void*
    {  // Use copy elision to allocate on heap with named constructor
        return new FixedEnvelope(
            FixedEnvelope::FromThreshold(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                         threshold,
                                         absolute,
                                         get_confs));
    });
}

void* setupTotalProbFixedEnvelope(void* iso,
                     double target_coverage,
                     bool optimize,
                     bool get_confs)
{
    return c_guard([&]() -> void*
    {  // Use copy elision to allocate on heap with named constructor
        return new FixedEnvelope(
            FixedEnvelope::FromTotalProb(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                         target_coverage,
                                         optimize,
                                         get_confs));
    });
}

void* setupStochasticFixedEnvelope(void* iso,
                    size_t no_molecules,
                    double precision,
                    double beta_bias,
                    bool get_confs)
{
    return c_guard([&]() -> void*
    {  // Use copy elision to allocate on heap with named constructor
        return new FixedEnvelope(
            FixedEnvelope::FromStochastic(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                          no_molecules,
                                          precision,
                                          beta_bias,
                                          get_confs));
    });
}


void* setupBinnedFixedEnvelope(void* iso,
                    double target_total_prob,
                    double bin_width,
                    double bin_middle)
{
    return c_guard([&]() -> void*
    {  // Use copy elision to allocate on heap with named constructor
        return new FixedEnvelope(
            FixedEnvelope::Binned(Iso(*reinterpret_cast<const Iso*>(iso), true),
                                  target_total_prob,
                                  bin_width,
                                  bin_middle));
    });
}

void* setupFixedEnvelope(double* masses, double* probs, size_t size, bool mass_sorted, bool prob_sorted, double total_prob)
{
    return c_guard([&]() -> void*
    {
        return new FixedEnvelope(masses, probs, size, mass_sorted, prob_sorted, total_prob);
    });
}

void* setupFixedEnvelopeWithConfs(double* masses, double* probs, int* confs, size_t size, int allDim, bool mass_sorted, bool prob_sorted, double total_prob)
{
    return c_guard([&]() -> void*
    {
        return new FixedEnvelope(masses, probs, confs, size, allDim, mass_sorted, prob_sorted, total_prob);
    });
}

void deleteFixedEnvelope(void* t, bool release_everything)
{
    FixedEnvelope* tt = reinterpret_cast<FixedEnvelope*>(t);
    if(release_everything)
        // Drops the arrays without freeing (and without the copy release_masses()
        // and friends would make of a non-malloc()'d one): this is the path used
        // for envelopes wrapped around buffers their real owner still holds.
        tt->release_everything();
    delete tt;
}

void* copyFixedEnvelope(void* other)
{
    return c_guard([&]() -> void*
    {
        return new FixedEnvelope(*reinterpret_cast<FixedEnvelope*>(other));
    });
}

const double* massesFixedEnvelope(void* tabulator)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator)->release_masses(); });
}

const double* probsFixedEnvelope(void* tabulator)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator)->release_probs(); });
}

const int* confsFixedEnvelope(void* tabulator)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator)->release_confs(); });
}

// Shared tail of the three *WithDeleter wrappers: publish the (size, deleter)
// pair through the out-parameters and return the pointer. Not a template, and
// not in an anonymous namespace, because everything here has C linkage.
static void* publish_released(const ReleasedArray& r, size_t* size_out, IsoSpecArrayDeleter* deleter_out)
{
    if(size_out != nullptr)
        *size_out = r.size;
    if(deleter_out != nullptr)
        *deleter_out = r.deleter;
    return r.ptr;
}

double* massesFixedEnvelopeWithDeleter(void* tabulator, size_t* size_out, IsoSpecArrayDeleter* deleter_out)
{
    return c_guard([&]{ return reinterpret_cast<double*>(publish_released(
        reinterpret_cast<FixedEnvelope*>(tabulator)->release_masses_with_deleter(), size_out, deleter_out)); });
}

double* probsFixedEnvelopeWithDeleter(void* tabulator, size_t* size_out, IsoSpecArrayDeleter* deleter_out)
{
    return c_guard([&]{ return reinterpret_cast<double*>(publish_released(
        reinterpret_cast<FixedEnvelope*>(tabulator)->release_probs_with_deleter(), size_out, deleter_out)); });
}

int* confsFixedEnvelopeWithDeleter(void* tabulator, size_t* size_out, IsoSpecArrayDeleter* deleter_out)
{
    return c_guard([&]{ return reinterpret_cast<int*>(publish_released(
        reinterpret_cast<FixedEnvelope*>(tabulator)->release_confs_with_deleter(), size_out, deleter_out)); });
}

void freeReleasedArrayWithDeleter(void* array, size_t size, IsoSpecArrayDeleter deleter)
{
    if(array != nullptr && deleter != nullptr)
        deleter(array, size);
}

size_t confs_noFixedEnvelope(void* tabulator)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator)->confs_no(); });
}

double empiricAverageMass(void* tabulator)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator)->empiric_average_mass(); });
}

double empiricVariance(void* tabulator)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator)->empiric_variance(); });
}

double empiricStddev(void* tabulator)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator)->empiric_stddev(); });
}

double wassersteinDistance(void* tabulator1, void* tabulator2)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator1)->WassersteinDistance(*reinterpret_cast<FixedEnvelope*>(tabulator2)); });
}

double orientedWassersteinDistance(void* tabulator1, void* tabulator2)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator1)->OrientedWassersteinDistance(*reinterpret_cast<FixedEnvelope*>(tabulator2)); });
}

double abyssalWassersteinDistance(void* tabulator1, void* tabulator2, double abyss_depth, double other_scale)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(tabulator1)->AbyssalWassersteinDistance(*reinterpret_cast<FixedEnvelope*>(tabulator2), abyss_depth, other_scale); });
}

struct ws_match_res wassersteinMatch(void* tabulator1, void* tabulator2, double flow_dist, double other_scale)
{
    return c_guard([&]() -> struct ws_match_res
    {
        struct ws_match_res res;
        auto tuple = reinterpret_cast<FixedEnvelope*>(tabulator1)->WassersteinMatch(*reinterpret_cast<FixedEnvelope*>(tabulator2), flow_dist, other_scale);
        res.res1 = std::get<0>(tuple);
        res.res2 = std::get<1>(tuple);
        res.flow = std::get<2>(tuple);
        return res;
    });
}

void* addEnvelopes(void* tabulator1, void* tabulator2)
{
    return c_guard([&]() -> void*
    {  //  Hopefully the compiler will do the copy elision...
        return new FixedEnvelope((*reinterpret_cast<FixedEnvelope*>(tabulator1))+(*reinterpret_cast<FixedEnvelope*>(tabulator2)));
    });
}

void* convolveEnvelopes(void* tabulator1, void* tabulator2)
{
    return c_guard([&]() -> void*
    {  //  Hopefully the compiler will do the copy elision...
        return new FixedEnvelope((*reinterpret_cast<FixedEnvelope*>(tabulator1))*(*reinterpret_cast<FixedEnvelope*>(tabulator2)));
    });
}

double getTotalProbOfEnvelope(void* envelope)
{
    return c_guard([&]{ return reinterpret_cast<FixedEnvelope*>(envelope)->get_total_prob(); });
}

void scaleEnvelope(void* envelope, double factor)
{
    return c_guard([&]{ reinterpret_cast<FixedEnvelope*>(envelope)->scale(factor); });
}

void normalizeEnvelope(void* envelope)
{
    return c_guard([&]{ reinterpret_cast<FixedEnvelope*>(envelope)->normalize(); });
}

void shiftMassEnvelope(void* envelope, double d_mass)
{
    return c_guard([&]{ reinterpret_cast<FixedEnvelope*>(envelope)->shift_mass(d_mass); });
}

void resampleEnvelope(void* envelope, size_t ionic_current, double beta_bias)
{
    return c_guard([&]{ reinterpret_cast<FixedEnvelope*>(envelope)->resample(ionic_current, beta_bias); });
}


void* binnedEnvelope(void* envelope, double width, double middle)
{
    return c_guard([&]() -> void*
    {  //  Again, counting on copy elision...
        return new FixedEnvelope(reinterpret_cast<FixedEnvelope*>(envelope)->bin(width, middle));
    });
}

void* linearCombination(void* const * const envelopes, const double* intensities, size_t count)
{
    return c_guard([&]() -> void*
    {  //  Same...
        return new FixedEnvelope(FixedEnvelope::LinearCombination(reinterpret_cast<const FixedEnvelope* const *>(envelopes), intensities, count));
    });
}

void sortEnvelopeByMass(void* envelope)
{
    return c_guard([&]{ reinterpret_cast<FixedEnvelope*>(envelope)->sort_by_mass(); });
}

void sortEnvelopeByProb(void* envelope)
{
    return c_guard([&]{ reinterpret_cast<FixedEnvelope*>(envelope)->sort_by_prob(); });
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
        array[ii] = (array[ii] * mul) + add;
#endif
}

void parse_fasta_c(const char* fasta, int atomCounts[6])
{
    // Same thing, only this time with C linkage
    parse_fasta(fasta, atomCounts);
}

namespace
{
// Owns a parsed composition's storage across the C ABI: the symbol pointers
// point into element_tables.cpp's static elem_table_symbol strings (program
// lifetime, never freed here), only the counts and the pointer array itself
// are this handle's own allocation.
struct ElementCompositionHandle
{
    std::vector<const char*> symbols;
    std::vector<int> counts;
};

const UnimodTable& resolve_unimod_table(const char* unimod_db_path)
{
    if(unimod_db_path == nullptr || unimod_db_path[0] == '\0')
        return embedded_unimod_table();
    return unimod_table_for_path(unimod_db_path);
}
}  // anonymous namespace

void* isoFromFastaWithMods(const char* sequence, bool use_nominal_masses, bool add_water, const char* unimod_db_path)
{
    return c_guard([&]() -> void*
    {
        return new Iso(Iso::FromFASTAWithMods(sequence, use_nominal_masses, add_water, unimod_db_path));
    });
}

void* parseFastaWithModsC(const char* sequence, const char* unimod_db_path)
{
    return c_guard([&]() -> void*
    {
        // Reused across calls (parse_fasta_with_mods_into only clears it,
        // never releases capacity) -- avoids reallocating the scratch
        // ElementComposition on every single call from Python/R, which is
        // the real hot path this matters for (millions of peptides). Plain
        // `static`, not `thread_local`: IsoSpec is single-threaded by design
        // (this repo's CLAUDE.md) so they'd behave identically here, and
        // thread_local's extra TLS teardown machinery is a real, observed
        // crash risk for a dynamically-loaded Python/R extension module
        // (segfault during interpreter shutdown, seen in CI on macOS) that
        // buys nothing in a library that never has a second thread. Nothing
        // below retains a pointer into `scratch` itself: `counts` is copied
        // out by value, and `symbols` is built from elem_table_symbol's
        // static string-literal pointers, not from `scratch`'s own storage --
        // safe to overwrite `scratch` on the next call.
        static ElementComposition scratch;
        parse_fasta_with_mods_into(sequence, scratch, resolve_unimod_table(unimod_db_path));
        std::unique_ptr<ElementCompositionHandle> handle(new ElementCompositionHandle());
        handle->counts = scratch.count;
        handle->symbols.reserve(scratch.element_first_index.size());
        for(int idx : scratch.element_first_index)
            handle->symbols.push_back(elem_table_symbol[idx]);
        return handle.release();
    });
}

size_t compositionSizeC(void* composition)
{
    return c_guard([&]() -> size_t
    {
        return reinterpret_cast<ElementCompositionHandle*>(composition)->counts.size();
    });
}

const char* const* compositionSymbolsC(void* composition)
{
    return c_guard([&]() -> const char* const*
    {
        return reinterpret_cast<ElementCompositionHandle*>(composition)->symbols.data();
    });
}

const int* compositionCountsC(void* composition)
{
    return c_guard([&]() -> const int*
    {
        return reinterpret_cast<ElementCompositionHandle*>(composition)->counts.data();
    });
}

void deleteCompositionC(void* composition)
{
    delete reinterpret_cast<ElementCompositionHandle*>(composition);
}

const char* activeSimdLevel()
{
    return active_simd_level();
}
}  //  extern "C" ends here
