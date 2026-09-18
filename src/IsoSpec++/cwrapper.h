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

#pragma once

#include "platform.h"

#define ISOSPEC_ALGO_LAYERED 0
#define ISOSPEC_ALGO_ORDERED 1
#define ISOSPEC_ALGO_THRESHOLD_ABSOLUTE 2
#define ISOSPEC_ALGO_THRESHOLD_RELATIVE 3
#define ISOSPEC_ALGO_LAYERED_ESTIMATE 4


// Error reporting across the C ABI: these functions never let a C++ exception
// escape (that would be undefined behaviour for a C / cffi caller). On bad input
// or out-of-memory they instead return an out-of-band error value:
//   * functions returning a handle/pointer return NULL,
//   * functions returning a double return NaN,
//   * functions returning an integer/bool return 0 / false.
// Callers should check the returned handle for NULL before using it.

#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

ISOSPEC_C_API void * setupIso(int             dimNumber,
                const int*      isotopeNumbers,
                const int*      atomCounts,
                const double*   isotopeMasses,
                const double*   isotopeProbabilities);

ISOSPEC_C_API void * isoFromFasta(const char* fasta, bool use_nominal_masses, bool add_water);

ISOSPEC_C_API double getLightestPeakMassIso(void* iso);
ISOSPEC_C_API double getLightestPeakLProbIso(void* iso);
ISOSPEC_C_API void getLightestPeakSignature(void* iso, int* space);
ISOSPEC_C_API double getHeaviestPeakMassIso(void* iso);
ISOSPEC_C_API double getHeaviestPeakLProbIso(void* iso);
ISOSPEC_C_API void getHeaviestPeakSignature(void* iso, int* space);
ISOSPEC_C_API double getMonoisotopicPeakMassIso(void* iso);
ISOSPEC_C_API double getMonoisotopicPeakLProbIso(void* iso);
ISOSPEC_C_API void getMonoisotopicPeakSignature(void* iso, int* space);
ISOSPEC_C_API double getModeLProbIso(void* iso);
ISOSPEC_C_API double getModeMassIso(void* iso);
ISOSPEC_C_API double getTheoreticalAverageMassIso(void* iso);
ISOSPEC_C_API double getIsoVariance(void* iso);
ISOSPEC_C_API double getIsoStddev(void* iso);
ISOSPEC_C_API double* getMarginalLogSizeEstimates(void* iso, double target_total_prob);


ISOSPEC_C_API void deleteIso(void* iso);


// ______________________________________________________THRESHOLD GENERATOR
ISOSPEC_C_API void* setupIsoThresholdGenerator(void* iso,
                                 double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize,
                                 bool reorder_marginals);
ISOSPEC_C_API double massIsoThresholdGenerator(void* generator);
ISOSPEC_C_API double lprobIsoThresholdGenerator(void* generator);
ISOSPEC_C_API double probIsoThresholdGenerator(void* generator);
ISOSPEC_C_API void get_conf_signatureIsoThresholdGenerator(void* generator, int* space);
ISOSPEC_C_API bool advanceToNextConfigurationIsoThresholdGenerator(void* generator);
ISOSPEC_C_API void deleteIsoThresholdGenerator(void* generator);


// ______________________________________________________LAYERED GENERATOR
ISOSPEC_C_API void* setupIsoLayeredGenerator(void* iso,
                               int _tabSize,
                               int _hashSize,
                               bool reorder_marginals,
                               double t_prob_hint);
ISOSPEC_C_API double massIsoLayeredGenerator(void* generator);
ISOSPEC_C_API double lprobIsoLayeredGenerator(void* generator);
ISOSPEC_C_API double probIsoLayeredGenerator(void* generator);
ISOSPEC_C_API void get_conf_signatureIsoLayeredGenerator(void* generator, int* space);
ISOSPEC_C_API bool advanceToNextConfigurationIsoLayeredGenerator(void* generator);
ISOSPEC_C_API void deleteIsoLayeredGenerator(void* generator);

// ______________________________________________________ORDERED GENERATOR
ISOSPEC_C_API void* setupIsoOrderedGenerator(void* iso,
                               int _tabSize,
                               int _hashSize);
ISOSPEC_C_API double massIsoOrderedGenerator(void* generator);
ISOSPEC_C_API double lprobIsoOrderedGenerator(void* generator);
ISOSPEC_C_API double probIsoOrderedGenerator(void* generator);
ISOSPEC_C_API void get_conf_signatureIsoOrderedGenerator(void* generator, int* space);
ISOSPEC_C_API bool advanceToNextConfigurationIsoOrderedGenerator(void* generator);
ISOSPEC_C_API void deleteIsoOrderedGenerator(void* generator);

// ______________________________________________________STOCHASTIC GENERATOR
ISOSPEC_C_API void* setupIsoStochasticGenerator(void* iso,
                                   size_t no_molecules,
                                   double precision,
                                   double beta_bias);
ISOSPEC_C_API double massIsoStochasticGenerator(void* generator);
ISOSPEC_C_API double lprobIsoStochasticGenerator(void* generator);
ISOSPEC_C_API double probIsoStochasticGenerator(void* generator);
ISOSPEC_C_API void get_conf_signatureIsoStochasticGenerator(void* generator, int* space);
ISOSPEC_C_API bool advanceToNextConfigurationIsoStochasticGenerator(void* generator);
ISOSPEC_C_API void deleteIsoStochasticGenerator(void* generator);

ISOSPEC_C_API void* setupThresholdFixedEnvelope(void* iso,
                              double threshold,
                              bool absolute,
                              bool get_confs);

ISOSPEC_C_API void* setupTotalProbFixedEnvelope(void* iso,
                              double taget_coverage,
                              bool optimize,
                              bool get_confs);

ISOSPEC_C_API void* setupStochasticFixedEnvelope(void* iso,
                              size_t no_molecules,
                              double precision,
                              double beta_bias,
                              bool get_confs);

ISOSPEC_C_API void* setupBinnedFixedEnvelope(void* iso,
                    double target_total_prob,
                    double bin_width,
                    double bin_middle);

ISOSPEC_C_API void freeReleasedArray(void* array);

ISOSPEC_C_API void array_add(double* array, size_t N, double what);
ISOSPEC_C_API void array_mul(double* array, size_t N, double what);
ISOSPEC_C_API void array_fma(double* array, size_t N, double mul, double add);

ISOSPEC_C_API void* setupFixedEnvelope(double* masses, double* probs, size_t size, bool mass_sorted, bool prob_sorted, double total_prob);
ISOSPEC_C_API void* setupFixedEnvelopeWithConfs(double* masses, double* probs, int* confs, size_t size, int allDim, bool mass_sorted, bool prob_sorted, double total_prob);
ISOSPEC_C_API void* copyFixedEnvelope(void* other);
ISOSPEC_C_API void deleteFixedEnvelope(void* tabulator, bool releaseEverything);

/* These three hand the envelope's array over to the caller, who must free() it
   (with freeReleasedArray) once done. That costs a copy of the whole array
   whenever the envelope did not allocate it with plain malloc() -- which is the
   normal case now that everything the library produces is SIMD-aligned, and one
   that grew large enough is mapped straight from the OS. Prefer the *WithDeleter
   variants below, which never copy. */
ISOSPEC_C_API const double* massesFixedEnvelope(void* tabulator);
ISOSPEC_C_API const double* probsFixedEnvelope(void* tabulator);
ISOSPEC_C_API const int*    confsFixedEnvelope(void* tabulator);
ISOSPEC_C_API size_t confs_noFixedEnvelope(void* tabulator);

/* How a buffer obtained from one of the *WithDeleter entry points below has to
   be given back: deleter(array, size), exactly once, with the very size that
   came out alongside the pointer. free() is not, in general, the right answer.
   Call it through freeReleasedArrayWithDeleter() rather than directly, so
   callers that cannot invoke a raw function pointer (cffi, say) need not. */
typedef void (*IsoSpecArrayDeleter)(void* array, size_t size);

/* Zero-copy counterparts of the three getters above: each yields the array's
   pointer and writes out the (size, deleter) pair needed to release it. Both
   out-parameters may be NULL if not wanted -- but then the array can only be
   leaked. A NULL return means either an empty envelope or an error. */
ISOSPEC_C_API double* massesFixedEnvelopeWithDeleter(void* tabulator, size_t* size_out, IsoSpecArrayDeleter* deleter_out);
ISOSPEC_C_API double* probsFixedEnvelopeWithDeleter(void* tabulator, size_t* size_out, IsoSpecArrayDeleter* deleter_out);
ISOSPEC_C_API int*    confsFixedEnvelopeWithDeleter(void* tabulator, size_t* size_out, IsoSpecArrayDeleter* deleter_out);

ISOSPEC_C_API void freeReleasedArrayWithDeleter(void* array, size_t size, IsoSpecArrayDeleter deleter);

ISOSPEC_C_API double empiricAverageMass(void* tabulator);
ISOSPEC_C_API double empiricVariance(void* tabulator);
ISOSPEC_C_API double empiricStddev(void* tabulator);

ISOSPEC_C_API double wassersteinDistance(void* tabulator1, void* tabulator2);
ISOSPEC_C_API double orientedWassersteinDistance(void* tabulator1, void* tabulator2);
ISOSPEC_C_API double abyssalWassersteinDistance(void* tabulator1, void* tabulator2, double abyss_depth, double other_scale);
// ISOSPEC_C_API double abyssalWassersteinDistanceGrad(void* const* envelopes, const double* scales, double* ret_gradient, size_t N, double abyss_depth_exp, double abyss_depth_the);

struct ws_match_res{
double res1;
double res2;
double flow;
};

ISOSPEC_C_API struct ws_match_res wassersteinMatch(void* tabulator1, void* tabulator2, double flow_dist, double other_scale);

ISOSPEC_C_API void* addEnvelopes(void* tabulator1, void* tabulator2);
ISOSPEC_C_API void* convolveEnvelopes(void* tabulator1, void* tabulator2);

ISOSPEC_C_API double getTotalProbOfEnvelope(void* envelope);
ISOSPEC_C_API void scaleEnvelope(void* envelope, double factor);
ISOSPEC_C_API void normalizeEnvelope(void* envelope);
ISOSPEC_C_API void shiftMassEnvelope(void* envelope, double d_mass);
ISOSPEC_C_API void resampleEnvelope(void* envelope, size_t ionic_current, double beta_bias);
ISOSPEC_C_API void* binnedEnvelope(void* envelope, double width, double middle);
ISOSPEC_C_API void* linearCombination(void* const * const envelopes, const double* intensities, size_t count);

ISOSPEC_C_API void sortEnvelopeByMass(void* envelope);
ISOSPEC_C_API void sortEnvelopeByProb(void* envelope);

ISOSPEC_C_API void parse_fasta_c(const char* fasta, int atomCounts[6]);

/* Unimod-modification-aware peptide sequence parsing (fasta_mods.h/unimod.h):
   recognizes [UNIMOD:<id>] brackets in N-terminal ("[id]-SEQ"), internal
   ("X[id]"), and C-terminal ("SEQ-[id]") placement -- see docs/ai/unimod.md.
   unimod_db_path == NULL or "" selects the compile-time embedded default
   table; any other path is loaded and cached (see unimod_table_for_path).
   NULL on error (malformed bracket, unknown/excluded id) -- callers must
   NULL-check like every other handle-returning function here. */
ISOSPEC_C_API void* isoFromFastaWithMods(const char* sequence, bool use_nominal_masses, bool add_water, const char* unimod_db_path);

/* Opaque parsed-composition handle: NULL_UNIMOD_db_path selects the embedded
   default table, same as above. NULL return on error. */
ISOSPEC_C_API void* parseFastaWithModsC(const char* sequence, const char* unimod_db_path);
ISOSPEC_C_API size_t compositionSizeC(void* composition);
/* Element symbols, e.g. "C", "H", "Se" -- static strings (elem_table_symbol),
   must NOT be freed individually; only the returned array belongs to the
   handle, freed by deleteCompositionC. */
ISOSPEC_C_API const char* const* compositionSymbolsC(void* composition);
ISOSPEC_C_API const int* compositionCountsC(void* composition);
ISOSPEC_C_API void deleteCompositionC(void* composition);

/* What the library's batched kernels are vectorised to on this machine, as a
   stable lowercase token: "scalar", "sse2", "avx", "avx2", "avx512", "neon", or
   "simd" for a vector unit with no name here. Never NULL; a static string, so
   it must NOT be freed. See active_simd_level() in isa_kernels.h. */
ISOSPEC_C_API const char* activeSimdLevel(void);


#ifdef __cplusplus
}
#endif
