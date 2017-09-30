
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


#ifndef CWRAPPER_H
#define CWRAPPER_H

#define ALGO_LAYERED 0
#define ALGO_ORDERED 1
#define ALGO_THRESHOLD_ABSOLUTE 2
#define ALGO_THRESHOLD_RELATIVE 3
#define ALGO_LAYERED_ESTIMATE 4


#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

void* setupIso( int             _dimNumber,
                const int*      _isotopeNumbers,
                const int*      _atomCounts,
                const double*   _isotopeMasses,
                const double*   _isotopeProbabilities);



// ================================================================


void* setupIsoThreshold( int             _dimNumber,
                         const int*      _isotopeNumbers,
                         const int*      _atomCounts,
                         const double*   _isotopeMasses,
                         const double*   _isotopeProbabilities,
                         const double    _threshold,
                         int             _absolute,
                         int             tabSize,
                         int             hashSize
);

// ATTENTION! BELOW THIS LINE MATTEO WAS CODING AND IT IS BETTER NOT TO COMPILE THAT
// BEAR IN MIND THAT THE ABOVE COMMENT WAS ALSO WRITTEN BY MATTEO
// AND THAT HE IS REALLY STRETCHING HIS ABILITY TO MAKE JOKES ABOUT
// HIMSELF TO THE LIMIT.

typedef struct MassSpectrum{double* masses; double* logprobs; int confs_no;} MassSpectrum;

#define C_HEADER(generatorType, dataType, method)\
dataType method##generatorType(void* generator);

#define C_HEADER_GET_CONF_SIGNATURE(generatorType)\
void method##generatorType(void* generator);

#define C_HEADERS(generatorType)\
C_HEADER(generatorType, double, mass) \
C_HEADER(generatorType, double, lprob) \
C_HEADER_GET_CONF_SIGNATURE(generatorType) \
C_HEADER(generatorType, bool, advanceToNextConfiguration) \
C_HEADER(generatorType, void, delete)


//______________________________________________________THRESHOLD GENERATOR
void* setupIsoThresholdGenerator(int dimNumber,
                                 const int* isotopeNumbers,
                                 const int* atomCounts,
                                 const double* isotopeMasses,
                                 const double* isotopeProbabilities,
                                 const double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize);
C_HEADERS(IsoThresholdGenerator)


//______________________________________________________LAYERED GENERATOR
void* setupIsoLayeredGenerator(int dimNumber,
                                 const int* isotopeNumbers,
                                 const int* atomCounts,
                                 const double* isotopeMasses,
                                 const double* isotopeProbabilities,
                                 const double threshold,
                                 bool _absolute,
                                 int _tabSize,
                                 int _hashSize);
C_HEADERS(IsoLayeredGenerator)

//______________________________________________________ORDERED GENERATOR
void* setupIsoOrderedGenerator(int dimNumber,
                               const int* isotopeNumbers,
                               const int* atomCounts,
                               const double* isotopeMasses,
                               const double* isotopeProbabilities,
                               int _tabSize,
                               int _hashSize);
C_HEADERS(IsoOrderedGenerator)



// Check if there is bool in CFFI
void* setupThresholdTabulator(void* generator,
                              bool  get_masses,
                              bool  get_probs,
                              bool  get_lprobs,
                              bool  get_confs);

void deleteThresholdTabulator(void* tabulator);

const double* massesThresholdTabulator(void* tabulator);
const double* lprobsThresholdTabulator(void* tabulator);
const double* probsThresholdTabulator(void* tabulator);
const int*    confsThresholdTabulator(void* tabulator);
int confs_noThresholdTabulator(void* tabulator);

#ifdef __cplusplus
}
#endif

#endif
