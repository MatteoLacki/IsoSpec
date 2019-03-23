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

#pragma once

#include <stdlib.h>

#include "isoSpec++.h"

#define ISOSPEC_INIT_TABLE_SIZE 1024

namespace IsoSpec
{


class Tabulator
{
protected:
    double* _masses;
    double* _lprobs;
    double* _probs;
    int*    _confs;
    size_t  _confs_no;
    int     allDim;
public:
    Tabulator();
    virtual ~Tabulator();

    inline double*   masses(bool release = false)   { double* ret = _masses; if(release) _masses = nullptr; return ret; };
    inline double*   lprobs(bool release = false)   { double* ret = _lprobs; if(release) _lprobs = nullptr; return ret; };
    inline double*   probs(bool release = false)    { double* ret = _probs;  if(release) _probs  = nullptr; return ret; };
    inline int*      confs(bool release = false)    { int*    ret = _confs;  if(release) _confs  = nullptr; return ret; };
    inline size_t    confs_no() const { return _confs_no; };
    inline int       getAllDim() const { return allDim; };
};

inline void reallocate(double **array, int new_size){
    if( *array != nullptr ){
        *array = (double *) realloc(*array, new_size);
    }
}


class ThresholdTabulator : public Tabulator
{
public:
    ThresholdTabulator(Iso&& iso, double threshold, bool absolute,
                       bool get_masses, bool get_probs,
                       bool get_lprobs, bool get_confs);

    virtual ~ThresholdTabulator();
};



class LayeredTabulator : public Tabulator
{
public:
    LayeredTabulator(Iso&& iso,
                     bool get_masses, bool get_probs,
                     bool get_lprobs, bool get_confs,
                     double _total_prob, bool _optimize = false);

    virtual ~LayeredTabulator();

private:
    void swap(size_t idx1, size_t idx2, int* conf_swapspace);
    double target_total_prob;
    size_t current_size;
    bool optimize;
    void addConf(IsoLayeredGenerator& generator);

    double* tmasses;
    double* tprobs;
    double* tlprobs;
    int* tconfs;
    size_t allDimSizeofInt;
};

} // namespace IsoSpec

