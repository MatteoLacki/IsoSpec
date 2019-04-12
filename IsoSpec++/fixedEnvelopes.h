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

class FixedEnvelope {
public:
    double* _masses;
    double* _lprobs;
    double* _probs;
    int*    _confs;
    size_t  _confs_no;
    int     allDim;

    FixedEnvelope() : _masses(nullptr),
        _lprobs(nullptr),
        _probs(nullptr),
        _confs(nullptr),
        _confs_no(0)
        {};

    virtual ~FixedEnvelope()
    {
        if( _masses != nullptr ) free(_masses);
        if( _lprobs != nullptr ) free(_lprobs);
        if( _probs  != nullptr ) free(_probs);
        if( _confs  != nullptr ) free(_confs);
    };

    inline size_t    confs_no() const { return _confs_no; };
    inline int       getAllDim() const { return allDim; };

    inline double*   lprobs(bool release = false)   { double* ret = _lprobs; if(release) _lprobs = nullptr; return ret; };
    inline double*   masses(bool release = false)   { double* ret = _masses; if(release) _masses = nullptr; return ret; };
    inline double*   probs(bool release = false)    { double* ret = _probs;  if(release) _probs  = nullptr; return ret; };
    inline int*      confs(bool release = false)    { int*    ret = _confs;  if(release) _confs  = nullptr; return ret; };

    inline double    mass(size_t i)  { return _masses[i]; };
    inline double    lprob(size_t i) { return _lprobs[i]; };
    inline double    prob(size_t i)  { return _probs[i];  };
    inline int*      conf(size_t i) { return _confs + i*allDim; };

protected:
    double* tmasses;
    double* tlprobs;
    double* tprobs;
    int*    tconfs;

    int allDimSizeofInt;

    template<typename T, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> ISOSPEC_FORCE_INLINE void store_conf(T& generator)
    {
        if constexpr(tgetlProbs) { *tlprobs = generator.lprob(); tlprobs++; };
        if constexpr(tgetMasses) { *tmasses = generator.mass();  tmasses++; };
        if constexpr(tgetProbs)  { *tprobs  = generator.prob();  tprobs++;  };
        if constexpr(tgetConfs)  { generator.get_conf_signature(tconfs); tconfs += allDim; };
    }

    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void reallocate_memory(size_t new_size)
    {
        // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...
        if constexpr(tgetlProbs) { _lprobs = (double*) realloc(_lprobs, new_size * sizeof(double)); tlprobs = _lprobs + _confs_no; }
        if constexpr(tgetMasses) { _masses = (double*) realloc(_masses, new_size * sizeof(double)); tmasses = _masses + _confs_no; }
        if constexpr(tgetProbs)  { _probs  = (double*) realloc(_probs,  new_size * sizeof(double));  tprobs  = _probs  + _confs_no; }
        if constexpr(tgetConfs)  { _confs  = (int*)    realloc(_confs,  new_size * allDimSizeofInt); tconfs = _confs + (allDim * _confs_no); }
    }

};

template<typename T> void call_init(T* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);

class ThresholdFixedEnvelope : public FixedEnvelope
{
    const double threshold;
    const bool absolute;
public:
    ThresholdFixedEnvelope(Iso&& iso, double _threshold, bool _absolute, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs) :
    FixedEnvelope(),
    threshold(_threshold),
    absolute(_absolute)
    {
        call_init<ThresholdFixedEnvelope>(this, std::move(iso), tgetlProbs, tgetMasses, tgetProbs, tgetConfs);
    }

    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void init(Iso&& iso)
    {
        IsoThresholdGenerator generator(std::move(iso), threshold, absolute);

        size_t tab_size = generator.count_confs();
        this->allDim = generator.getAllDim();
        this->allDimSizeofInt = this->allDim * sizeof(int);

        this->reallocate_memory<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(tab_size);

        while(generator.advanceToNextConfiguration())
            store_conf<IsoThresholdGenerator, tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(generator);

        this->_confs_no = tab_size;
    }

    virtual ~ThresholdFixedEnvelope() {};

    template<typename T> friend void call_init(T* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);

};


class TotalProbFixedEnvelope : public FixedEnvelope
{
    const bool optimize;
public:
    TotalProbFixedEnvelope(Iso&& iso, double _target_total_prob, bool _optimize, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs) :
    FixedEnvelope(),
    optimize(_optimize),
    target_total_prob(_target_total_prob >= 1.0 ? std::numeric_limits<double>::infinity() : _target_total_prob),
    current_size(ISOSPEC_INIT_TABLE_SIZE)
    {
        if(_target_total_prob <= 0.0)
            return;

        call_init(this, std::move(iso), tgetlProbs, tgetMasses, tgetProbs || _optimize, tgetConfs);

        if(!tgetProbs && optimize)
        {
            delete[] _probs;
            _probs = nullptr;
        }
    }

private:

    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void init(Iso&& iso);

public:
    virtual ~TotalProbFixedEnvelope() {};

private:
    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void swap([[maybe_unused]] size_t idx1, [[maybe_unused]] size_t idx2, [[maybe_unused]] int* conf_swapspace)
    {
        if constexpr(tgetlProbs) std::swap<double>(this->_lprobs[idx1], this->_lprobs[idx2]);
        if constexpr(tgetProbs)  std::swap<double>(this->_probs[idx1],  this->_probs[idx2]);
        if constexpr(tgetMasses) std::swap<double>(this->_masses[idx1], this->_masses[idx2]);
        if constexpr(tgetConfs)
        {
            int* c1 = this->_confs + (idx1*this->allDim);
            int* c2 = this->_confs + (idx2*this->allDim);
            memcpy(conf_swapspace, c1, this->allDimSizeofInt);
            memcpy(c1, c2, this->allDimSizeofInt);
            memcpy(c2, conf_swapspace, this->allDimSizeofInt);
        }
    }


    double target_total_prob;
    size_t current_size;


public:
    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void addConf(IsoLayeredGenerator& generator)
    {
        if(this->_confs_no == this->current_size)
        {
            this->current_size *= 2;
            this->template reallocate_memory<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(this->current_size);
        }

        this->template store_conf<IsoLayeredGenerator, tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(generator);
        this->_confs_no++;
    }

    template<typename T> friend void call_init(T* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);

};


} // namespace IsoSpec

