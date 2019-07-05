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

class ISOSPEC_EXPORT_SYMBOL FixedEnvelope {
protected:
    double* _masses;
    double* _lprobs;
    double* _probs;
    int*    _confs;
    size_t  _confs_no;
    int     allDim;

public:
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

    inline const double*   lprobs()   { return _lprobs; };
    inline const double*   masses()   { return _masses; };
    inline const double*   probs()    { return _probs; };
    inline const int*      confs()    { return _confs; };

    inline double*   release_lprobs()   { double* ret = _lprobs; _lprobs = nullptr; return ret; };
    inline double*   release_masses()   { double* ret = _masses; _masses = nullptr; return ret; };
    inline double*   release_probs()    { double* ret = _probs;  _probs  = nullptr; return ret; };
    inline int*      release_confs()    { int*    ret = _confs;  _confs  = nullptr; return ret; };


    inline double     mass(size_t i)  { return _masses[i]; };
    inline double     lprob(size_t i) { return _lprobs[i]; };
    inline double     prob(size_t i)  { return _probs[i];  };
    inline const int* conf(size_t i)  { return _confs + i*allDim; };

protected:
    double* tmasses;
    double* tlprobs;
    double* tprobs;
    int*    tconfs;

    int allDimSizeofInt;

    template<typename T, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> ISOSPEC_FORCE_INLINE void store_conf(T& generator)
    {
        constexpr_if(tgetlProbs) { *tlprobs = generator.lprob(); tlprobs++; };
        constexpr_if(tgetMasses) { *tmasses = generator.mass();  tmasses++; };
        constexpr_if(tgetProbs)  { *tprobs  = generator.prob();  tprobs++;  };
        constexpr_if(tgetConfs)  { generator.get_conf_signature(tconfs); tconfs += allDim; };
    }

    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void reallocate_memory(size_t new_size);
};

template<typename T> void call_init(T* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);

class ISOSPEC_EXPORT_SYMBOL ThresholdFixedEnvelope : public FixedEnvelope
{
    const double threshold;
    const bool absolute;
public:
    ThresholdFixedEnvelope(Iso&& iso, double _threshold, bool _absolute, bool tgetConfs = false, bool tgetlProbs = false, bool tgetMasses = true, bool tgetProbs = true) :
    FixedEnvelope(),
    threshold(_threshold),
    absolute(_absolute)
    {
        call_init<ThresholdFixedEnvelope>(this, std::move(iso), tgetlProbs, tgetMasses, tgetProbs, tgetConfs);
    }

    inline ThresholdFixedEnvelope(const Iso& iso, double _threshold, bool _absolute, bool tgetConfs = false, bool tgetlProbs = false, bool tgetMasses = true, bool tgetProbs = true) :
    ThresholdFixedEnvelope(Iso(iso, false), _threshold, _absolute, tgetConfs, tgetlProbs, tgetMasses, tgetProbs) {};

    virtual ~ThresholdFixedEnvelope() {};

private:
    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void init(Iso&& iso);

    template<typename T> friend void call_init(T* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);
};


class ISOSPEC_EXPORT_SYMBOL TotalProbFixedEnvelope : public FixedEnvelope
{
    const bool optimize;
    double target_total_prob;
    size_t current_size;

public:
    TotalProbFixedEnvelope(Iso&& iso, double _target_total_prob, bool _optimize, bool tgetConfs = false, bool tgetlProbs = false, bool tgetMasses = true, bool tgetProbs = true) :
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
            free(_probs);
            _probs = nullptr;
        }
    }

    inline TotalProbFixedEnvelope(const Iso& iso, double _target_total_prob, bool _optimize, bool tgetConfs = false, bool tgetlProbs = false, bool tgetMasses = true, bool tgetProbs = true) :
    TotalProbFixedEnvelope(Iso(iso, false), _target_total_prob, _optimize, tgetConfs, tgetlProbs, tgetMasses, tgetProbs) {};

    virtual ~TotalProbFixedEnvelope() {};

private:

    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void init(Iso&& iso);

    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void swap([[maybe_unused]] size_t idx1, [[maybe_unused]] size_t idx2, [[maybe_unused]] int* conf_swapspace)
    {
        constexpr_if(tgetlProbs) std::swap<double>(this->_lprobs[idx1], this->_lprobs[idx2]);
        constexpr_if(tgetProbs)  std::swap<double>(this->_probs[idx1],  this->_probs[idx2]);
        constexpr_if(tgetMasses) std::swap<double>(this->_masses[idx1], this->_masses[idx2]);
        constexpr_if(tgetConfs)
        {
            int* c1 = this->_confs + (idx1*this->allDim);
            int* c2 = this->_confs + (idx2*this->allDim);
            memcpy(conf_swapspace, c1, this->allDimSizeofInt);
            memcpy(c1, c2, this->allDimSizeofInt);
            memcpy(c2, conf_swapspace, this->allDimSizeofInt);
        }
    }


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

