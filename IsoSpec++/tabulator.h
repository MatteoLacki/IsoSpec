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

class TabulatorParentCls {
public:
    double* _masses;
    double* _lprobs;
    double* _probs;
    int*    _confs;
    size_t  _confs_no;
    int     allDim;

    TabulatorParentCls() : _masses(nullptr),
        _lprobs(nullptr),
        _probs(nullptr),
        _confs(nullptr),
        _confs_no(0)
        {};

    virtual ~TabulatorParentCls()
    {
        if( _masses != nullptr ) free(_masses);
        if( _lprobs != nullptr ) free(_lprobs);
        if( _probs  != nullptr ) free(_probs);
        if( _confs  != nullptr ) free(_confs);
    };

    inline size_t    confs_no() const { return _confs_no; };
    inline int       getAllDim() const { return allDim; };

    // These are deliberately non-virtual
    inline double*   lprobs(bool release = false)   { double* ret = _lprobs; if(release) _lprobs = nullptr; return ret; };
    inline double*   masses(bool release = false)   { double* ret = _masses; if(release) _masses = nullptr; return ret; };
    inline double*   probs(bool release = false)    { double* ret = _probs;  if(release) _probs  = nullptr; return ret; };
    inline int*      confs(bool release = false)    { int*    ret = _confs;  if(release) _confs  = nullptr; return ret; };
};

class Tabulator : public TabulatorParentCls
{
public:
    Tabulator() : TabulatorParentCls() {};

    virtual ~Tabulator() {}

/*
    inline double*   lprobs(bool release = false)   { if constexpr(tgetlProbs) { double* ret = _lprobs; if(release) _lprobs = nullptr; return ret; } else throw std::logic_error("Logprobs requested, and yet tgetLprobs template argument was false"); };
    inline double*   masses(bool release = false)   { if constexpr(tgetMasses) { double* ret = _masses; if(release) _masses = nullptr; return ret; } else throw std::logic_error("Masses requested, and yet tgetMasses template argument was false"); };
    inline double*   probs(bool release = false)    { if constexpr(tgetProbs)  { double* ret = _probs;  if(release) _probs  = nullptr; return ret; } else throw std::logic_error("Probabilities requested, and yet tgetProbs template argument was false"); };
    inline int*      confs(bool release = false)    { if constexpr(tgetConfs)  { int*    ret = _confs;  if(release) _confs  = nullptr; return ret; } else throw std::logic_error("Configurations requested, and yet tgetConfs template argument was false"); };
*/
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

template<typename T> void call_init(T* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs)
{
    if(tgetlProbs)
    {
        if(tgetMasses)
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<true, true, true, true>(std::move(iso));
                else
                    tabulator->template init<true, true, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<true, true, false, true>(std::move(iso));
                else
                    tabulator->template init<true, true, false, false>(std::move(iso));
            }
        }
        else
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<true, false, true, true>(std::move(iso));
                else
                    tabulator->template init<true, false, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<true, false, false, true>(std::move(iso));
                else
                    tabulator->template init<true, false, false, false>(std::move(iso));
            }
        }
    }
    else
    {
        if(tgetMasses)
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<false, true, true, true>(std::move(iso));
                else
                    tabulator->template init<false, true, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<false, true, false, true>(std::move(iso));
                else
                    tabulator->template init<false, true, false, false>(std::move(iso));
            }
        }
        else
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<false, false, true, true>(std::move(iso));
                else
                    tabulator->template init<false, false, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<false, false, false, true>(std::move(iso));
                else
                    tabulator->template init<false, false, false, false>(std::move(iso));
            }
        }
    }
}

class ThresholdTabulator : public Tabulator
{
    const double threshold;
    const bool absolute;
public:
    ThresholdTabulator(Iso&& iso, double _threshold, bool _absolute, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs) :
    Tabulator(),
    threshold(_threshold),
    absolute(_absolute)
    {
        call_init<ThresholdTabulator>(this, std::move(iso), tgetlProbs, tgetMasses, tgetProbs, tgetConfs);
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

    virtual ~ThresholdTabulator() {};
};


class LayeredTabulator : public Tabulator
{
    const bool optimize;
public:
    LayeredTabulator(Iso&& iso, double _target_total_prob, bool _optimize, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs) : 
    Tabulator(),
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

    template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void init(Iso&& iso)
    {
        if(optimize && !tgetProbs)
        // If we want to optimize, we need the probs
            throw std::logic_error("Cannot perform quicktrim if we're not computing probabilities");

        IsoLayeredGenerator generator(std::move(iso), 1000, 1000, std::max<double>(target_total_prob, 1.0));

        this->allDim = generator.getAllDim();
        this->allDimSizeofInt = this->allDim*sizeof(int);


        this->reallocate_memory<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(ISOSPEC_INIT_TABLE_SIZE);

        size_t last_switch = 0;
        double prob_at_last_switch = 0.0;
        double prob_so_far = 0.0;
        do
        { // Store confs until we accumulate more prob than needed - and, if optimizing,
          // store also the rest of the last layer
            while(generator.advanceToNextConfigurationWithinLayer())
            {
                this->template addConf<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(generator);
                prob_so_far += generator.prob();
                if(prob_so_far >= target_total_prob)
                {
                    if (optimize)
                    {
                        while(generator.advanceToNextConfigurationWithinLayer())
                            this->template addConf<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(generator);
                        break;
                    }
                    else
                        return;
                }
            }
            if(prob_so_far >= target_total_prob)
                break;

            last_switch = this->_confs_no;
            prob_at_last_switch = prob_so_far;
        } while(generator.nextLayer(-3.0));

        if(!optimize || prob_so_far <= target_total_prob)
            return;

        // Right. We have extra configurations and we have been asked to produce an optimal p-set, so
        // now we shall trim unneeded configurations, using an algorithm dubbed "quicktrim"
        // - similar to the quickselect algorithm, except that we use the cumulative sum of elements
        // left of pivot to decide whether to go left or right, instead of the positional index.
        // We'll be sorting by the prob array, permuting the other ones in parallel.

        int* conf_swapspace = nullptr;
        if constexpr(tgetConfs)
            conf_swapspace = (int*) malloc(this->allDimSizeofInt);

        size_t start = last_switch;
        size_t end = this->_confs_no;
        double sum_to_start = prob_at_last_switch;

        while(start < end)
        {
            // Partition part
            size_t len = end - start;
#if ISOSPEC_BUILDING_R
            size_t pivot = len/2 + start;
#else
            size_t pivot = rand() % len + start;
#endif
            double pprob = this->_probs[pivot];
            swap<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(pivot, end-1, conf_swapspace);

            double new_csum = sum_to_start;

            size_t loweridx = start;
            for(size_t ii=start; ii<end-1; ii++)
                if(this->_probs[ii] > pprob)
                {
                    swap<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(ii, loweridx, conf_swapspace);
                    new_csum += this->_probs[loweridx];
                    loweridx++;
                }

            swap<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(end-1, loweridx, conf_swapspace);

            // Selection part
            if(new_csum < target_total_prob)
            {
                start = loweridx + 1;
                sum_to_start = new_csum + this->_probs[loweridx];
            }
            else
                end = loweridx;
        }

        if constexpr(tgetConfs)
            free(conf_swapspace);

        if(end <= current_size/2)
            // Overhead in memory of 2x or more, shrink to fit
            this->template reallocate_memory<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(end);

        this->_confs_no = end;
    }


    virtual ~LayeredTabulator() {};

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

};


} // namespace IsoSpec

