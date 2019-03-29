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

#include "tabulator.h"


namespace IsoSpec
{

Tabulator::Tabulator() :
_masses(nullptr),
_lprobs(nullptr),
_probs(nullptr),
_confs(nullptr),
_confs_no(0),
mem_size(0)
{}

Tabulator::~Tabulator()
{
    if( _masses != nullptr ) free(_masses);
    if( _lprobs != nullptr ) free(_lprobs);
    if( _probs  != nullptr ) free(_probs);
    if( _confs  != nullptr ) free(_confs);
}


ThresholdTabulator::ThresholdTabulator(Iso&& iso, double threshold, bool absolute,
                                       bool get_masses, bool get_probs,
                                       bool get_lprobs, bool get_confs) :
Tabulator()
{
    IsoThresholdGenerator generator(std::move(iso), threshold, absolute);

    _confs_no = generator.count_confs();
    allDim = generator.getAllDim();
    allDimSizeofInt = allDim * sizeof(int);

    reallocate_memory(get_lprobs, get_probs, get_masses, get_confs, _confs_no);

    while(generator.advanceToNextConfiguration())
    {
        if(_masses != nullptr) { *_masses = generator.mass();  _masses++; }
        if(_lprobs != nullptr) { *_lprobs = generator.lprob(); _lprobs++; }
        if(_probs  != nullptr) { *_probs  = generator.prob();  _probs++;  }
        if(_confs  != nullptr){
            generator.get_conf_signature(_confs);
            _confs += allDim;
        }
    }

    if(_masses != nullptr) { _masses -= _confs_no; }
    if(_lprobs != nullptr) { _lprobs -= _confs_no; }
    if(_probs  != nullptr) { _probs -= _confs_no; }
    if(_confs  != nullptr) { _confs -= _confs_no*allDim; }

}

ThresholdTabulator::~ThresholdTabulator() {}



ISOSPEC_FORCE_INLINE void LayeredTabulator::swap(size_t idx1, size_t idx2, int* conf_swapspace)
{
    std::swap(_probs[idx1], _probs[idx2]);
    if(_lprobs != nullptr) std::swap<double>(_lprobs[idx1], _lprobs[idx2]);
    if(_masses != nullptr) std::swap<double>(_masses[idx1], _masses[idx2]);
    if(_confs != nullptr)
    {
        int* c1 = _confs + (idx1*allDim);
        int* c2 = _confs + (idx2*allDim);
        memcpy(conf_swapspace, c1, allDimSizeofInt);
        memcpy(c1, c2, allDimSizeofInt);
        memcpy(c2, conf_swapspace, allDimSizeofInt);
    }
}


template<bool t_get_lprobs, bool t_get_probs, bool t_get_masses, bool t_get_confs> double LayeredTabulator::layered_main_loop(IsoLayeredGenerator& generator, size_t& last_switch, double& prob_at_last_switch)
{
    double prob_so_far = 0.0;
    do
    { // Store confs until we accumulate more prob than needed - and, if optimizing,
      // store also the rest of the last layer
        while(generator.advanceToNextConfigurationWithinLayer())
        {
            addConf<t_get_lprobs, t_get_probs, t_get_masses, t_get_confs>(generator);
            prob_so_far += generator.prob();
            if(prob_so_far >= target_total_prob)
            {
                if (optimize)
                {
                    // We're not at the end of the layer: extract the rest of it. No need to keep storing prob though
                    while(generator.advanceToNextConfigurationWithinLayer())
                        addConf<t_get_lprobs, t_get_probs, t_get_masses, t_get_confs>(generator);
                    break;
                }
                return prob_so_far;
            }
        }
        last_switch = _confs_no;
        prob_at_last_switch = prob_so_far;
    } while(generator.nextLayer(-3.0));

    return prob_so_far;
}

LayeredTabulator::LayeredTabulator(Iso&& iso,
                     bool get_masses, bool get_probs,
                     bool get_lprobs, bool get_confs,
                     double _target_total_prob, bool _optimize) :
Tabulator(),
target_total_prob(_target_total_prob >= 1.0 ? std::numeric_limits<double>::infinity() : _target_total_prob),
current_size(ISOSPEC_INIT_TABLE_SIZE),
optimize(_optimize)
{
    if(_target_total_prob <= 0.0)
        return;

    IsoLayeredGenerator generator(std::move(iso));

    allDim = generator.getAllDim();
    allDimSizeofInt = allDim*sizeof(int);

    bool user_wants_probs = get_probs;
    if(optimize)
    // If we want to optimize, we need the probs
        get_probs = true;

    reallocate_memory(get_lprobs, get_probs, get_masses, get_confs, ISOSPEC_INIT_TABLE_SIZE);

    tmasses = _masses;
    tlprobs = _lprobs;
    tprobs = _probs;
    tconfs = _confs;

    size_t last_switch = 0;
    double prob_at_last_switch = 0.0;

    double prob_so_far;
    
    if(get_lprobs)
    {
        if(get_probs)
        {
            if(get_masses)
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<true, true, true, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<true, true, true, false>(generator, last_switch, prob_at_last_switch);
            }
            else
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<true, true, false, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<true, true, false, false>(generator, last_switch, prob_at_last_switch);
            }
        }
        else
        {
            if(get_masses)
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<true, false, true, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<true, false, true, false>(generator, last_switch, prob_at_last_switch);
            }
            else
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<true, false, false, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<true, false, false, false>(generator, last_switch, prob_at_last_switch);
            }
        }
    }

    else

    {
        if(get_probs)
        {
            if(get_masses)
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<false, true, true, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<false, true, true, false>(generator, last_switch, prob_at_last_switch);
            }
            else
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<false, true, false, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<false, true, false, false>(generator, last_switch, prob_at_last_switch);
            }
        }
        else
        {
            if(get_masses)
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<false, false, true, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<false, false, true, false>(generator, last_switch, prob_at_last_switch);
            }
            else
            {
                if(get_confs)
                    prob_so_far = layered_main_loop<false, false, false, true>(generator, last_switch, prob_at_last_switch);
                else
                    prob_so_far = layered_main_loop<false, false, false, false>(generator, last_switch, prob_at_last_switch);
            }
        }
    }



    if(!optimize || prob_so_far <= target_total_prob)
        return;

    // Right. We have extra configurations and we have been asked to produce an optimal p-set, so
    // now we shall trim unneeded configurations, using an algorithm dubbed "quicktrim"
    // - similar to the quickselect algorithm, except that we use the cumulative sum of elements
    // left of pivot to decide whether to go left or right, instead of the positional index.
    // We'll be sorting by the prob array, permuting the other ones in parallel.

    int* conf_swapspace = nullptr;
    if(get_confs)
        conf_swapspace = (int*) malloc(allDimSizeofInt);

    size_t start = last_switch;
    size_t end = _confs_no;
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
        double pprob = _probs[pivot];
        swap(pivot, end-1, conf_swapspace);

        double new_csum = sum_to_start;

        size_t loweridx = start;
        for(size_t ii=start; ii<end-1; ii++)
            if(_probs[ii] > pprob)
            {
                swap(ii, loweridx, conf_swapspace);
                new_csum += _probs[loweridx];
                loweridx++;
            }

        swap(end-1, loweridx, conf_swapspace);

        // Selection part
        if(new_csum < target_total_prob)
        {
            start = loweridx + 1;
            sum_to_start = new_csum + _probs[loweridx];
        }
        else
            end = loweridx;
    }

    if(get_confs)
        free(conf_swapspace);

    if(!user_wants_probs)
    {
        free(_probs);
        _probs = nullptr;
    }

    _confs_no = end;

    if(end <= current_size/2)
        // Overhead in memory of 2x or more, shrink to fit
        reallocate_memory(get_lprobs, user_wants_probs, get_masses, get_confs, end);
}



template<bool t_get_lprobs, bool t_get_probs, bool t_get_masses, bool t_get_confs> void LayeredTabulator::addConf(IsoLayeredGenerator& generator)
{
    if( _confs_no == current_size )
    {
        current_size *= 2;

        reallocate_memory(t_get_lprobs, t_get_probs, t_get_masses, t_get_confs, current_size);
    }

    store_conf<t_get_lprobs, t_get_probs, t_get_masses, t_get_confs, IsoLayeredGenerator>(generator);

    _confs_no++;
}


LayeredTabulator::~LayeredTabulator() {}


void Tabulator::reallocate_memory(bool t_get_lprobs, bool t_get_probs, bool t_get_masses, bool t_get_confs, size_t new_size)
{
    // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...
    if(t_get_masses) { _masses = (double*) realloc(_masses, new_size * sizeof(double)); tmasses = _masses + _confs_no; }
    if(t_get_lprobs) { _lprobs = (double*) realloc(_lprobs, new_size * sizeof(double)); tlprobs = _lprobs + _confs_no; }
    if(t_get_probs)  { _probs  = (double*) realloc(_probs,  new_size * sizeof(double));  tprobs  = _probs  + _confs_no; }
    if(t_get_confs)  { _confs  = (int*)    realloc(_confs,  new_size * allDimSizeofInt); tconfs = _confs + (allDim * _confs_no); }
}


} // namespace IsoSpec
