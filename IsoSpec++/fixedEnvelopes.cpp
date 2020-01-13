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

#include "fixedEnvelopes.h"

namespace IsoSpec
{

FixedEnvelope::FixedEnvelope(const FixedEnvelope& other) :
_masses(array_copy<double>(other._masses, other._confs_no)),
_lprobs(array_copy<double>(other._lprobs, other._confs_no)),
_probs(array_copy<double>(other._probs, other._confs_no)),
_confs(array_copy<int>(other._confs, other._confs_no*other.allDim)),
_confs_no(other._confs_no),
allDim(other.allDim),
sorted_by_mass(other.sorted_by_mass),
sorted_by_prob(other.sorted_by_prob),
total_prob(other.total_prob)
{}

FixedEnvelope::FixedEnvelope(FixedEnvelope&& other) :
_masses(other._masses),
_lprobs(other._lprobs),
_probs(other._probs),
_confs(other._confs),
_confs_no(other._confs_no),
allDim(other.allDim),
sorted_by_mass(other.sorted_by_mass),
sorted_by_prob(other.sorted_by_prob),
total_prob(other.total_prob)
{
other._masses = nullptr;
other._lprobs = nullptr;
other._probs  = nullptr;
other._confs  = nullptr;
other._confs_no = 0;
other.total_prob = 0.0;
}

FixedEnvelope::FixedEnvelope(double* masses, double* probs, size_t confs_no, bool masses_sorted, bool probs_sorted, double _total_prob) :
_masses(masses),
_lprobs(nullptr),
_probs(probs),
_confs(nullptr),
_confs_no(confs_no),
allDim(0),
sorted_by_mass(masses_sorted),
sorted_by_prob(probs_sorted),
total_prob(_total_prob)
{}

FixedEnvelope FixedEnvelope::operator+(const FixedEnvelope& other) const
{
    if(_confs_no > 0)
        if(_masses == nullptr || _probs == nullptr)
            throw std::logic_error("Probabilities and masses must be available for spectrum addition to be meaningful");
    if(other._confs_no > 0)
        if(other._masses == nullptr || other._probs == nullptr)
            throw std::logic_error("Probabilities and masses must be available for spectrum addition to be meaningful");

    double* nprobs  = (double*) malloc(sizeof(double) * (_confs_no+other._confs_no));
    double* nmasses = (double*) malloc(sizeof(double) * (_confs_no+other._confs_no));

    memcpy(nprobs,  _probs,  sizeof(double) * _confs_no);
    memcpy(nmasses, _masses, sizeof(double) * _confs_no);

    memcpy(nprobs+_confs_no,  other._probs,  sizeof(double) * other._confs_no);
    memcpy(nmasses+_confs_no, other._masses, sizeof(double) * other._confs_no);

    return FixedEnvelope(nmasses, nprobs, _confs_no + other._confs_no);
}

FixedEnvelope FixedEnvelope::operator*(const FixedEnvelope& other) const
{
    if(_confs_no > 0)
        if(_masses == nullptr || _probs == nullptr)
            throw std::logic_error("Probabilities and masses must be available for spectrum convolution to be meaningful");
    if(other._confs_no > 0)
        if(other._masses == nullptr || other._probs == nullptr)
            throw std::logic_error("Probabilities and masses must be available for spectrum convolution to be meaningful");

    double* nprobs =  (double*) malloc(sizeof(double) * _confs_no * other._confs_no);
    double* nmasses = (double*) malloc(sizeof(double) * _confs_no * other._confs_no);

    size_t tgt_idx = 0;

    for(size_t ii=0; ii<_confs_no; ii++)
        for(size_t jj=0; jj<other._confs_no; jj++)
        {
            nprobs[tgt_idx]  = _probs[ii]  * other._probs[jj];
            nmasses[tgt_idx] = _masses[ii] + other._masses[jj];
            tgt_idx++;
        }

    return FixedEnvelope(nmasses, nprobs, _confs_no + other._confs_no);
}

void FixedEnvelope::sort_by_mass()
{
    if(_masses == nullptr)
        throw std::logic_error("Can't sort by masses if masses have not been computed");

    if(sorted_by_mass)
        return;

    if((_probs == nullptr) && (_lprobs == nullptr) && (_confs == nullptr))
        std::sort(_masses, _masses + _confs_no);
    else
        sort_by(_masses);

    sorted_by_mass = true;
    sorted_by_prob = false;
}


void FixedEnvelope::sort_by_prob()
{
    if((_probs == nullptr) && (_lprobs == nullptr))
        throw std::logic_error("Can't sort by probabilities if neither probs nor logprobs have not been computed");

    if(sorted_by_prob)
        return;

    if((_masses == nullptr) && (_confs == nullptr))
    {
        if(_probs != nullptr)
            std::sort(_probs, _probs + _confs_no);
        if(_lprobs != nullptr)
            std::sort(_lprobs, _lprobs + _confs_no);
        return;
    }

    if(_probs == nullptr)
        sort_by(_lprobs);
    else
        sort_by(_probs);

    sorted_by_prob = true;
    sorted_by_mass = false;
}

template<typename T> void reorder_array(T* arr, size_t* order, size_t size, bool can_destroy = false)
{
    if(!can_destroy)
    {
        size_t* order_c = new size_t[size];
        memcpy(order_c, order, sizeof(size_t)*size);
        order = order_c;
    }

    for(size_t ii=0; ii<size; ii++)
        while(order[ii] != ii)
        {
            std::swap(arr[ii], arr[order[ii]]);
            std::swap(order[order[ii]], order[ii]);
        }

    if(!can_destroy)
        delete[] order;
}

void FixedEnvelope::sort_by(double* order)
{
    size_t* indices = new size_t[_confs_no];

    for(size_t ii=0; ii<_confs_no; ii++)
        indices[ii] = ii;

    std::sort<size_t*>(indices, indices + _confs_no, TableOrder<double>(order));

    size_t* inverse = new size_t[_confs_no];

    for(size_t ii=0; ii<_confs_no; ii++)
        inverse[indices[ii]] = ii;

    delete[] indices;

    if(_masses != nullptr) reorder_array(_masses, inverse, _confs_no);
    if(_probs  != nullptr) reorder_array(_probs,  inverse, _confs_no);
    if(_lprobs != nullptr) reorder_array(_lprobs, inverse, _confs_no);
    if(_confs  != nullptr)
    {
        const int allDimSizeofInt = sizeof(int) * allDim;
        int* swapspace = new int[allDim];
        for(size_t ii = 0; ii < _confs_no; ii++)
            while(order[ii] != ii)
            {
                memcpy(swapspace, &_confs[ii*allDim], allDimSizeofInt);
                memcpy(&_confs[ii*allDim], &_confs[inverse[ii]*allDim], allDimSizeofInt);
                memcpy(&_confs[inverse[ii]*allDim], swapspace, allDimSizeofInt);
            }
        delete[] swapspace;
    }
    delete[] inverse;
}


double FixedEnvelope::get_total_prob()
{
    if(_probs == nullptr)
        throw std::logic_error("Cannot compute total probability if probabilities have not been computed");
    if(std::isnan(total_prob))
    {
        total_prob = 0.0;
        for(size_t ii=0; ii<_confs_no; ii++)
            total_prob += _probs[ii];
    }
    return total_prob;
}

void FixedEnvelope::scale(double factor)
{
    for(size_t ii = 0; ii<_confs_no; ii++)
        _probs[ii] *= factor;
    total_prob *= factor;
}

void FixedEnvelope::normalize()
{
    double tp = get_total_prob();
    if(tp != 1.0)
    {
        scale(1.0/tp);
        total_prob = 1.0;
    }
}

FixedEnvelope FixedEnvelope::LinearCombination(const std::vector<const FixedEnvelope*>& spectra, const std::vector<double>& intensities)
{
    return LinearCombination(spectra.data(), intensities.data(), spectra.size());
}

FixedEnvelope FixedEnvelope::LinearCombination(const FixedEnvelope* const * spectra, const double* intensities, size_t size)
{
    size_t ret_size = 0;
    for(size_t ii = 0; ii < size; ii++)
        ret_size += spectra[ii]->_confs_no;

    double* newprobs = (double*) malloc(sizeof(double)*ret_size);
    double* newmasses = (double*) malloc(sizeof(double)*ret_size);

    size_t cntr = 0;
    for(size_t ii = 0; ii < size; ii++)
    {
        double mul = intensities[ii];
        for(size_t jj = 0; jj < spectra[ii]->_confs_no; jj++)
            newprobs[jj+cntr] = spectra[ii]->_probs[jj] * mul;
        memcpy(newmasses + cntr, spectra[ii]->_masses, sizeof(double) * spectra[ii]->_confs_no);
        cntr += spectra[ii]->_confs_no;
    }
    return FixedEnvelope(newmasses, newprobs, cntr);
}

double FixedEnvelope::WassersteinDistance(FixedEnvelope& other)
{
    double ret = 0.0;
    if((get_total_prob()*0.999>other.get_total_prob()) || (other.get_total_prob()>get_total_prob()*1.001))
        throw std::logic_error("Spectra must be normalized before computing Wasserstein Distance");

    if(_confs_no == 0 || other._confs_no == 0)
        return 0.0;

    sort_by_mass();
    other.sort_by_mass();

    size_t idx_this = 0;
    size_t idx_other = 0;

    double acc_prob = 0.0;
    double last_point = 0.0;


    while(idx_this < _confs_no && idx_other < other._confs_no)
    {
        if(_masses[idx_this] < other._masses[idx_other])
        {
            ret += (_masses[idx_this] - last_point) * std::abs(acc_prob);
            acc_prob += _probs[idx_this];
            last_point = _masses[idx_this];
            idx_this++;
        }
        else
        {
            ret += (other._masses[idx_other] - last_point) * std::abs(acc_prob);
            acc_prob -= other._probs[idx_other];
            last_point = other._masses[idx_other];
            idx_other++;
        }
    }

    while(idx_this < _confs_no)
    {
        ret += (_masses[idx_this] - last_point) * std::abs(acc_prob);
        acc_prob += _probs[idx_this];
        last_point = _masses[idx_this];
        idx_this++;
    }

    while(idx_other < other._confs_no)
    {
        ret += (other._masses[idx_other] - last_point) * std::abs(acc_prob);
        acc_prob -= other._probs[idx_other];
        last_point = other._masses[idx_other];
        idx_other++;
    }

    return ret;
}

FixedEnvelope FixedEnvelope::bin(double bin_width, double middle)
{
    sort_by_mass();

    FixedEnvelope ret;

    if(_confs_no == 0)
        return ret;

    ret.reallocate_memory<false, true, true, false>(ISOSPEC_INIT_TABLE_SIZE);
    ret.current_size = ISOSPEC_INIT_TABLE_SIZE;

    size_t ii = 0;

    double half_width = 0.5*bin_width;
    double hwmm = half_width-middle;

    while(ii < _confs_no)
    {
        double current_bin_middle = floor((_masses[ii]+hwmm)/bin_width)*bin_width + middle;
        double current_bin_end = current_bin_middle + half_width;
        double bin_prob = 0.0;

        while(ii < _confs_no && _masses[ii] <= current_bin_end)
        {
            bin_prob += _probs[ii];
            ii++;
        }
        ret.store_conf(current_bin_middle, bin_prob);
    }

    return ret;
}

template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void FixedEnvelope::reallocate_memory(size_t new_size)
{
    // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...
    constexpr_if(tgetlProbs) { _lprobs = (double*) realloc(_lprobs, new_size * sizeof(double)); tlprobs = _lprobs + _confs_no; }
    constexpr_if(tgetMasses) { _masses = (double*) realloc(_masses, new_size * sizeof(double)); tmasses = _masses + _confs_no; }
    constexpr_if(tgetProbs)  { _probs  = (double*) realloc(_probs,  new_size * sizeof(double)); tprobs  = _probs  + _confs_no; }
    constexpr_if(tgetConfs)  { _confs  = (int*)    realloc(_confs,  new_size * allDimSizeofInt); tconfs = _confs + (allDim * _confs_no); }
}

void FixedEnvelope::slow_reallocate_memory(size_t new_size)
{
    // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...
    if(_lprobs != nullptr) { _lprobs = (double*) realloc(_lprobs, new_size * sizeof(double)); tlprobs = _lprobs + _confs_no; }
    if(_masses != nullptr) { _masses = (double*) realloc(_masses, new_size * sizeof(double)); tmasses = _masses + _confs_no; }
    if(_probs  != nullptr) { _probs  = (double*) realloc(_probs,  new_size * sizeof(double)); tprobs  = _probs  + _confs_no; }
    if(_confs  != nullptr) { _confs  = (int*)    realloc(_confs,  new_size * allDimSizeofInt); tconfs = _confs + (allDim * _confs_no); }
}

template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void ThresholdFixedEnvelope::init(Iso&& iso)
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


template<bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs> void TotalProbFixedEnvelope::init(Iso&& iso)
{
    if(optimize && !tgetProbs)
    // If we want to optimize, we need the probs
        throw std::logic_error("Cannot perform quicktrim if we're not computing probabilities");

    IsoLayeredGenerator generator(std::move(iso), 1000, 1000, true, std::min<double>(target_total_prob, 0.9999));

    this->allDim = generator.getAllDim();
    this->allDimSizeofInt = this->allDim*sizeof(int);


    this->reallocate_memory<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(ISOSPEC_INIT_TABLE_SIZE);

    size_t last_switch = 0;
    double prob_at_last_switch = 0.0;
    double prob_so_far = 0.0;
    double layer_delta;

    const double sum_above = log1p(-target_total_prob) - 2.3025850929940455; // log(0.1);

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

        layer_delta = sum_above - log1p(-prob_so_far);
        layer_delta = std::max(std::min(layer_delta, -0.1), -5.0);
    } while(generator.nextLayer(layer_delta));

    if(!optimize || prob_so_far <= target_total_prob)
        return;

    // Right. We have extra configurations and we have been asked to produce an optimal p-set, so
    // now we shall trim unneeded configurations, using an algorithm dubbed "quicktrim"
    // - similar to the quickselect algorithm, except that we use the cumulative sum of elements
    // left of pivot to decide whether to go left or right, instead of the positional index.
    // We'll be sorting by the prob array, permuting the other ones in parallel.

    int* conf_swapspace = nullptr;
    constexpr_if(tgetConfs)
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

    constexpr_if(tgetConfs)
        free(conf_swapspace);

    if(end <= current_size/2)
        // Overhead in memory of 2x or more, shrink to fit
        this->template reallocate_memory<tgetlProbs, tgetMasses, tgetProbs, tgetConfs>(end);

    this->_confs_no = end;
}

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

template void call_init<TotalProbFixedEnvelope>(TotalProbFixedEnvelope* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);
template void call_init<ThresholdFixedEnvelope>(ThresholdFixedEnvelope* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);


} // namespace IsoSpec
