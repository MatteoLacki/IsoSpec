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
#include <limits>
#include <memory>
#include <cassert>
#include "isoMath.h"
#include "aligned_ptr.h"

namespace IsoSpec
{

FixedEnvelope::FixedEnvelope(const FixedEnvelope& other) :
_masses(array_copy_malloc<double>(other._masses, other._confs_no)),
_probs(array_copy_malloc<double>(other._probs, other._confs_no)),
_confs(array_copy_nptr_malloc<int>(other._confs, other._confs_no*other.allDim)),
_confs_no(other._confs_no),
allDim(other.allDim),
sorted_by_mass(other.sorted_by_mass),
sorted_by_prob(other.sorted_by_prob),
total_prob(other.total_prob),
current_size(other._confs_no),
allDimSizeofInt(other.allDimSizeofInt)
{}

FixedEnvelope::FixedEnvelope(FixedEnvelope&& other) :
_masses(other._masses),
_probs(other._probs),
_confs(other._confs),
_confs_no(other._confs_no),
allDim(other.allDim),
sorted_by_mass(other.sorted_by_mass),
sorted_by_prob(other.sorted_by_prob),
total_prob(other.total_prob),
current_size(other.current_size),
allDimSizeofInt(other.allDimSizeofInt)
{
other._masses = nullptr;
other._probs  = nullptr;
other._confs  = nullptr;
other._confs_no = 0;
other.total_prob = 0.0;
other.current_size = 0;
}

FixedEnvelope& FixedEnvelope::operator=(const FixedEnvelope& other)
{
    if(this == &other)
        return *this;
    // Copy-and-swap: construct the copy first so that if any allocation inside
    // the copy constructor throws, *this is left untouched.  Only after tmp is
    // fully built do we pilfer its pointers; tmp then destructs holding the old
    // ones, which frees them without any possibility of failure.
    FixedEnvelope tmp(other);
    std::swap(_masses,         tmp._masses);
    std::swap(_probs,          tmp._probs);
    std::swap(_confs,          tmp._confs);
    _confs_no       = tmp._confs_no;
    allDim          = tmp.allDim;
    allDimSizeofInt = tmp.allDimSizeofInt;
    sorted_by_mass  = tmp.sorted_by_mass;
    sorted_by_prob  = tmp.sorted_by_prob;
    total_prob      = tmp.total_prob;
    current_size    = tmp.current_size;
    return *this;
}

FixedEnvelope& FixedEnvelope::operator=(FixedEnvelope&& other)
{
    if(this == &other)
        return *this;
    free(_masses);
    free(_probs);
    free(_confs);
    _masses         = other._masses;         other._masses      = nullptr;
    _probs          = other._probs;          other._probs       = nullptr;
    _confs          = other._confs;          other._confs       = nullptr;
    _confs_no       = other._confs_no;       other._confs_no    = 0;
    allDim          = other.allDim;
    allDimSizeofInt = other.allDimSizeofInt;
    sorted_by_mass  = other.sorted_by_mass;
    sorted_by_prob  = other.sorted_by_prob;
    total_prob      = other.total_prob;      other.total_prob   = 0.0;
    current_size    = other.current_size;    other.current_size = 0;
    return *this;
}

FixedEnvelope::FixedEnvelope(double* in_masses, double* in_probs, size_t in_confs_no, bool masses_sorted, bool probs_sorted, double _total_prob) :
_masses(in_masses),
_probs(in_probs),
_confs(nullptr),
_confs_no(in_confs_no),
allDim(0),
sorted_by_mass(masses_sorted),
sorted_by_prob(probs_sorted),
total_prob(_total_prob),
allDimSizeofInt(0)
{}

FixedEnvelope::FixedEnvelope(double* in_masses, double* in_probs, int* in_confs, size_t in_confs_no, int _allDim, bool masses_sorted, bool probs_sorted, double _total_prob) :
_masses(in_masses),
_probs(in_probs),
_confs(in_confs),
_confs_no(in_confs_no),
allDim(_allDim),
sorted_by_mass(masses_sorted),
sorted_by_prob(probs_sorted),
total_prob(_total_prob),
allDimSizeofInt(_allDim * sizeof(int))
{}

FixedEnvelope FixedEnvelope::operator+(const FixedEnvelope& other) const
{
    double* nprobs  = reinterpret_cast<double*>(malloc(sizeof(double) * (_confs_no+other._confs_no)));
    if(nprobs == nullptr)
        throw std::bad_alloc();
    double* nmasses = reinterpret_cast<double*>(malloc(sizeof(double) * (_confs_no+other._confs_no)));
    if(nmasses == nullptr)
    {
        free(nprobs);
        throw std::bad_alloc();
    }

    // An empty envelope has null arrays, and memcpy() is declared nonnull: a
    // zero-length copy from it is undefined behaviour, so skip it entirely.
    if(_confs_no > 0)
    {
        memcpy(nprobs,  _probs,  sizeof(double) * _confs_no);
        memcpy(nmasses, _masses, sizeof(double) * _confs_no);
    }

    if(other._confs_no > 0)
    {
        memcpy(nprobs+_confs_no,  other._probs,  sizeof(double) * other._confs_no);
        memcpy(nmasses+_confs_no, other._masses, sizeof(double) * other._confs_no);
    }

    return FixedEnvelope(nmasses, nprobs, _confs_no + other._confs_no);
}

FixedEnvelope FixedEnvelope::operator*(const FixedEnvelope& other) const
{
    double* nprobs =  reinterpret_cast<double*>(malloc(sizeof(double) * _confs_no * other._confs_no));
    if(nprobs == nullptr)
        throw std::bad_alloc();
    //  deepcode ignore CMemoryLeak: It's not a memleak: the memory is passed to FixedEnvelope which
    //  deepcode ignore CMemoryLeak: takes ownership of it, and will properly free() it in destructor.
    double* nmasses = reinterpret_cast<double*>(malloc(sizeof(double) * _confs_no * other._confs_no));
    if(nmasses == nullptr)
    {
        free(nprobs);
        throw std::bad_alloc();
    }

    size_t tgt_idx = 0;

    for(size_t ii = 0; ii < _confs_no; ii++)
        for(size_t jj = 0; jj < other._confs_no; jj++)
        {
            nprobs[tgt_idx]  = _probs[ii]  * other._probs[jj];
            nmasses[tgt_idx] = _masses[ii] + other._masses[jj];
            tgt_idx++;
        }

    return FixedEnvelope(nmasses, nprobs, tgt_idx);
}

void FixedEnvelope::sort_by_mass()
{
    if(sorted_by_mass)
        return;

    sort_by(_masses);

    sorted_by_mass = true;
    sorted_by_prob = false;
}


void FixedEnvelope::sort_by_prob()
{
    if(sorted_by_prob)
        return;

    sort_by(_probs);

    sorted_by_prob = true;
    sorted_by_mass = false;
}

template<typename T> void reorder_array(T* arr, size_t* order, size_t size, bool can_destroy = false)
{
    std::unique_ptr<size_t[]> order_owned;
    if(!can_destroy)
    {
        order_owned = std::make_unique<size_t[]>(size);
        memcpy(order_owned.get(), order, sizeof(size_t)*size);
        order = order_owned.get();
    }

    for(size_t ii = 0; ii < size; ii++)
        while(order[ii] != ii)
        {
            std::swap(arr[ii], arr[order[ii]]);
            std::swap(order[order[ii]], order[ii]);
        }
}

void FixedEnvelope::sort_by(double* order)
{
    if(_confs_no <= 1)
        return;

    size_t* indices = new size_t[_confs_no];

    for(size_t ii = 0; ii < _confs_no; ii++)
        indices[ii] = ii;

    std::sort<size_t*>(indices, indices + _confs_no, TableOrder<double>(order));

    size_t* inverse = new size_t[_confs_no];

    for(size_t ii = 0; ii < _confs_no; ii++)
        inverse[indices[ii]] = ii;

    delete[] indices;

    reorder_array(_masses, inverse, _confs_no);
    reorder_array(_probs,  inverse, _confs_no, _confs == nullptr);
    if(_confs != nullptr)
    {
        int* swapspace = new int[allDim];
        for(size_t ii = 0; ii < _confs_no; ii++)
            while(inverse[ii] != ii)
            {
                memcpy(swapspace, &_confs[ii*allDim], allDimSizeofInt);
                memcpy(&_confs[ii*allDim], &_confs[inverse[ii]*allDim], allDimSizeofInt);
                memcpy(&_confs[inverse[ii]*allDim], swapspace, allDimSizeofInt);
                std::swap(inverse[inverse[ii]], inverse[ii]);
            }
        delete[] swapspace;
    }
    delete[] inverse;
}


double FixedEnvelope::get_total_prob()
{
    if(std::isnan(total_prob))
    {
        total_prob = 0.0;
        for(size_t ii = 0; ii < _confs_no; ii++)
            total_prob += _probs[ii];
    }
    return total_prob;
}

void FixedEnvelope::scale(double factor)
{
    for(size_t ii = 0; ii < _confs_no; ii++)
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

void FixedEnvelope::shift_mass(double value)
{
    for(size_t ii = 0; ii < _confs_no; ii++)
        _masses[ii] += value;
}

void FixedEnvelope::resample(size_t samples, double beta_bias)
{
    if(_confs_no == 0)
        throw std::logic_error("Resample called on an empty spectrum");

    // The probabilities are about to be replaced by molecule counts, so any
    // cached sum is stale; invalidate it and let get_total_prob() rescan.
    total_prob = NAN;

    double pprob = 0.0;
    double cprob = 0.0;
    size_t pidx = -1; // Overflows - but it doesn't matter.

    // Sentinel: prevents the inner while(pprob < cprob) from walking off the end
    // if floating-point rounding leaves the accumulated probability sum slightly
    // below 1.0. Safe to overwrite in-place: rdvariate_binom guards succ_prob>=1.0
    // so infinity is handled correctly, and the sentinel is always cleaned up before
    // return — either by the loop zeroing the slot when it advances pidx to here,
    // or by the memset below when the loop terminates earlier.
    _probs[_confs_no-1] = (std::numeric_limits<double>::max)();

    while(samples > 0)
    {
        pprob += _probs[++pidx];
        _probs[pidx] = 0.0;
        double covered_part = (pprob - cprob) / (1.0 - cprob);
        while(samples * covered_part < beta_bias && samples > 0)
        {
            cprob += rdvariate_beta_1_b(samples) * (1.0 - cprob);
            while(pprob < cprob)
            {
                pprob += _probs[++pidx];
                _probs[pidx] = 0.0;
            }
            _probs[pidx] += 1.0;
            samples--;
            covered_part = (pprob - cprob) / (1.0 - cprob);
        }
        if(samples <= 0)
                break;
        size_t nrtaken = rdvariate_binom(samples, covered_part);
        _probs[pidx] += static_cast<double>(nrtaken);
        samples -= nrtaken;
        cprob = pprob;
    }

    pidx++;
    memset(_probs + pidx, 0, sizeof(double)*(_confs_no - pidx));
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

    double* newprobs  = reinterpret_cast<double*>(malloc(sizeof(double)*ret_size));
    if(newprobs == nullptr)
        throw std::bad_alloc();
    double* newmasses = reinterpret_cast<double*>(malloc(sizeof(double)*ret_size));
    if(newmasses == nullptr)
    {
        free(newprobs);
        throw std::bad_alloc();
    }

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
    if((get_total_prob()*0.999 > other.get_total_prob()) || (other.get_total_prob() > get_total_prob()*1.001))
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

    acc_prob = std::abs(acc_prob);

    while(idx_this < _confs_no)
    {
        ret += (_masses[idx_this] - last_point) * acc_prob;
        acc_prob -= _probs[idx_this];
        last_point = _masses[idx_this];
        idx_this++;
    }

    while(idx_other < other._confs_no)
    {
        ret += (other._masses[idx_other] - last_point) * acc_prob;
        acc_prob -= other._probs[idx_other];
        last_point = other._masses[idx_other];
        idx_other++;
    }

    return ret;
}


double FixedEnvelope::OrientedWassersteinDistance(FixedEnvelope& other)
{
    double ret = 0.0;
    if((get_total_prob()*0.999 > other.get_total_prob()) || (other.get_total_prob() > get_total_prob()*1.001))
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
            ret += (_masses[idx_this] - last_point) * acc_prob;
            acc_prob += _probs[idx_this];
            last_point = _masses[idx_this];
            idx_this++;
        }
        else
        {
            ret += (other._masses[idx_other] - last_point) * acc_prob;
            acc_prob -= other._probs[idx_other];
            last_point = other._masses[idx_other];
            idx_other++;
        }
    }

    // acc_prob is the signed CDF difference (this minus other), so peaks of
    // *this* add to it — in the main loop above and here in the tail alike.
    // (The unsigned version folds the sign into abs() before its tails and so
    // subtracts in both; this one must not.)
    while(idx_this < _confs_no)
    {
        ret += (_masses[idx_this] - last_point) * acc_prob;
        acc_prob += _probs[idx_this];
        last_point = _masses[idx_this];
        idx_this++;
    }

    while(idx_other < other._confs_no)
    {
        ret += (other._masses[idx_other] - last_point) * acc_prob;
        acc_prob -= other._probs[idx_other];
        last_point = other._masses[idx_other];
        idx_other++;
    }

    return ret;
}

double FixedEnvelope::AbyssalWassersteinDistance(FixedEnvelope& other, double abyss_depth, double other_scale)
{
    sort_by_mass();
    other.sort_by_mass();

    std::vector<std::pair<double, double>> carried;

    size_t idx_this = 0;
    size_t idx_other = 0;

    //std::cout << "AAA" << std::endl;

    auto finished = [&]() -> bool { return idx_this >= _confs_no && idx_other >= other._confs_no; };
    auto next = [&]() -> std::pair<double, double> {
                            if(idx_this >= _confs_no || (idx_other < other._confs_no && _masses[idx_this] > other._masses[idx_other]))
                            {
                                std::pair<double, double> res = std::pair<double, double>(other._masses[idx_other], other._probs[idx_other]*other_scale);
                                idx_other++;
                                return res;
                            }
                            else
                            {
                                std::pair<double, double> res = std::pair<double, double>(_masses[idx_this], -_probs[idx_this]);
                                idx_this++;
                                return res;
                            }
                        };
    double accd = 0.0;
    double condemned = 0.0;

    while(!finished())
    {
        auto pair = next();
        double m = pair.first;
        double p = pair.second;
        if(!carried.empty() && carried[0].second * p > 0.0)
        {
            carried.emplace_back(m, p);
            continue;
        }

        while(!carried.empty())
        {
            double cm = carried.back().first;
            double cp = carried.back().second;
            if(m - cm >= abyss_depth)
            {
                for(auto it = carried.cbegin(); it != carried.cend(); it++)
                    condemned += fabs(it->second);
                carried.clear();
                break;
            }
            if((cp+p)*p > 0.0)
            {
                accd += fabs((m-cm)*cp);
                p += cp;
                carried.pop_back();
            }
            else
            {
                // The incoming peak is exhausted against part of the carried
                // one; write the remainder back (cp is a copy of the carried
                // amount, so updating it alone would lose the match and the
                // matched mass would be counted a second time — as transport
                // against a later peak, or as condemned mass at the end).
                accd += fabs((m-cm)*p);
                cp += p;
                if(cp == 0.0)
                    carried.pop_back();
                else
                    carried.back().second = cp;
                p = 0.0;
                break;
            }
        }
        if(p != 0.0)
            carried.emplace_back(m, p);
        //std::cout << m << " " << p << std::endl;
    }

    for(auto it = carried.cbegin(); it != carried.cend(); it++)
        condemned += fabs(it->second);

    return accd + condemned * abyss_depth * 0.5;
}

std::tuple<double, double, double> FixedEnvelope::WassersteinMatch(FixedEnvelope& other, double flow_distance, double other_scale)
{
    if(_confs_no == 0)
        return {0.0, other.get_total_prob() * other_scale, 0.0};

    double unmatched1 = 0.0;
    double unmatched2 = 0.0;
    double massflow = 0.0;

    sort_by_mass();
    other.sort_by_mass();

    size_t idx_this = 0;
    size_t idx_other = 0;
    double used_prob_this = 0.0;
    double used_prob_other = 0.0;

    while(idx_this < _confs_no && idx_other < other._confs_no)
    {
        bool moved = true;
        while(moved && idx_this < _confs_no && idx_other < other._confs_no)
        {
            moved = false;
            if(_masses[idx_this] < other._masses[idx_other] - flow_distance)
            {
                unmatched1 += _probs[idx_this] - used_prob_this;
                used_prob_this = 0.0;
                idx_this++;
                moved = true;
            }
            // The branch above may have just consumed the last peak of *this*;
            // without the bound check the read below runs off the end of the
            // masses array (the leftover peaks of `other` are accounted for by
            // the tail loop after the outer while).
            if(idx_this < _confs_no && other._masses[idx_other] < _masses[idx_this] - flow_distance)
            {
                unmatched2 += other._probs[idx_other]*other_scale - used_prob_other;
                used_prob_other = 0.0;
                idx_other++;
                moved = true;
            }
        }
        if(idx_this < _confs_no && idx_other < other._confs_no)
        {
            assert(_probs[idx_this] - used_prob_this >= 0.0);
            assert(other._probs[idx_other]*other_scale - used_prob_other >= 0.0);

            if(_probs[idx_this] - used_prob_this < other._probs[idx_other]*other_scale - used_prob_other)
            {
                massflow += _probs[idx_this] - used_prob_this;
                used_prob_other += _probs[idx_this] - used_prob_this;
                assert(used_prob_other >= 0.0);
                used_prob_this = 0.0;
                idx_this++;
            }
            else
            {
                massflow += other._probs[idx_other]*other_scale - used_prob_other;
                used_prob_this += other._probs[idx_other]*other_scale - used_prob_other;
                assert(used_prob_this >= 0.0);
                used_prob_other = 0.0;
                idx_other++;
            }
        }
    }

    unmatched1 -= used_prob_this;
    unmatched2 -= used_prob_other;

    for(; idx_this < _confs_no; idx_this++)
        unmatched1 += _probs[idx_this];
    for(; idx_other < other._confs_no; idx_other++)
        unmatched2 += other._probs[idx_other]*other_scale;

    return {unmatched1, unmatched2, massflow};
}

FixedEnvelope FixedEnvelope::bin(double bin_width, double middle)
{
    sort_by_mass();

    FixedEnvelope ret;

    if(_confs_no == 0)
        return ret;

    ret.reallocate_memory<false>(ISOSPEC_INIT_TABLE_SIZE);

    if(bin_width == 0)
    {
        double curr_mass = _masses[0];
        double accd_prob = _probs[0];
        for(size_t ii = 1; ii<_confs_no; ii++)
        {
            if(curr_mass != _masses[ii])
            {
                ret.store_conf(curr_mass, accd_prob);
                curr_mass = _masses[ii];
                accd_prob = _probs[ii];
            }
            else
                accd_prob += _probs[ii];
        }
        ret.store_conf(curr_mass, accd_prob);
        return ret;
    }

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

template<bool tgetConfs> void FixedEnvelope::reallocate_memory(size_t new_size)
{
    if(new_size > SIZE_MAX / sizeof(double))
        throw std::bad_alloc();
    double* tmp_masses = reinterpret_cast<double*>(realloc(_masses, new_size * sizeof(double)));
    if(tmp_masses == nullptr)
        throw std::bad_alloc();
    _masses = tmp_masses;
    tmasses = _masses + _confs_no;

    double* tmp_probs = reinterpret_cast<double*>(realloc(_probs,  new_size * sizeof(double)));
    if(tmp_probs == nullptr)
        throw std::bad_alloc();
    _probs = tmp_probs;
    tprobs  = _probs  + _confs_no;

    constexpr_if(tgetConfs)
    {
        if(allDimSizeofInt > 0 && new_size > SIZE_MAX / static_cast<size_t>(allDimSizeofInt))
            throw std::bad_alloc();
        int* tmp_confs = reinterpret_cast<int*>(realloc(_confs,  new_size * allDimSizeofInt));
        if(tmp_confs == nullptr)
            throw std::bad_alloc();
        _confs = tmp_confs;
        tconfs = _confs + (allDim * _confs_no);
    }
    current_size = new_size;
}

template<bool tgetConfs> void FixedEnvelope::aligned_allocate_memory(size_t new_size)
{
    // aligned_unique_ptr's small-allocation backend is aligned_alloc-based
    // (overflow-checked) and release() hands back a plain, free()-compatible
    // T* -- exactly what aligned_alloc()+free() already gave this function,
    // just without duplicating the overflow-safe rounding logic here. This
    // never grows afterwards (see the declaration in fixedEnvelopes.h), so
    // there's nothing to gain from keeping the aligned_unique_ptr itself
    // around past this one allocation.
    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> masses_buf(new_size);
    _masses = masses_buf.release();
    tmasses = _masses + _confs_no;

    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> probs_buf(new_size);
    _probs = probs_buf.release();
    tprobs  = _probs  + _confs_no;

    constexpr_if(tgetConfs)
    {
        aligned_unique_ptr<int, DOUBLE_SIMD_ALIGNMENT> confs_buf(new_size * static_cast<size_t>(allDim));
        _confs = confs_buf.release();
        tconfs = _confs + (allDim * _confs_no);
    }
    current_size = new_size;
}

void FixedEnvelope::slow_reallocate_memory(size_t new_size)
{
    if(new_size > SIZE_MAX / sizeof(double))
        throw std::bad_alloc();
    double* tmp_masses = reinterpret_cast<double*>(realloc(_masses, new_size * sizeof(double)));
    if(tmp_masses == nullptr)
        throw std::bad_alloc();
    _masses = tmp_masses;
    tmasses = _masses + _confs_no;

    double* tmp_probs = reinterpret_cast<double*>(realloc(_probs,  new_size * sizeof(double)));
    if(tmp_probs == nullptr)
        throw std::bad_alloc();
    _probs = tmp_probs;
    tprobs  = _probs  + _confs_no;

    if(_confs != nullptr)
    {
        if(allDimSizeofInt > 0 && new_size > SIZE_MAX / static_cast<size_t>(allDimSizeofInt))
            throw std::bad_alloc();
        int* tmp_confs = reinterpret_cast<int*>(realloc(_confs,  new_size * allDimSizeofInt));
        if(tmp_confs == nullptr)
            throw std::bad_alloc();
        _confs = tmp_confs;
        tconfs = _confs + (allDim * _confs_no);
    }
    current_size = new_size;
}

template<bool tgetConfs> void FixedEnvelope::threshold_init(Iso&& iso, double threshold, bool absolute)
{
    IsoThresholdGenerator generator(std::move(iso), threshold, absolute);

    size_t tab_size = generator.count_confs();
    this->allDim = generator.getAllDim();
    this->allDimSizeofInt = this->allDim * sizeof(int);

    this->aligned_allocate_memory<tgetConfs>(tab_size);

    double* ttmasses = this->_masses;
    double* ttprobs = this->_probs;
    ISOSPEC_MAYBE_UNUSED int* ttconfs;
    constexpr_if(tgetConfs)
        ttconfs = _confs;

    constexpr_if(tgetConfs) {
        while(generator.advanceToNextConfiguration())
        {
            *ttmasses = generator.mass(); ttmasses++;
            *ttprobs = generator.prob(); ttprobs++;
            generator.get_conf_signature(ttconfs); ttconfs += allDim;
        }
    }
    else
    {
#if ISOSPEC_HAS_SIMD
        // SIMD fill. Each marginal-0 run (one per higher-dimensional carry state) starts
        // at index 0 and descends until it drops below the cutoff; we batch the bulk of a
        // run W-wide and drain the < W tail scalar.
        //
        // Convention bridge: simd_massprobs()/advanceToNextConfiguration_no_carry() use
        // "lProbs_ptr points at the last-emitted config" (advance, then emit), matching the
        // generator's initial position of one-before-index-0. carry() instead leaves
        // lProbs_ptr *at* the index-0 config to be emitted. So after each successful carry
        // we emit that index-0 config scalar (it is always above the cutoff when carry()
        // succeeds, exactly as the scalar path relies on), which restores the last-emitted
        // convention and lets the next simd_massprobs resume cleanly from index 1.
        simd_double simd_masses;
        simd_double simd_probs;
        do {
            while(generator.simd_massprobs(simd_masses, simd_probs))
            {
                // Output is packed contiguously and the scalar tail advances the pointers
                // by a non-multiple of W each run, so the store target is not W-aligned in
                // general -> element_aligned (unaligned) store.
                simd_masses.copy_to(ttmasses, simd_ns::element_aligned); ttmasses += simd_masses.size();
                simd_probs.copy_to(ttprobs, simd_ns::element_aligned); ttprobs += simd_probs.size();
            }
            while(generator.advanceToNextConfiguration_no_carry())
            {
                *ttmasses = generator.mass(); ttmasses++;
                *ttprobs = generator.prob(); ttprobs++;
            }
            if(!generator.carry())
                break;
            *ttmasses = generator.mass(); ttmasses++;
            *ttprobs = generator.prob(); ttprobs++;
        } while(true);
#else
        while(generator.advanceToNextConfiguration())
        {
            *ttmasses = generator.mass(); ttmasses++;
            *ttprobs = generator.prob(); ttprobs++;
        }
#endif
    }

    // The count_confs pre-pass fixes tab_size to the exact number of above-threshold
    // configurations; the fill above must emit precisely that many. A mismatch means
    // the fill and the pre-pass disagree (e.g. a carry-convention or off-by-one bug),
    // which would leave _confs_no inconsistent with the buffer contents.
    ISOSPEC_IMPOSSIBLE(static_cast<size_t>(ttmasses - this->_masses) != tab_size);

    this->_confs_no = tab_size;
}

template void FixedEnvelope::threshold_init<true>(Iso&& iso, double threshold, bool absolute);
template void FixedEnvelope::threshold_init<false>(Iso&& iso, double threshold, bool absolute);


template<bool tgetConfs> void FixedEnvelope::total_prob_init(Iso&& iso, double target_total_prob, bool optimize)
{
    if(target_total_prob <= 0.0)
        return;

    if(target_total_prob >= 1.0)
    {
        threshold_init<tgetConfs>(std::move(iso), 0.0, true);
        return;
    }

    IsoLayeredGenerator generator(std::move(iso), 1000, 1000, true, std::min<double>(target_total_prob, 0.9999));

    this->allDim = generator.getAllDim();
    this->allDimSizeofInt = this->allDim*sizeof(int);


    this->reallocate_memory<tgetConfs>(ISOSPEC_INIT_TABLE_SIZE);

    size_t last_switch = 0;
    double prob_at_last_switch = 0.0;
    double prob_so_far = 0.0;
    double layer_delta;

    const double sum_above = log1p(-target_total_prob) - 2.3025850929940455;  // log(0.1);

    do
    {  // Store confs until we accumulate more prob than needed - and, if optimizing,
       // store also the rest of the last layer
        while(generator.advanceToNextConfigurationWithinLayer())
        {
            this->template addConfILG<tgetConfs>(generator);
            prob_so_far += *(tprobs-1);  // The just-stored probability
            if(prob_so_far >= target_total_prob)
            {
                if (optimize)
                {
                    while(generator.advanceToNextConfigurationWithinLayer())
                        this->template addConfILG<tgetConfs>(generator);
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
        layer_delta = (std::max)((std::min)(layer_delta, -0.1), -5.0);
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
        conf_swapspace = reinterpret_cast<int*>(malloc(this->allDimSizeofInt));

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
        size_t pivot = random_gen() % len + start;  // Using Mersenne twister directly - we don't
                                                    // need a very uniform distribution just for pivot
                                                    // selection
#endif
        double pprob = this->_probs[pivot];
        swap<tgetConfs>(pivot, end-1, conf_swapspace);

        double new_csum = sum_to_start;

        size_t loweridx = start;
        for(size_t ii = start; ii < end-1; ii++)
            if(this->_probs[ii] > pprob)
            {
                swap<tgetConfs>(ii, loweridx, conf_swapspace);
                new_csum += this->_probs[loweridx];
                loweridx++;
            }

        swap<tgetConfs>(end-1, loweridx, conf_swapspace);

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
        this->template reallocate_memory<tgetConfs>(end);

    this->_confs_no = end;
}

template void FixedEnvelope::total_prob_init<true>(Iso&& iso, double target_total_prob, bool optimize);
template void FixedEnvelope::total_prob_init<false>(Iso&& iso, double target_total_prob, bool optimize);

template<bool tgetConfs> void FixedEnvelope::stochastic_init(Iso&& iso, size_t _no_molecules, double _precision, double _beta_bias)
{
    IsoStochasticGenerator generator(std::move(iso), _no_molecules, _precision, _beta_bias);

    this->allDim = generator.getAllDim();
    this->allDimSizeofInt = this->allDim * sizeof(int);

    this->reallocate_memory<tgetConfs>(ISOSPEC_INIT_TABLE_SIZE);

    while(generator.advanceToNextConfiguration())
        addConfILG<tgetConfs, IsoStochasticGenerator>(generator);
}

template void FixedEnvelope::stochastic_init<true>(Iso&& iso, size_t _no_molecules, double _precision, double _beta_bias);
template void FixedEnvelope::stochastic_init<false>(Iso&& iso, size_t _no_molecules, double _precision, double _beta_bias);

double FixedEnvelope::empiric_average_mass()
{
    double ret = 0.0;
    for(size_t ii = 0; ii < _confs_no; ii++)
    {
        ret += _masses[ii] * _probs[ii];
    }
    return ret / get_total_prob();
}

double FixedEnvelope::empiric_variance()
{
    // Single sweep replacing the previous implementation, which swept the arrays
    // up to three times (empiric_average_mass(), the centered-square loop, and
    // get_total_prob()). Masses are shifted by a reference m0 (the first mass)
    // before accumulating, so the running sums stay O(spread) rather than
    // O(mass^2); this keeps the final subtraction between same-magnitude
    // quantities and avoids the catastrophic cancellation that a raw
    // Sum(m^2 p) - avg^2 would suffer for the large absolute masses seen here.
    const double m0 = (_confs_no > 0) ? _masses[0] : 0.0;
    double sum_p   = 0.0;  // Sum p
    double sum_dp  = 0.0;  // Sum p (m - m0)
    double sum_d2p = 0.0;  // Sum p (m - m0)^2
    for(size_t ii = 0; ii < _confs_no; ii++)
    {
        const double p  = _probs[ii];
        const double d  = _masses[ii] - m0;
        const double dp = d * p;
        sum_p   += p;
        sum_dp  += dp;
        sum_d2p += d * dp;
    }

    // Match the memoised get_total_prob(): reuse the cached value if present,
    // otherwise the total we just computed (and cache it, as get_total_prob would).
    if(std::isnan(total_prob))
        total_prob = sum_p;

    // avg - m0, then Var = Sum_ii p_ii ((m_ii - m0) - (avg - m0))^2 / total,
    // expanded so the loop above is the only pass over the arrays.
    const double avg_d = (sum_dp + m0 * sum_p) / total_prob - m0;
    return (sum_d2p - 2.0 * avg_d * sum_dp + avg_d * avg_d * sum_p) / total_prob;
}

// Map a mass to its bin index.  Multiplying by the precomputed reciprocal avoids
// a per-configuration division on this hot path.  The idx_min bound in Binned()
// MUST use this exact same expression: since the bin index is monotonic in mass
// and every configuration mass is >= min_mass, that guarantees bin_idx >= idx_min
// and keeps the scatter in bounds (acc is rebased by -idx_min).  (Reciprocal vs
// true division can disagree by <=1 ULP, which at an exact bin boundary could
// move a peak to an adjacent bin, but never out of the allocated range.)
static ISOSPEC_FORCE_INLINE std::ptrdiff_t bin_index(double mass, double hwmm, double inv_bin_width)
{
    return static_cast<std::ptrdiff_t>(floor((mass + hwmm) * inv_bin_width));
}

// The seeding contract for both fillers below: nonzero_idx just needs to be the
// bin of *some* populated bin, to seed the outward compaction scan in Binned().
// The scan walks both ways from the seed until it hits >10 Da of empty bins, so
// *any* nonzero bin works equally well (it need not be the smallest, leftmost, or
// the peak) -- we take whichever configuration the generator yields first, which
// is the cheapest choice and lands on the mode.  The one requirement is that the
// seed be nonzero: an empty seed lying >10 Da from the support would let the scan
// stop before ever reaching it.  Hence zero-probability configurations are skipped.

// Scatter a generator's entire output into the dense bin accumulator `acc` (used
// for the target>=1 full-enumeration path).  Returns false if empty.
template<typename GenType>
static bool fill_bins_full(GenType& generator,
                           double* acc,
                           double hwmm,
                           double inv_bin_width,
                           std::ptrdiff_t& nonzero_idx)
{
    bool non_empty;
    while((non_empty = generator.advanceToNextConfiguration()) && generator.prob() == 0.0)
    {}

    if(!non_empty)
        return false;

    nonzero_idx = bin_index(generator.mass(), hwmm, inv_bin_width);
    acc[nonzero_idx] = generator.prob();

    while(generator.advanceToNextConfiguration())
        acc[bin_index(generator.mass(), hwmm, inv_bin_width)] += generator.prob();

    return true;
}

// Scatter a layered generator's output into `acc`, stopping once target_total_prob
// is reached.  Drives the layers with a step adaptively tuned to the remaining
// probability (as total_prob_init does), rather than advanceToNextConfiguration()'s
// fixed -2.0 nat step, to cut marginal over-expansion and layer-transition
// (extend()) overhead.  Returns false if the distribution is empty.
static bool fill_bins_to_prob(IsoLayeredGenerator& generator,
                              double* acc,
                              double hwmm,
                              double inv_bin_width,
                              double target_total_prob,
                              std::ptrdiff_t& nonzero_idx)
{
    double prob_so_far = 0.0;
    double layer_delta;
    bool seeded = false;
    const double sum_above = log1p(-target_total_prob) - 2.3025850929940455;  // log(0.1)

    do
    {
        while(generator.advanceToNextConfigurationWithinLayer())
        {
            double prob = generator.prob();
            if(prob == 0.0)
                continue;

            std::ptrdiff_t bin_idx = bin_index(generator.mass(), hwmm, inv_bin_width);
            acc[bin_idx] += prob;

            if(!seeded)
            {
                nonzero_idx = bin_idx;
                seeded = true;
            }

            prob_so_far += prob;
            if(prob_so_far >= target_total_prob)
                return true;
        }

        layer_delta = sum_above - log1p(-prob_so_far);
        // Cap the widest step.  total_prob_init uses -5.0, but it trims afterwards
        // so overshoot is harmless there; Binned keeps everything it generates, so
        // a step that overshoots the target is pure waste.  -3.0 measured as a
        // uniform ~16-22% win over the old fixed -2.0 step across small and large
        // molecules with no overshoot, whereas -5.0 regressed on large ones.
        // Overridable at build time for retuning.
#ifndef ISOSPEC_BINNED_LAYER_MAXSTEP
#  define ISOSPEC_BINNED_LAYER_MAXSTEP (-3.0)
#endif
        layer_delta = (std::max)((std::min)(layer_delta, -0.1), ISOSPEC_BINNED_LAYER_MAXSTEP);
    } while(generator.nextLayer(layer_delta));

    return seeded;
}

FixedEnvelope FixedEnvelope::Binned(Iso&& iso, double target_total_prob, double bin_width, double bin_middle)
{
    FixedEnvelope ret;

    if(target_total_prob <= 0.0)
        return ret;

    double min_mass = iso.getLightestPeakMass();
    double range_len = iso.getHeaviestPeakMass() - min_mass;
    size_t no_bins = static_cast<size_t>(range_len / bin_width) + 2;
    double half_width = 0.5*bin_width;
    double hwmm = half_width-bin_middle;
    double inv_bin_width = 1.0/bin_width;
    // Bin indices are signed: for masses below bin_middle - half_width they
    // are legitimately negative.  Using size_t here would UB on the floor()
    // cast and underflow the acc-=idx_min pointer arithmetic below.
    std::ptrdiff_t idx_min = bin_index(min_mass, hwmm, inv_bin_width);
    std::ptrdiff_t idx_max = idx_min + static_cast<std::ptrdiff_t>(no_bins);

    double* acc;
# if ISOSPEC_GOT_MMAN
    acc = reinterpret_cast<double*>(mmap(nullptr, sizeof(double)*no_bins, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
#else
    // This will probably crash for large molecules and high resolutions...
    acc = reinterpret_cast<double*>(calloc(no_bins, sizeof(double)));
#endif
    if(acc == NULL)
        throw std::bad_alloc();

#if ISOSPEC_GOT_MMAN && defined(MADV_HUGEPAGE)
    // For a large accumulator the scatter (acc[bin_idx] += prob, at essentially
    // random indices) thrashes the DTLB with 4 KiB pages.  Advise transparent
    // huge pages so the touched region is backed by 2 MiB pages.  Gated by size:
    // for a small array a single on-fault 2 MiB page would just waste memory and
    // add zeroing latency with no TLB benefit.  Overridable at build time.
    // Measured (Alder Lake + Opteron 6380, THP active): ~2-9% on fine-binned large
    // molecules, no regression on coarse/small; needs a genuinely THP-capable host
    // (no effect where the kernel won't hand out anon huge pages).
#  ifndef ISOSPEC_BINNED_HUGEPAGE_MIN_BYTES
#    define ISOSPEC_BINNED_HUGEPAGE_MIN_BYTES (size_t(64) << 20)  // 64 MiB
#  endif
    if(sizeof(double)*no_bins >= ISOSPEC_BINNED_HUGEPAGE_MIN_BYTES)
        madvise(acc, sizeof(double)*no_bins, MADV_HUGEPAGE);
#endif

    acc -= idx_min;

    std::ptrdiff_t nonzero_idx = 0;
    bool non_empty;

    if(target_total_prob >= 1.0)
    {
        // The whole distribution is requested: a threshold generator with a zero
        // cutoff enumerates everything and is cheaper than the layered generator
        // descending layer-by-layer to the least-probable peak.  Mirrors the
        // >= 1.0 fast path in total_prob_init.
        IsoThresholdGenerator generator(std::move(iso), 0.0, true);
        non_empty = fill_bins_full<IsoThresholdGenerator>(
                        generator, acc, hwmm, inv_bin_width, nonzero_idx);
    }
    else
    {
        // Pass the requested probability as the layer-sizing hint (as
        // total_prob_init does) so the marginal ordering is tuned to the amount
        // of mass actually needed; the default 0.99 hint over- or under-shoots.
        IsoLayeredGenerator generator(std::move(iso), 1000, 1000, true,
                                      std::min<double>(target_total_prob, 0.9999));
        non_empty = fill_bins_to_prob(
                        generator, acc, hwmm, inv_bin_width, target_total_prob, nonzero_idx);
    }

    if(non_empty)
    {

        // Making the assumption that there won't be gaps of more than 10 Da in the spectrum. This is true for all
        // molecules made of natural elements.
        // FIXME: this has to be computed from the actual molecule, because
        // there are also people that hijack IsoSpec for statistical calculations with arbitrary distributions
        // transcribed onto artificial "elements".
        size_t distance_10da = static_cast<size_t>(10.0/bin_width) + 1;

        size_t empty_steps = 0;

        ret.reallocate_memory<false>(ISOSPEC_INIT_TABLE_SIZE);

        for(std::ptrdiff_t ii = nonzero_idx; empty_steps < distance_10da; )
        {
            if(acc[ii] > 0.0)
            {
                empty_steps = 0;
                ret.store_conf(static_cast<double>(ii)*bin_width + bin_middle, acc[ii]);
            }
            else
                empty_steps++;
            if(ii == idx_min) break;
            ii--;
        }

        empty_steps = 0;
        for(std::ptrdiff_t ii = nonzero_idx+1; ii < idx_max && empty_steps < distance_10da; ii++)
        {
            if(acc[ii] > 0.0)
            {
                empty_steps = 0;
                ret.store_conf(static_cast<double>(ii)*bin_width + bin_middle, acc[ii]);
            }
            else
                empty_steps++;
        }
    }

    acc += idx_min;

# if ISOSPEC_GOT_MMAN
    munmap(acc, sizeof(double)*no_bins);
#else
    free(acc);
#endif

    return ret;
}

}  // namespace IsoSpec
