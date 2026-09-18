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
#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>
#include "isoMath.h"
#include "aligned_ptr.h"

namespace IsoSpec
{

// Copies of an envelope are always made into aligned arrays of our own, even
// when the source adopted foreign, arbitrarily aligned ones. The size is known
// up front, so these never grow and stay on the small backend.
FixedEnvelope::FixedEnvelope(const FixedEnvelope& other) :
_masses(nullptr),
_probs(nullptr),
_confs(nullptr),
_confs_no(other._confs_no),
allDim(other.allDim),
sorted_by_mass(other.sorted_by_mass),
sorted_by_prob(other.sorted_by_prob),
total_prob(other.total_prob),
current_size(other._confs_no),
allDimSizeofInt(other.allDimSizeofInt)
{
    if(other._masses != nullptr)
    {
        _masses_owner.reset(_confs_no);
        _masses = _masses_owner.get();
        if(_confs_no > 0)
            memcpy(_masses, other._masses, sizeof(double) * _confs_no);
    }

    if(other._probs != nullptr)
    {
        _probs_owner.reset(_confs_no);
        _probs = _probs_owner.get();
        if(_confs_no > 0)
            memcpy(_probs, other._probs, sizeof(double) * _confs_no);
    }

    if(other._confs != nullptr)
    {
        const size_t no_ints = _confs_no * static_cast<size_t>(allDim);
        _confs_owner.reset(no_ints);
        _confs = _confs_owner.get();
        if(no_ints > 0)
            memcpy(_confs, other._confs, sizeof(int) * no_ints);
    }
}

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
allDimSizeofInt(other.allDimSizeofInt),
_masses_owner(std::move(other._masses_owner)),
_probs_owner(std::move(other._probs_owner)),
_confs_owner(std::move(other._confs_owner))
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
    std::swap(_masses_owner,   tmp._masses_owner);
    std::swap(_probs_owner,    tmp._probs_owner);
    std::swap(_confs_owner,    tmp._confs_owner);
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
    free_array(_masses_owner, _masses);
    free_array(_probs_owner,  _probs);
    free_array(_confs_owner,  _confs);
    _masses_owner   = std::move(other._masses_owner);
    _probs_owner    = std::move(other._probs_owner);
    _confs_owner    = std::move(other._confs_owner);
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
current_size(in_confs_no),
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
current_size(in_confs_no),
allDimSizeofInt(_allDim * sizeof(int))
{}

FixedEnvelope::FixedEnvelope(aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT>&& in_masses,
                             aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT>&& in_probs,
                             size_t in_confs_no, bool masses_sorted, bool probs_sorted, double _total_prob) :
_masses(in_masses.get()),
_probs(in_probs.get()),
_confs(nullptr),
_confs_no(in_confs_no),
allDim(0),
sorted_by_mass(masses_sorted),
sorted_by_prob(probs_sorted),
total_prob(_total_prob),
current_size(in_confs_no),
allDimSizeofInt(0),
_masses_owner(std::move(in_masses)),
_probs_owner(std::move(in_probs))
{}

FixedEnvelope FixedEnvelope::operator+(const FixedEnvelope& other) const
{
    const size_t new_size = _confs_no + other._confs_no;

    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> nprobs;
    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> nmasses;
    nprobs.reset(new_size);
    nmasses.reset(new_size);

    // An empty envelope has null arrays, and memcpy() is declared nonnull: a
    // zero-length copy from it is undefined behaviour, so skip it entirely.
    if(_confs_no > 0)
    {
        memcpy(nprobs.get(),  _probs,  sizeof(double) * _confs_no);
        memcpy(nmasses.get(), _masses, sizeof(double) * _confs_no);
    }

    if(other._confs_no > 0)
    {
        memcpy(nprobs.get()+_confs_no,  other._probs,  sizeof(double) * other._confs_no);
        memcpy(nmasses.get()+_confs_no, other._masses, sizeof(double) * other._confs_no);
    }

    return FixedEnvelope(std::move(nmasses), std::move(nprobs), new_size);
}

FixedEnvelope FixedEnvelope::operator*(const FixedEnvelope& other) const
{
    if(other._confs_no != 0 && _confs_no > SIZE_MAX / other._confs_no)
        throw std::bad_alloc();

    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> nprobs;
    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> nmasses;
    nprobs.reset(_confs_no * other._confs_no);
    nmasses.reset(_confs_no * other._confs_no);

    size_t tgt_idx = 0;

    for(size_t ii = 0; ii < _confs_no; ii++)
        for(size_t jj = 0; jj < other._confs_no; jj++)
        {
            nprobs[tgt_idx]  = _probs[ii]  * other._probs[jj];
            nmasses[tgt_idx] = _masses[ii] + other._masses[jj];
            tgt_idx++;
        }

    return FixedEnvelope(std::move(nmasses), std::move(nprobs), tgt_idx);
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
    size_t pidx = -1;  // Overflows - but it doesn't matter.

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

    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> newprobs;
    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> newmasses;
    newprobs.reset(ret_size);
    newmasses.reset(ret_size);

    size_t cntr = 0;
    for(size_t ii = 0; ii < size; ii++)
    {
        double mul = intensities[ii];
        for(size_t jj = 0; jj < spectra[ii]->_confs_no; jj++)
            newprobs[jj+cntr] = spectra[ii]->_probs[jj] * mul;
        if(spectra[ii]->_confs_no > 0)
            memcpy(newmasses.get() + cntr, spectra[ii]->_masses, sizeof(double) * spectra[ii]->_confs_no);
        cntr += spectra[ii]->_confs_no;
    }
    return FixedEnvelope(std::move(newmasses), std::move(newprobs), cntr);
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
        for(size_t ii = 1; ii < _confs_no; ii++)
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
    // aligned_unique_ptr::realloc() does the overflow-checked size arithmetic
    // and throws std::bad_alloc on failure, so there is nothing to check here.
    grow_array(_masses_owner, _masses, new_size, current_size);
    tmasses = _masses + _confs_no;

    grow_array(_probs_owner, _probs, new_size, current_size);
    tprobs  = _probs  + _confs_no;

    constexpr_if(tgetConfs)
    {
        if(allDim > 0 && new_size > SIZE_MAX / static_cast<size_t>(allDim))
            throw std::bad_alloc();
        const size_t new_ints = new_size * static_cast<size_t>(allDim);
        grow_array(_confs_owner, _confs, new_ints, current_size * static_cast<size_t>(allDim));
        tconfs = _confs + (allDim * _confs_no);
    }
    current_size = new_size;
}

template<bool tgetConfs> void FixedEnvelope::allocate_memory(size_t new_size)
{
    fresh_array(_masses_owner, _masses, new_size);
    tmasses = _masses + _confs_no;

    fresh_array(_probs_owner, _probs, new_size);
    tprobs  = _probs  + _confs_no;

    constexpr_if(tgetConfs)
    {
        if(allDim > 0 && new_size > SIZE_MAX / static_cast<size_t>(allDim))
            throw std::bad_alloc();
        fresh_array(_confs_owner, _confs, new_size * static_cast<size_t>(allDim));
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

    // count_confs() has already fixed the exact output size, so this allocates
    // it once and never grows.
    this->allocate_memory<tgetConfs>(tab_size);

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
        // Batched fill. Each marginal-0 run (one per higher-dimensional carry state)
        // starts at index 0 and descends until it drops below the cutoff; the kernel
        // takes the bulk of a run W at a time and the < W tail is drained scalar.
        // Which build of the kernel runs is decided from the CPU, once -- see
        // isa_kernels.h. When the toolchain has no <experimental/simd> at all the
        // kernel emits nothing and the scalar drain below does the whole job.
        //
        // Convention bridge: the kernel and advanceToNextConfiguration_no_carry()
        // both use "lProbs_ptr points at the last-emitted config" (advance, then
        // emit), matching the generator's initial position of one-before-index-0.
        // carry() instead leaves lProbs_ptr *at* the index-0 config to be emitted.
        // So after each successful carry we emit that index-0 config scalar (it is
        // always above the cutoff when carry() succeeds, exactly as the scalar path
        // relies on), which restores the last-emitted convention and lets the next
        // batch resume cleanly from index 1.
        const isa_fill_run_fn fill_run = isa_kernels().fill_run;
        do {
            const size_t batched = generator.batch_fill_run(fill_run, ttmasses, ttprobs);
            ttmasses += batched;
            ttprobs  += batched;

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


// Left deliberately unbatched, unlike threshold_init()'s fill: the layered
// generator's runs over marginal 0 are far too short for it. A layer admits
// only configurations whose log-probability falls in a thin band, and for any
// fixed choice of the higher marginals very few of marginal 0's fall in it --
// measured at 1.36 configurations per run (37889 configurations over 27863 runs,
// C520H817N139O147S8 to p=0.99999), where a W=4 batch needs 5. Batching it
// anyway fires on 0.0% of configurations and costs 1-4% in the probe it adds per
// run. Widening the layers does not rescue it either: at a -50 nat cap (vs the
// -5 of ISOSPEC_TOTALPROB_LAYER_MAXSTEP) the batch still covers only 2.8% of
// configurations, while generating 13% more of them.

template<bool tgetConfs> void FixedEnvelope::store_layer(IsoLayeredGenerator& generator)
{
    while(generator.advanceToNextConfigurationWithinLayer())
        this->template addConfILG<tgetConfs>(generator);
}

template<bool tgetConfs> bool FixedEnvelope::store_layer_to_prob(IsoLayeredGenerator& generator,
                                                                 double& prob_so_far,
                                                                 double target_total_prob)
{
    while(generator.advanceToNextConfigurationWithinLayer())
    {
        this->template addConfILG<tgetConfs>(generator);
        prob_so_far += *(tprobs-1);  // The just-stored probability
        if(prob_so_far >= target_total_prob)
            return true;
    }
    return false;
}

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
        if(this->template store_layer_to_prob<tgetConfs>(generator, prob_so_far, target_total_prob))
        {
            if(!optimize)
                return;
            this->template store_layer<tgetConfs>(generator);
            break;
        }

        last_switch = this->_confs_no;
        prob_at_last_switch = prob_so_far;

        layer_delta = sum_above - log1p(-prob_so_far);
        layer_delta = (std::max)((std::min)(layer_delta, -0.1), ISOSPEC_TOTALPROB_LAYER_MAXSTEP);
    } while(generator.nextLayer(layer_delta));

    if(!optimize || prob_so_far <= target_total_prob)
        return;

    // Right. We have extra configurations and we have been asked to produce an optimal p-set, so
    // now we shall trim unneeded configurations, using an algorithm dubbed "quicktrim"
    // - similar to the quickselect algorithm, except that we use the cumulative sum of elements
    // left of pivot to decide whether to go left or right, instead of the positional index.
    // We'll be sorting by the prob array, permuting the other ones in parallel.
    //
    // The partition is Hoare-style, with the left side's sum accumulated inside
    // the partition scans themselves. This pass is memory-bound, so what counts
    // is passes over the data and stores into it: the fused sum keeps it at one
    // pass per level (a separate summing pass gives back most of the win), and
    // Hoare swaps only misplaced *pairs* where the previous Lomuto partition
    // swapped once per every element above the pivot. Measured on the harvested
    // real trim inputs: 1.65x faster than Lomuto on Piledriver, 1.3x on Alder
    // Lake, never slower, including L2-resident segment sizes.

    int* conf_swapspace = nullptr;
    constexpr_if(tgetConfs)
        conf_swapspace = reinterpret_cast<int*>(malloc(this->allDimSizeofInt));

    size_t start = last_switch;
    size_t end = this->_confs_no;
    double sum_to_start = prob_at_last_switch;

    while(start < end)
    {
        size_t len = end - start;
        if(len == 1)
        {
            // A single element cannot be partitioned (both Hoare sides must be
            // nonempty): it is kept iff the prefix sum still needs it.
            if(sum_to_start >= target_total_prob)
                end = start;
            break;
        }
#if ISOSPEC_BUILDING_R
        size_t pivot = len/2 + start;
#else
        size_t pivot = random_gen() % len + start;  // Using Mersenne twister directly - we don't
                                                    // need a very uniform distribution just for pivot
                                                    // selection
#endif
        const double pprob = this->_probs[pivot];

        // Descending Hoare partition of [start, end) around the value pprob:
        // afterwards [start, split) all >= pprob and [split, end) all <= pprob.
        // left_sum accumulates the values that end up on the left: everything
        // the i-scan walks past stays there, and each swap moves the j-side
        // value into a left slot.
        std::ptrdiff_t ii = static_cast<std::ptrdiff_t>(start) - 1;
        std::ptrdiff_t jj = static_cast<std::ptrdiff_t>(end);
        double left_sum = 0.0;
        while(true)
        {
            while(true)
            {
                ii++;
                if(!(this->_probs[ii] > pprob))
                    break;
                left_sum += this->_probs[ii];
            }
            do { jj--; } while(this->_probs[jj] < pprob);
            if(ii >= jj)
                break;
            swap<tgetConfs>(ii, jj, conf_swapspace);
            left_sum += this->_probs[ii];
        }
        // The scans cross with at most one element between them (it then equals
        // the pivot). jj == ii is the case where that element sits on the left
        // without the i-scan ever having counted it.
        if(jj == ii)
            left_sum += this->_probs[ii];

        size_t split = static_cast<size_t>(jj) + 1;
        if(split == end)
        {
            // Right side came out empty; then _probs[end-1] == pprob (the j-scan
            // stopped on it immediately), so splitting just before it is a valid
            // partition too, and keeps both sides strictly shrinking.
            split = end - 1;
            left_sum -= this->_probs[end-1];
        }

        // Selection part
        double new_csum = sum_to_start + left_sum;
        if(new_csum < target_total_prob)
        {
            sum_to_start = new_csum;
            start = split;
        }
        else
            end = split;
    }

    constexpr_if(tgetConfs)
        free(conf_swapspace);

    // Set before the shrink below, so that the tmasses/tprobs/tconfs cursors it
    // recomputes land inside the smaller buffer rather than past its end.
    this->_confs_no = end;

    if(end <= current_size/2)
        // Overhead in memory of 2x or more, shrink to fit
        this->template reallocate_memory<tgetConfs>(end);
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
static bool fill_bins_full(IsoThresholdGenerator& generator,
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

    // Same batching as threshold_init()'s fill, including the carry convention:
    // the seeding loop above left the pointer on the configuration it emitted,
    // and IsoThresholdGenerator::carry() likewise leaves it on the index-0
    // configuration of the next run, so that one is scattered by hand before
    // batching resumes at index 1.
    const isa_bin_run_fn bin_run = isa_kernels().bin_run;
    do
    {
        generator.batch_bin_run(bin_run, acc, hwmm, inv_bin_width);
        while(generator.advanceToNextConfiguration_no_carry())
            acc[bin_index(generator.mass(), hwmm, inv_bin_width)] += generator.prob();
        if(!generator.carry())
            break;
        acc[bin_index(generator.mass(), hwmm, inv_bin_width)] += generator.prob();
    } while(true);

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
        // Unbatched, deliberately: see the note on store_layer() -- the layered
        // generator's runs over marginal 0 are ~1.4 configurations long, so the
        // batched kernel never fires here. (Binned()'s other filler, over the
        // threshold generator, does batch and does gain from it.)
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
        // Cap the widest step; see ISOSPEC_BINNED_LAYER_MAXSTEP in fixedEnvelopes.h.
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
        non_empty = fill_bins_full(
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
