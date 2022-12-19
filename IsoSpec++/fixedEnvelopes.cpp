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
#include "isoMath.h"

namespace IsoSpec
{

FixedEnvelope::FixedEnvelope(const FixedEnvelope& other) :
_masses(array_copy<double>(other._masses, other._confs_no)),
_probs(array_copy<double>(other._probs, other._confs_no)),
_confs(array_copy_nptr<int>(other._confs, other._confs_no*other.allDim)),
_confs_no(other._confs_no),
allDim(other.allDim),
sorted_by_mass(other.sorted_by_mass),
sorted_by_prob(other.sorted_by_prob),
total_prob(other.total_prob)
{}

FixedEnvelope::FixedEnvelope(FixedEnvelope&& other) :
_masses(other._masses),
_probs(other._probs),
_confs(other._confs),
_confs_no(other._confs_no),
allDim(other.allDim),
sorted_by_mass(other.sorted_by_mass),
sorted_by_prob(other.sorted_by_prob),
total_prob(other.total_prob)
{
other._masses = nullptr;
other._probs  = nullptr;
other._confs  = nullptr;
other._confs_no = 0;
other.total_prob = 0.0;
}

FixedEnvelope::FixedEnvelope(double* in_masses, double* in_probs, size_t in_confs_no, bool masses_sorted, bool probs_sorted, double _total_prob) :
_masses(in_masses),
_probs(in_probs),
_confs(nullptr),
_confs_no(in_confs_no),
allDim(0),
sorted_by_mass(masses_sorted),
sorted_by_prob(probs_sorted),
total_prob(_total_prob)
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

    memcpy(nprobs,  _probs,  sizeof(double) * _confs_no);
    memcpy(nmasses, _masses, sizeof(double) * _confs_no);

    memcpy(nprobs+_confs_no,  other._probs,  sizeof(double) * other._confs_no);
    memcpy(nmasses+_confs_no, other._masses, sizeof(double) * other._confs_no);

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
    if(!can_destroy)
    {
        size_t* order_c = new size_t[size];
        memcpy(order_c, order, sizeof(size_t)*size);
        order = order_c;
    }

    for(size_t ii = 0; ii < size; ii++)
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

    double pprob = 0.0;
    double cprob = 0.0;
    size_t pidx = -1; // Overflows - but it doesn't matter.

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
                accd += fabs((m-cm)*p);
                cp += p;
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

#if 0
double FixedEnvelope::ScaledAbyssalWassersteinDistance(FixedEnvelope * const * others, double abyss_depth, const double* other_scales, const size_t N)
{
    sort_by_mass();

    std::priority_queue<std::pair<double, size_t>> PQ;
    std::unique_ptr<size_t[]> indices = std::make_unique<size_t[]>(N);
    memset(indices.get(), 0, sizeof(size_t)*N);

    for(size_t ii=0; ii<N; ii++)
    {
        others[ii]->sort_by_mass();
        if(others[ii]->_confs_no > 0)
            PQ.emplace_back({others._probs[0], ii});
    }




    std::vector<std::pair<double, double>> carried;

    size_t idx_this = 0;
    size_t idx_other = 0;

    //std::cout << "AAA" << std::endl;

    auto finished = [&]() -> bool { return idx_this >= _confs_no && PQ.empty(); };
    auto next = [&]() -> std::tuple<double, double, size_t> {
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
        auto [m, p] = next();
        if(!carried.empty() && carried[0].second * p > 0.0)
        {
            carried.emplace_back(m, p);
            continue;
        }

        while(!carried.empty())
        {
            auto& [cm, cp] = carried.back();
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
                accd += fabs((m-cm)*p);
                cp += p;
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

#endif

#if 0
double AbyssalWassersteinDistanceGrad(FixedEnvelope* const* envelopes, const double* scales, double* ret_gradient, size_t N, double abyss_depth_exp, double abyss_depth_the)
{
return 0.0;
    std::unique_ptr<size_t[]> env_idx = std::make_unique<size_t[]>(N+1);
    memset(env_idx.get(), 0, (N+1)*sizeof(size_t));
    memset(ret_gradient, 0, (N+1)*sizeof(double));

    auto pq_cmp = [](std::pair<double, size_t>& p1, std::pair<double, size_t>& p2) { return p1.first > p2.first; };
    std::priority_queue<std::pair<double, size_t>, std::vector<std::pair<double, size_t>>, decltype(pq_cmp)> PQ(pq_cmp);

    for(size_t ii=0; ii<=N; ii++)
    {
        envelopes[ii]->sort_by_mass();
        if(envelopes[ii]->_confs_no > 0)
        {
            std::cout << "Initial push: " << envelopes[ii]->_masses[0] << " " << ii << '\n';
            PQ.push({envelopes[ii]->_masses[0], ii});
        }
    }

    std::vector<std::tuple<double, double, size_t>> carried;

    auto next = [&]() -> std::tuple<double, double, size_t> {
                            auto [mass, eidx] = PQ.top();
                            double prob = envelopes[eidx]->_probs[env_idx[eidx]];
                            PQ.pop();
                            if(eidx == 0)
                                prob = -prob;
                            else
                                prob = prob * scales[eidx];
                            std::cout << "Next: " << mass << " " << prob << " " << eidx << '\n';
                            env_idx[eidx]++;
                            if(env_idx[eidx] < envelopes[eidx]->_confs_no)
                                PQ.push({envelopes[eidx]->_masses[env_idx[eidx]], eidx});

                            return {mass, prob, eidx};
                        };
    double accd = 0.0;
    double condemned = 0.0;
    //double flow;
    const double max_flow_dist = abyss_depth_exp + abyss_depth_the;
    max_flow_dist *= 2.0;

    while(!PQ.empty())
    {
        auto [m, p, eidx] = next();
        if(!carried.empty() && std::get<1>(carried[0]) * p > 0.0)
        {
            carried.emplace_back(m, p, eidx);
            continue;
        }

        while(!carried.empty())
        {
            auto& [cm, cp, ceidx] = carried.back();
            if(m - cm >= max_flow_dist)
            {
                for(auto it = carried.cbegin(); it != carried.cend(); it++)
                    condemned += fabs(std::get<1>(*it));
                carried.clear();
                break;
            }
            std::cout << "accing\n";
            if((cp+p)*p > 0.0)
            {
                accd += fabs((m-cm)*cp);
                p += cp;
                carried.pop_back();
                std::cout << "cprob:" << cp << '\n';
            }
            else
            {
                accd += fabs((m-cm)*p);
                cp += p;
                std::cout << "prob:" << p << '\n';
                p = 0.0;
                break;
            }
        }
        if(p != 0.0)
            carried.emplace_back(m, p, eidx);
        //std::cout << m << " " << p << std::endl;
    }

    for(auto it = carried.cbegin(); it != carried.cend(); it++)
        condemned += fabs(std::get<1>(*it));

    std::cout << "ISO:" << accd << " " << condemned << '\n';
    return accd + condemned * max_flow_dist * 0.5;
    while(!PQ.empty())
    {
        auto [m, p, eidx] = next();
        if(!carried.empty() && (std::get<1>(carried[0]) * p > 0.0))
        {
            carried.emplace_back(m, p, eidx);
            continue;
        }

        while(!carried.empty())
        {
            auto& [cm, cp, ceidx] = carried.back();
            if(m - cm >= max_flow_dist)
            {
                for(auto it = carried.cbegin(); it != carried.cend(); it++)
                {
                    flow = fabs(std::get<1>(*it));
                    const size_t target_idx = std::get<2>(*it);
                    flow *= target_idx == 0 ? abyss_depth_exp : abyss_depth_the;
                    ret_gradient[target_idx] += flow;
                    condemned += flow;
                    std::cout << "condemnin': " << m << " " << p << " " << eidx << " | " << cm << " " << cp << " " << ceidx << '\n';
                }
                carried.clear();
                break;
            }
            if((cp+p)*p > 0.0)
            {
                flow = fabs((m-cm)*cp);
                accd += flow;
                p += cp;
                ret_gradient[ceidx] += flow;
                carried.pop_back();
            }
            else
            {
                flow = fabs((m-cm)*p);
                accd += flow;
                cp += p;
                ret_gradient[eidx] += flow;
                p = 0.0;
                break;
            }
        }
        if(p != 0.0)
            carried.emplace_back(m, p, eidx);
        //std::cout << m << " " << p << std::endl;
    }

    for(auto it = carried.cbegin(); it != carried.cend(); it++)
        condemned += fabs(std::get<1>(*it));


    return accd + condemned * (abyss_depth_exp + abyss_depth_the) * 0.5;
}
#endif


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
            if(other._masses[idx_other] < _masses[idx_this] - flow_distance)
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
    current_size = new_size;
    // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...
    _masses = reinterpret_cast<double*>(realloc(_masses, new_size * sizeof(double)));
    if(_masses == nullptr)
        throw std::bad_alloc();
    tmasses = _masses + _confs_no;

    _probs  = reinterpret_cast<double*>(realloc(_probs,  new_size * sizeof(double)));
    if(_probs == nullptr)
        throw std::bad_alloc();
    tprobs  = _probs  + _confs_no;

    constexpr_if(tgetConfs)
    {
        _confs  = reinterpret_cast<int*>(realloc(_confs,  new_size * allDimSizeofInt));
        if(_confs == nullptr)
            throw std::bad_alloc();
        tconfs = _confs + (allDim * _confs_no);
    }
}

void FixedEnvelope::slow_reallocate_memory(size_t new_size)
{
    current_size = new_size;
    // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...
    _masses = reinterpret_cast<double*>(realloc(_masses, new_size * sizeof(double)));
    if(_masses == nullptr)
        throw std::bad_alloc();
    tmasses = _masses + _confs_no;

    _probs  = reinterpret_cast<double*>(realloc(_probs,  new_size * sizeof(double)));
    if(_probs == nullptr)
        throw std::bad_alloc();
    tprobs  = _probs  + _confs_no;

    if(_confs != nullptr)
    {
        _confs  = reinterpret_cast<int*>(realloc(_confs,  new_size * allDimSizeofInt));
        if(_confs == nullptr)
            throw std::bad_alloc();
        tconfs = _confs + (allDim * _confs_no);
    }
}

template<bool tgetConfs> void FixedEnvelope::threshold_init(Iso&& iso, double threshold, bool absolute)
{
    IsoThresholdGenerator generator(std::move(iso), threshold, absolute);

    size_t tab_size = generator.count_confs();
    this->allDim = generator.getAllDim();
    this->allDimSizeofInt = this->allDim * sizeof(int);

    this->reallocate_memory<tgetConfs>(tab_size);

    double* ttmasses = this->_masses;
    double* ttprobs = this->_probs;
    ISOSPEC_MAYBE_UNUSED int* ttconfs;
    constexpr_if(tgetConfs)
        ttconfs = _confs;

    while(generator.advanceToNextConfiguration())
    {
        *ttmasses = generator.mass(); ttmasses++;
        *ttprobs = generator.prob(); ttprobs++;
        constexpr_if(tgetConfs)  { generator.get_conf_signature(ttconfs); ttconfs += allDim; }
    }

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
    double ret = 0.0;
    double avg = empiric_average_mass();
    for(size_t ii = 0; ii < _confs_no; ii++)
    {
        double msq = _masses[ii] - avg;
        ret += msq * msq * _probs[ii];
    }

    return ret / get_total_prob();
}

FixedEnvelope FixedEnvelope::Binned(Iso&& iso, double target_total_prob, double bin_width, double bin_middle)
{
    FixedEnvelope ret;

    double min_mass = iso.getLightestPeakMass();
    double range_len = iso.getHeaviestPeakMass() - min_mass;
    size_t no_bins = static_cast<size_t>(range_len / bin_width) + 2;
    double half_width = 0.5*bin_width;
    double hwmm = half_width-bin_middle;
    size_t idx_min = floor((min_mass + hwmm) / bin_width);
    size_t idx_max = idx_min + no_bins;

    double* acc;
# if ISOSPEC_GOT_MMAN
    acc = reinterpret_cast<double*>(mmap(nullptr, sizeof(double)*no_bins, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
#else
    // This will probably crash for large molecules and high resolutions...
    acc = reinterpret_cast<double*>(calloc(no_bins, sizeof(double)));
#endif
    if(acc == NULL)
        throw std::bad_alloc();

    acc -= idx_min;

    IsoLayeredGenerator ITG(std::move(iso));


    bool non_empty;
    while((non_empty = ITG.advanceToNextConfiguration()) && ITG.prob() == 0.0)
    {}

    if(non_empty)
    {
        double accum_prob = ITG.prob();
        size_t nonzero_idx = floor((ITG.mass() + hwmm)/bin_width);
        acc[nonzero_idx] = accum_prob;

        while(ITG.advanceToNextConfiguration() && accum_prob < target_total_prob)
        {
            double prob = ITG.prob();
            accum_prob += prob;
            size_t bin_idx = floor((ITG.mass() + hwmm)/bin_width);
            acc[bin_idx] += prob;
        }

        // Making the assumption that there won't be gaps of more than 10 Da in the spectrum. This is true for all
        // molecules made of natural elements.
        size_t distance_10da = static_cast<size_t>(10.0/bin_width) + 1;

        size_t empty_steps = 0;

        ret.reallocate_memory<false>(ISOSPEC_INIT_TABLE_SIZE);

        for(size_t ii = nonzero_idx; ii >= idx_min && empty_steps < distance_10da; ii--)
        {
            if(acc[ii] > 0.0)
            {
                empty_steps = 0;
                ret.store_conf(static_cast<double>(ii)*bin_width + bin_middle, acc[ii]);
            }
            else
                empty_steps++;
        }

        empty_steps = 0;
        for(size_t ii = nonzero_idx+1; ii < idx_max && empty_steps < distance_10da; ii++)
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
