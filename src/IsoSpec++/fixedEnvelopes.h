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

#pragma once

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <utility>

#include "isoSpec++.h"

#ifdef DEBUG
#define ISOSPEC_INIT_TABLE_SIZE 16
#else
#define ISOSPEC_INIT_TABLE_SIZE 1024
#endif

namespace IsoSpec
{

class FixedEnvelope;

double AbyssalWassersteinDistanceGrad(FixedEnvelope* const* envelopes, const double* scales, double* ret_gradient, size_t N, double abyss_depth_exp, double abyss_depth_the);


class ISOSPEC_EXPORT_SYMBOL FixedEnvelope {
 protected:
    double* _masses;
    double* _probs;
    int*    _confs;
    size_t  _confs_no;
    int     allDim;
    bool sorted_by_mass;
    bool sorted_by_prob;
    double total_prob;
    size_t current_size;
    double* tmasses;
    double* tprobs;
    int*    tconfs;
    int allDimSizeofInt;

 public:
    ISOSPEC_FORCE_INLINE FixedEnvelope() : _masses(nullptr),
        _probs(nullptr),
        _confs(nullptr),
        _confs_no(0),
        allDim(0),
        sorted_by_mass(false),
        sorted_by_prob(false),
        total_prob(0.0),
        current_size(0),
        allDimSizeofInt(0)
        // Deliberately not initializing tmasses, tprobs, tconfs
        {};

    FixedEnvelope(const FixedEnvelope& other);
    FixedEnvelope(FixedEnvelope&& other);

    FixedEnvelope(double* masses, double* probs, size_t confs_no, bool masses_sorted = false, bool probs_sorted = false, double _total_prob = NAN);

    virtual ~FixedEnvelope()
    {
        free(_masses);
        free(_probs);
        free(_confs);
    };

    FixedEnvelope operator+(const FixedEnvelope& other) const;
    FixedEnvelope operator*(const FixedEnvelope& other) const;

    inline size_t    confs_no()  const { return _confs_no; }
    inline int       getAllDim() const { return allDim; }

    inline const double*   masses() const { return _masses; }
    inline const double*   probs()  const { return _probs; }
    inline const int*      confs()  const { return _confs; }

    inline double*   release_masses()     { double* ret = _masses; _masses = nullptr; return ret; }
    inline double*   release_probs()      { double* ret = _probs;  _probs  = nullptr; return ret; }
    inline int*      release_confs()      { int*    ret = _confs;  _confs  = nullptr; return ret; }
    inline void      release_everything() { _confs = nullptr; _probs = _masses = nullptr; }


    inline double     mass(size_t i)  const { return _masses[i]; }
    inline double     prob(size_t i)  const { return _probs[i];  }
    inline const int* conf(size_t i)  const { return _confs + i*allDim; }

    void sort_by_mass();
    void sort_by_prob();

    double get_total_prob();
    void scale(double factor);
    void normalize();
    void shift_mass(double shift);
    void resample(size_t ionic_current, double beta_bias = 1.0);

    double empiric_average_mass();
    double empiric_variance();
    double empiric_stddev() { return sqrt(empiric_variance()); }

    double WassersteinDistance(FixedEnvelope& other);
    double OrientedWassersteinDistance(FixedEnvelope& other);
    double AbyssalWassersteinDistance(FixedEnvelope& other, double abyss_depth, double other_scale = 1.0);
    std::tuple<double, double, double> WassersteinMatch(FixedEnvelope& other, double flow_distance, double other_scale = 1.0);


    static FixedEnvelope LinearCombination(const std::vector<const FixedEnvelope*>& spectra, const std::vector<double>& intensities);
    static FixedEnvelope LinearCombination(const FixedEnvelope* const * spectra, const double* intensities, size_t size);


    FixedEnvelope bin(double bin_width = 1.0, double middle = 0.0);

 private:
    void sort_by(double* order);


 protected:
    template<typename T, bool tgetConfs> ISOSPEC_FORCE_INLINE void store_conf(const T& generator)
    {
        *tmasses = generator.mass(); tmasses++;
        *tprobs  = generator.prob(); tprobs++;
        constexpr_if(tgetConfs) { generator.get_conf_signature(tconfs); tconfs += allDim; }
    }

    ISOSPEC_FORCE_INLINE void store_conf(double _mass, double _prob)
    {
        if(_confs_no == current_size)
            reallocate_memory<false>(current_size*2);

        *tprobs = _prob;
        *tmasses = _mass;
        tprobs++;
        tmasses++;

        _confs_no++;
    }

    template<bool tgetConfs> ISOSPEC_FORCE_INLINE void swap(size_t idx1, size_t idx2, ISOSPEC_MAYBE_UNUSED int* conf_swapspace)
    {
        std::swap<double>(this->_probs[idx1],  this->_probs[idx2]);
        std::swap<double>(this->_masses[idx1], this->_masses[idx2]);
        constexpr_if(tgetConfs)
        {
            int* c1 = this->_confs + (idx1*this->allDim);
            int* c2 = this->_confs + (idx2*this->allDim);
            memcpy(conf_swapspace, c1, this->allDimSizeofInt);
            memcpy(c1, c2, this->allDimSizeofInt);
            memcpy(c2, conf_swapspace, this->allDimSizeofInt);
        }
    }

    template<bool tgetConfs> void reallocate_memory(size_t new_size);
    void slow_reallocate_memory(size_t new_size);

 public:
    template<bool tgetConfs> void threshold_init(Iso&& iso, double threshold, bool absolute);

    template<bool tgetConfs, typename GenType = IsoLayeredGenerator> void addConfILG(const GenType& generator)
    {
        if(this->_confs_no == this->current_size)
            this->template reallocate_memory<tgetConfs>(this->current_size*2);

        this->template store_conf<GenType, tgetConfs>(generator);
        this->_confs_no++;
    }

    template<bool tgetConfs> void total_prob_init(Iso&& iso, double target_prob, bool trim);

    static FixedEnvelope FromThreshold(Iso&& iso, double threshold, bool absolute, bool tgetConfs = false)
    {
        FixedEnvelope ret;

        if(tgetConfs)
            ret.threshold_init<true>(std::move(iso), threshold, absolute);
        else
            ret.threshold_init<false>(std::move(iso), threshold, absolute);
        return ret;
    }

    inline static FixedEnvelope FromThreshold(const Iso& iso, double _threshold, bool _absolute, bool tgetConfs = false)
    {
        return FromThreshold(Iso(iso, false), _threshold, _absolute, tgetConfs);
    }

    static FixedEnvelope FromTotalProb(Iso&& iso, double target_total_prob, bool optimize, bool tgetConfs = false)
    {
        FixedEnvelope ret;

        if(tgetConfs)
            ret.total_prob_init<true>(std::move(iso), target_total_prob, optimize);
        else
            ret.total_prob_init<false>(std::move(iso), target_total_prob, optimize);

        return ret;
    }

    inline static FixedEnvelope FromTotalProb(const Iso& iso, double _target_total_prob, bool _optimize, bool tgetConfs = false)
    {
        return FromTotalProb(Iso(iso, false), _target_total_prob, _optimize, tgetConfs);
    }

    template<bool tgetConfs> void stochastic_init(Iso&& iso, size_t _no_molecules, double _precision, double _beta_bias);

    inline static FixedEnvelope FromStochastic(Iso&& iso, size_t _no_molecules, double _precision = 0.9999, double _beta_bias = 5.0, bool tgetConfs = false)
    {
        FixedEnvelope ret;

        if(tgetConfs)
            ret.stochastic_init<true>(std::move(iso), _no_molecules, _precision, _beta_bias);
        else
            ret.stochastic_init<false>(std::move(iso), _no_molecules, _precision, _beta_bias);

        return ret;
    }

    static FixedEnvelope FromStochastic(const Iso& iso, size_t _no_molecules, double _precision = 0.9999, double _beta_bias = 5.0, bool tgetConfs = false)
    {
        return FromStochastic(Iso(iso, false), _no_molecules, _precision, _beta_bias, tgetConfs);
    }

    static FixedEnvelope Binned(Iso&& iso, double target_total_prob, double bin_width, double bin_middle = 0.0);
    static FixedEnvelope Binned(const Iso& iso, double target_total_prob, double bin_width, double bin_middle = 0.0)
    {
        return Binned(Iso(iso, false), target_total_prob, bin_width, bin_middle);
    }

    friend double AbyssalWassersteinDistanceGrad(FixedEnvelope* const* envelopes, const double* scales, double* ret_gradient, size_t N, double abyss_depth_exp, double abyss_depth_the);
};

}  // namespace IsoSpec
