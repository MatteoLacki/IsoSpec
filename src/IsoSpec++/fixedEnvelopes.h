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
#include <cstring>
#include <algorithm>
#include <vector>
#include <utility>
#include <tuple>

#include "aligned_ptr.h"
#include "isoSpec++.h"

#ifdef DEBUG
#define ISOSPEC_INIT_TABLE_SIZE 16
#else
#define ISOSPEC_INIT_TABLE_SIZE 1024
#endif

// Widest layer step Binned() will take.  total_prob_init uses -5.0, but it trims
// afterwards so overshoot is harmless there; Binned keeps everything it
// generates, so a step that overshoots the target is pure waste.  -3.0 measured
// as a uniform ~16-22% win over the fixed -2.0 step advanceToNextConfiguration()
// would take, across small and large molecules with no overshoot, whereas -5.0
// regressed on large ones.  Overridable at build time for retuning; declared
// here rather than in the .cpp so the tests can mirror the layer stepping.
#ifndef ISOSPEC_BINNED_LAYER_MAXSTEP
#  define ISOSPEC_BINNED_LAYER_MAXSTEP (-3.0)
#endif

// The same knob for total_prob_init, which can afford a wider step than Binned:
// it trims afterwards, so overshooting the target costs only the work of
// generating what gets thrown away.
#ifndef ISOSPEC_TOTALPROB_LAYER_MAXSTEP
#  define ISOSPEC_TOTALPROB_LAYER_MAXSTEP (-5.0)
#endif

namespace IsoSpec
{

class FixedEnvelope;

//! One of FixedEnvelope's arrays, handed over to the caller together with the
//! means of freeing it.
/*!
    release_masses() and friends promise a pointer that plain free() accepts,
    which costs a copy whenever the array is not actually a malloc()'d one (see
    aligned_unique_ptr::release()). These are the copy-free counterpart: the
    caller must invoke deleter(ptr, size) exactly once, and nothing else --
    in particular not free(), unless the deleter happens to be that.
*/
struct ISOSPEC_EXPORT_SYMBOL ReleasedArray {
    void*       ptr;                              /*!< nullptr if the envelope had no such array. */
    std::size_t size;                             /*!< Byte count the deleter expects; meaningless on its own. */
    void (*deleter)(void*, std::size_t) noexcept; /*!< Never null when ptr isn't. */
};

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

    // Ownership records for the three arrays above. The arrays themselves are
    // always reached through the plain pointers, so every hot loop indexes them
    // directly with no branch; these only say how each one has to be given back.
    // An array comes from exactly one of two places, and the two free differently:
    //   * allocated by us -- the matching aligned_unique_ptr is engaged and holds
    //     it. Always DOUBLE_SIMD_ALIGNMENT-aligned, so the SIMD fills can store
    //     to it; a buffer that grows repeatedly also becomes VM-backed past
    //     ISOSPEC_ALIGNED_PTR_VM_THRESHOLD, so growing it stops copying.
    //   * adopted across the C ABI -- the aligned_unique_ptr is empty and the
    //     plain pointer is a malloc()'d buffer somebody else built (a cffi array
    //     taken over zero-copy). Arbitrarily aligned; goes back to free().
    // Invariant: the owner is engaged iff the pointer beside it is ours, and
    // then the two are equal.
    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> _masses_owner;
    aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> _probs_owner;
    aligned_unique_ptr<int,    DOUBLE_SIMD_ALIGNMENT> _confs_owner;

    static void free_deleter(void* p, std::size_t) noexcept { free(p); }

    //! Free whichever of the two backends holds ptr, and null it.
    template<typename T, std::size_t A>
    static void free_array(aligned_unique_ptr<T, A>& owner, T*& ptr) noexcept
    {
        if(owner)
            owner.reset(0);
        else
            free(ptr);
        ptr = nullptr;
    }

    //! Hand ptr back as something free() accepts, copying only if it has to.
    template<typename T, std::size_t A>
    static T* release_array(aligned_unique_ptr<T, A>& owner, T*& ptr)
    {
        T* ret = owner ? owner.release() : ptr;
        ptr = nullptr;
        return ret;
    }

    //! Hand ptr back with its deleter. Never copies.
    template<typename T, std::size_t A>
    static ReleasedArray release_array_with_deleter(aligned_unique_ptr<T, A>& owner, T*& ptr) noexcept
    {
        ReleasedArray ret;
        if(owner)
        {
            typename aligned_unique_ptr<T, A>::release_result r = owner.release_with_deleter();
            ret = ReleasedArray{static_cast<void*>(r.ptr), r.size, r.deleter};
        }
        else
            ret = ReleasedArray{static_cast<void*>(ptr), 0, &free_deleter};
        ptr = nullptr;
        return ret;
    }

    //! Forget ptr without freeing it -- whoever asked for this now owns it.
    template<typename T, std::size_t A>
    static void disown_array(aligned_unique_ptr<T, A>& owner, T*& ptr) noexcept
    {
        if(owner)
            (void)owner.release_with_deleter();  // discards the deleter, but never copies
        ptr = nullptr;
    }

    //! Allocate new_n elements, keeping nothing. See ISOSPEC_ALIGNED_PTR_VM_THRESHOLD
    //! for why a one-shot allocation of a known size is not the same operation
    //! as growing into one.
    template<typename T, std::size_t A>
    static void fresh_array(aligned_unique_ptr<T, A>& owner, T*& ptr, size_t new_n)
    {
        if(!owner && ptr != nullptr)
        {
            free(ptr);
            ptr = nullptr;
        }
        owner.reset(new_n);
        ptr = owner.get();
    }

    //! Resize to new_n elements, preserving min(old_n, new_n) of them.
    template<typename T, std::size_t A>
    static void grow_array(aligned_unique_ptr<T, A>& owner, T*& ptr, size_t new_n, size_t old_n)
    {
        if(!owner && ptr != nullptr)
        {
            // An adopted buffer did not come from the aligned allocator and so
            // cannot be handed to it for resizing: move the contents across into
            // one that did, then give the original back to free(). Only reachable
            // if somebody grows an envelope built from foreign arrays.
            aligned_unique_ptr<T, A> fresh;
            fresh.realloc(new_n);
            const size_t keep = (std::min)(old_n, new_n);
            if(keep > 0)
                memcpy(fresh.get(), ptr, sizeof(T) * keep);
            free(ptr);
            owner = std::move(fresh);
        }
        else
            // Self-allocated: for a VM-backed buffer this is an in-place remap
            // rather than a copy.
            owner.realloc(new_n);
        ptr = owner.get();
    }

 public:
    ISOSPEC_FORCE_INLINE FixedEnvelope() : _masses(nullptr),
        _probs(nullptr),
        _confs(nullptr),
        _confs_no(0),
        allDim(0),
        sorted_by_mass(false),
        sorted_by_prob(false),
        total_prob(NAN),
        current_size(0),
        allDimSizeofInt(0)
        // Deliberately not initializing tmasses, tprobs, tconfs
        {};

    FixedEnvelope(const FixedEnvelope& other);
    FixedEnvelope(FixedEnvelope&& other);
    FixedEnvelope& operator=(const FixedEnvelope& other);
    FixedEnvelope& operator=(FixedEnvelope&& other);

    //! Adopt foreign, malloc()'d arrays (see the ownership note above).
    FixedEnvelope(double* masses, double* probs, size_t confs_no, bool masses_sorted = false, bool probs_sorted = false, double _total_prob = NAN);
    FixedEnvelope(double* masses, double* probs, int* confs, size_t confs_no, int _allDim, bool masses_sorted = false, bool probs_sorted = false, double _total_prob = NAN);

    //! Internal counterpart to the two adopting constructors above: takes over
    //! arrays this library allocated itself, so that everything the library
    //! produces stays aligned instead of being demoted to a plain malloc()'d,
    //! foreign buffer.
    FixedEnvelope(aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT>&& masses,
                  aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT>&& probs,
                  size_t confs_no, bool masses_sorted = false, bool probs_sorted = false, double _total_prob = NAN);

    virtual ~FixedEnvelope()
    {
        free_array(_masses_owner, _masses);
        free_array(_probs_owner,  _probs);
        free_array(_confs_owner,  _confs);
    };

    FixedEnvelope operator+(const FixedEnvelope& other) const;
    FixedEnvelope operator*(const FixedEnvelope& other) const;

    inline size_t    confs_no()  const { return _confs_no; }
    inline int       getAllDim() const { return allDim; }

    inline const double*   masses() const { return _masses; }
    inline const double*   probs()  const { return _probs; }
    inline const int*      confs()  const { return _confs; }

    //! Hand an array over as a pointer the caller must free() itself.
    /*! Copies the array first if it is not already a malloc()'d one -- see
        ReleasedArray and the release_*_with_deleter() overloads below for the
        variant that never copies. */
    inline double*   release_masses()     { return release_array(_masses_owner, _masses); }
    inline double*   release_probs()      { return release_array(_probs_owner,  _probs);  }
    inline int*      release_confs()      { return release_array(_confs_owner,  _confs);  }

    //! Hand an array over together with its deleter. Never copies.
    inline ReleasedArray release_masses_with_deleter() { return release_array_with_deleter(_masses_owner, _masses); }
    inline ReleasedArray release_probs_with_deleter()  { return release_array_with_deleter(_probs_owner,  _probs);  }
    inline ReleasedArray release_confs_with_deleter()  { return release_array_with_deleter(_confs_owner,  _confs);  }

    //! Drop all three arrays without freeing any of them: somebody else owns them.
    inline void      release_everything() { disown_array(_masses_owner, _masses); disown_array(_probs_owner, _probs); disown_array(_confs_owner, _confs); }


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

    //! Resize the three arrays to hold new_size configurations, keeping the
    //! _confs_no already stored -- the growth path, for fills that do not know
    //! how much they will produce.
    template<bool tgetConfs> void reallocate_memory(size_t new_size);

    //! Allocate the three arrays for exactly new_size configurations, keeping
    //! nothing -- for fills that counted their output up front.
    template<bool tgetConfs> void allocate_memory(size_t new_size);

    //! Store every remaining configuration of the generator's current layer.
    template<bool tgetConfs> void store_layer(IsoLayeredGenerator& generator);

    //! Store configurations of the generator's current layer, adding each one's
    //! probability to prob_so_far, until that reaches target_total_prob or the
    //! layer runs out; returns whether it was reached.
    template<bool tgetConfs> bool store_layer_to_prob(IsoLayeredGenerator& generator,
                                                      double& prob_so_far,
                                                      double target_total_prob);

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
        return FromThreshold(Iso(iso, true), _threshold, _absolute, tgetConfs);
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
        return FromTotalProb(Iso(iso, true), _target_total_prob, _optimize, tgetConfs);
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
        return FromStochastic(Iso(iso, true), _no_molecules, _precision, _beta_bias, tgetConfs);
    }

    static FixedEnvelope Binned(Iso&& iso, double target_total_prob, double bin_width, double bin_middle = 0.0);
    static FixedEnvelope Binned(const Iso& iso, double target_total_prob, double bin_width, double bin_middle = 0.0)
    {
        return Binned(Iso(iso, true), target_total_prob, bin_width, bin_middle);
    }
};

}  // namespace IsoSpec
