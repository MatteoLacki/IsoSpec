/*
 *   Copyright (C) 2015-2018 Mateusz Łącki and Michał Startek.
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

#include <string.h>
#include "conf.h"
#include "isoMath.h"
#include "misc.h"

namespace IsoSpec
{

//! The hash function class.
/*!
    Needed for the unordered-map.
*/
class KeyHasher
{
private:
    int dim;
public:
    //! Constructor.
    /*!
        \param dim the number of the ints that make up a configuration.
    */
    KeyHasher(int dim);

    //! The __call__ operator.
    /*!
        \param conf An array of integer counts.
        \return The hash for counts.
    */
    inline std::size_t operator()(const int* conf) const
    {
        // Following Boost... like, what the fuck are they think they doing????
        std::size_t seed = 0;
        for(int i = 0; i < dim; ++i )
            seed ^= conf[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    };
};


//! The equality of configurations operator.
/*!
    Needed for the unordered-map.
*/
class ConfEqual
{
private:
    int size;
public:
    //! Constructor.
    /*!
        \param dim the number of the ints that make up a configuration.
    */
    ConfEqual(int dim);

    //! The __call__ operator.
    /*!
        Let us quote the sacred MAN of memcmp:
        "The memcmp() function returns zero if the two strings are identical,
        otherwise returns the difference between the first two differing bytes
        (treated as unsigned char values, so that `\200' is greater than `\0',
        for example). Zero-length strings are always identical.  This behavior
        is not required by C and portable code should only depend on the sign of
        the returned value."

        \param conf1 An array of integer counts.
        \param conf2 An array of integer counts.
        \return Are conf1 and conf2 the same configuration?
    */
    inline bool operator()(const int* conf1, const int* conf2) const
    {
        // The memcmp() function returns zero if the two strings are identical, oth-
        // erwise returns the difference between the first two differing bytes
        // (treated as unsigned char values, so that `\200' is greater than `\0',
        // for example).  Zero-length strings are always identical.  This behavior
        // is not required by C and portable code should only depend on the sign of
        // the returned value.
        //                                          sacred man of memcmp.
        return memcmp(conf1, conf2, size) == 0;
    }
};


//! The class used for comparing the position of configurations in the order of descending probabilities.
/*!
    Needed for the priority queue.
*/
class ConfOrder
{
public:
    inline bool operator()(void* conf1,void* conf2) const
    {
        return *reinterpret_cast<double*>(conf1) < *reinterpret_cast<double*>(conf2);
    };
};


//! The class used for comparing the position of subisotopologues in the order of descending probabilities.
/*!
    Needed for the priority queue.
*/
class ConfOrderMarginal
{
    const double*  logProbs;
    int dim;
public:
    //! Constructor.
    /*!
        \param logProbs
        \param dim The number of isotopes.
    */
    ConfOrderMarginal(const double* logProbs, int dim);

    //! Constructor.
    /*!
        \param conf1 An array of integer counts.
        \param conf2 An array of integer counts.
        \return True if conf1 is less probable than conf2.
    */
    inline bool operator()(const Conf conf1, const Conf conf2)
    {
        return unnormalized_logProb(conf1,logProbs,dim) < unnormalized_logProb(conf2,logProbs,dim);
    };
};


//! The class used for comparing the position of subisotopologues in the order of descending probabilities.
/*!
    Needed for the priority queue.
*/
class ConfOrderMarginalDescending
{
//configurations comparator
    const double*  logProbs;
    int dim;
public:
    //! Contstructor.
    /*!
        \param logProbs
        \param dim The number of isotopes.
    */
    ConfOrderMarginalDescending(const double* logProbs, int dim);

    //! Constructor.
    /*!
        \param conf1 An array of integer counts.
        \param conf2 An array of integer counts.
        \return True if conf1 is more probable than conf2.
    */
    inline bool operator()(const Conf conf1, const Conf conf2)
    {
        return unnormalized_logProb(conf1,logProbs,dim) > unnormalized_logProb(conf2,logProbs,dim);
    };
};



template<typename T> class ReverseOrder
{
public:
    inline ReverseOrder() {};
    inline bool operator()(const T a,const T b) const { return a > b; };
};

template<typename T> class TableOrder
{
	const T* tbl;
public:
	inline TableOrder(const T* _tbl) : tbl(_tbl) {};
	inline bool operator()(unsigned int i, unsigned int j) { return tbl[i] < tbl[j]; };
};

} // namespace IsoSpec

#include "marginalTrek++.h"

namespace IsoSpec
{

class OrderMarginalsBySizeDecresing
{
    PrecalculatedMarginal const* const* const T;
public:
    inline OrderMarginalsBySizeDecresing(PrecalculatedMarginal const* const * _T) : T(_T) {};
    inline bool operator()(int m1, int m2)
    {
        return T[m1]->get_no_confs() > T[m2]->get_no_confs();
    }
};


} // namespace IsoSpec

