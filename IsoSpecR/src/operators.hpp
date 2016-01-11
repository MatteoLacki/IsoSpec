/*
 *   Copyright (C) 2015 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License
 *   version 3, as published by the Free Software Foundation.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with IsoSpec.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPERATORS_HPP
#define OPERATORS_HPP

#include <string.h>
#include "conf.hpp"
#include "logFactorial.hpp"
#include "misc.hpp"



class KeyHasher
{
private:
    int dim;
public:
    KeyHasher(int dim);

    inline std::size_t operator()(const int* conf) const
    {
        // Following Boost...
        std::size_t seed = 0;
        for(int i = 0; i < dim; ++i )
            seed ^= conf[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    };
};


class ConfEqual
{
private:
    int size;
public:
    ConfEqual(int dim);

    inline bool operator()(const int* conf1, const int* conf2) const
    {
        return not memcmp(conf1, conf2, size);
    }
};


class ConfOrder
{
//configurations comparator
public:
    inline bool operator()(void* conf1,void* conf2) const
    {
        return *reinterpret_cast<double*>(conf1) < *reinterpret_cast<double*>(conf2);
    };
};



class ConfOrderMarginal
{
//configurations comparator
    const double*  logProbs;
    int dim;
public:
    ConfOrderMarginal(const double* logProbs, int dim);

    inline bool operator()(const Conf conf1, const Conf conf2)
    {// Return true if conf1 is less probable than conf2.
        return logProb(conf1,logProbs,dim) < logProb(conf2,logProbs,dim);
    };
};



#endif
