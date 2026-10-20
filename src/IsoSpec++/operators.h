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

#include <string.h>
#include "platform.h"
#include "conf.h"
#include "isoMath.h"
#include "misc.h"

namespace IsoSpec
{

class ConfOrder
{
// configurations comparator
 public:
    inline bool operator()(void* conf1, void* conf2) const
    {
        return *reinterpret_cast<double*>(conf1) < *reinterpret_cast<double*>(conf2);
    };
};

template<typename T> class TableOrder
{
        const T* tbl;
 public:
        inline explicit TableOrder(const T* _tbl) : tbl(_tbl) {}
        inline bool operator()(size_t i, size_t j) { return tbl[i] < tbl[j]; }
};

}  // namespace IsoSpec

#include "marginalTrek++.h"

class PrecalculatedMarginal;  // In case marginalTrek++.h us including us, and can't be included again...

namespace IsoSpec
{

template<typename T> class OrderMarginalsBySizeDecresing
{
    T const* const* const MT;
 public:
    explicit OrderMarginalsBySizeDecresing(T const* const * const MT_) : MT(MT_) {};
    inline bool operator()(size_t m1, size_t m2) { return MT[m1]->get_no_confs() > MT[m2]->get_no_confs(); }
};


}  // namespace IsoSpec
