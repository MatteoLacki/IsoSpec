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


#include "misc.h"
#include <utility>
#include "platform.h"
#include "isoMath.h"



namespace IsoSpec
{

void* quickselect(void ** array, size_t n, size_t start, size_t end)
{
    if(start == end)
        return array[start];

    while(true)
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
        void* pval = array[pivot];
        double pprob = getLProb(pval);
        std::swap(array[pivot], array[end-1]);
        size_t loweridx = start;
        for(size_t i = start; i < end-1; i++)
        {
            if(getLProb(array[i]) < pprob)
            {
                std::swap(array[i], array[loweridx]);
                loweridx++;
            }
        }
        std::swap(array[end-1], array[loweridx]);

        // Selection part
        if(n == loweridx)
            return array[n];
        if(n < loweridx)
            end = loweridx;
        else
            start = loweridx+1;
    };
}

}  // namespace IsoSpec
