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



#include "misc.hpp"

#define mswap(x, y) swapspace = x; x = y; y=swapspace;



void* quickselect(void** array, int n, int start, int end)
{
    void* swapspace;

    if(start == end)
        return array[start];

    while(true)
    {
        // Partition part
        int len = end - start;
        int pivot = rand() % len + start;
        void* pval = array[pivot];
        double pprob = getLProb(pval);
        mswap(array[pivot], array[end-1]);
        int loweridx = start;
        for(int i=start; i<end-1; i++)
        {
            if(getLProb(array[i]) < pprob)
            {
                mswap(array[i], array[loweridx]);
                loweridx++;
            }
        }
        mswap(array[end-1], array[loweridx]);

        // Selection part
        if(n==loweridx)
            return array[n];
        if(n<loweridx)
            end = loweridx;
        else
            start = loweridx+1;
    };
}

