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

#ifndef SUMMATOR_HPP
#define SUMMATOR_HPP
#include <cmath>

class SSummator
{
    // Schewchuk algorithm
    std::vector<double> partials;
    int maxpart;
public:
    inline SSummator(int plen = 10000)
    {
        maxpart = 0;
    }
    inline SSummator(SSummator& other)
    {
        this->partials = other.partials;
        this->maxpart = other.maxpart;
    }
    inline void add(double x)
    {
        int i=0;
        for(int pidx=0; pidx<maxpart; pidx++)
        {
            double y = partials[pidx];
            if(std::abs(x) < std::abs(y))
                std::swap(x, y);
            double hi = x+y;
            double lo = y-(hi-x);
            if(lo != 0.0)
            {
                partials[i] = lo;
                i += 1;
            }
            x = hi;
        }
        while(partials.size() <= i)
            partials.push_back(0.0);
        partials[i] = x;
        maxpart = i+1;
    }
    inline double get()
    {
        double ret = 0.0;
        for(int i=0; i<maxpart; i++)
            ret += partials[i];
        return ret;
    }
};







class Summator{
    // Kahan algorithm
   double sum = 0.0;
   double c = 0.0;

public:
    inline void add(double what)
    {
        double y = what - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    inline double get()
    {
        return sum;
    }
};

class TSummator
{
    // Tirival algorithm
    double sum = 0.0;
public:
    inline void add(double what)
    {
        sum += what;
    }
    inline double get()
    {
    	return sum;
    }
};

#endif

