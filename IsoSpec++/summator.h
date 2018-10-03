/*
 *   Copyright (C) 2015-2016 Mateusz Łącki and Michał Startek.
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

#ifndef SUMMATOR_HPP
#define SUMMATOR_HPP
#include <cmath>
#include <atomic>

class SSummator
{
    // Shewchuk algorithm
    std::vector<double> partials;
    int maxpart;
public:
    inline SSummator()
    { maxpart = 0; }

    inline SSummator(SSummator& other)
    {
        this->partials = other.partials;
        this->maxpart = other.maxpart;
    }
    inline void add(double x)
    {
        unsigned int i=0;
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
   double sum;
   double c;

public:
    inline Summator()
    { sum = 0.0; c = 0.0;}

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
    double sum;
public:
    inline TSummator()
    { sum = 0.0; }

    inline void add(double what)
    {
        sum += what;
    }
    inline double get()
    {
    	return sum;
    }
};

class ThreadSummator
{
    // Trivial but thread-safe summator
    std::atomic<double> sum;
public:
    inline void add(double what)
    {
        double previous = sum.load(std::memory_order_relaxed);
        while(!sum.compare_exchange_weak(previous, previous+what, std::memory_order_relaxed)) {};
    }
    inline double get()
    {
        return sum.load(std::memory_order_relaxed);
    }
};

#endif

