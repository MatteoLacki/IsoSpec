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

#include <cmath>
#include <atomic>

namespace IsoSpec
{

//! A class implementing the Shewchuk's summator algorithm for a more numerically stable summation of floating point numbers.
class SSummator
{
    std::vector<double> partials;
    int maxpart;
public:
    //! Constructor (sum defaults to zero).
    inline SSummator()
    { maxpart = 0; }

    //! Copy constructor.
    inline SSummator(SSummator& other)
    {
        this->partials = other.partials;
        this->maxpart = other.maxpart;
    }

    //! Add a number to the existing sum.
    /*!
        \param x A double floating point number to add.
    */
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

    //! Get the current value of the sum of the added floating point numbers.
    inline double get()
    {
        double ret = 0.0;
        for(int i=0; i<maxpart; i++)
            ret += partials[i];
        return ret;
    }
};



//! A class implementing the Kahan's summator algorithm for a more numerically stable summation of floating point numbers.
class Summator{
    // Kahan algorithm
   double sum;
   double c;

public:
    //! Constructor (sum defaults to zero).
    inline Summator()
    { sum = 0.0; c = 0.0;}

    //! Add a number to the existing sum.
    /*!
        \param x A double floating point number to add.
    */
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


//! A class implementing the trivial summator of floating point numbers.
class TSummator
{
    // Trivial algorithm, for testing only
    double sum;
public:
    //! Constructor (sum defaults to zero).
    inline TSummator()
    { sum = 0.0; }

    //! Add a number to the existing sum.
    /*!
        \param x A double floating point number to add.
    */
    inline void add(double what)
    {
        sum += what;
    }

    //! Get the current value of the sum of the added floating point numbers.
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
    //! Constructor (sum defaults to zero).
    inline ThreadSummator() : sum(0.0) {};

    //! Add a number to the existing sum.
    /*!
        \param x A double floating point number to add.
    */
    inline void add(double what)
    {
        double previous = sum.load(std::memory_order_relaxed);
        while(!sum.compare_exchange_weak(previous, previous+what, std::memory_order_relaxed)) {};
    }

    //! Get the current value of the sum of the added floating point numbers.
    inline double get()
    {
        return sum.load(std::memory_order_relaxed);
    }
};

} // namespace IsoSpec

