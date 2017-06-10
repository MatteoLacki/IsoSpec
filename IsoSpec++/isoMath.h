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


#ifndef ISOMATH_HPP
#define ISOMATH_HPP

#include <cmath>

static inline double logFactorial(int n) { return lgamma(n+1); }
double NormalCDFInverse(double p);
double NormalCDFInverse(double p, double mean, double stdev);
double NormalCDF(double x, double mean, double stdev);
double NormalPDF(double x, double mean = 0.0, double stdev = 1.0);

inline unsigned int next_pow2(unsigned int base)
{
	// Rounds up a number to the next power of 2.
	// Base has to be >=1.
	// TODO: replace with one of those magic super-fast
	// bit operations.
	unsigned int ret = 1;
	while (ret <= base)
	    ret *= 2;
	return ret;
}

inline unsigned int prev_pow2(unsigned int arg)
{
	// computes the largest power of 2 which is lower than arg. 
	// Arg has to be > 1.
	unsigned int prev = 1;
	unsigned int ret = 1;
	while(ret <= arg)
	{
	    ret *= 2;
	    prev = ret;
	}
	return prev;
}

inline unsigned int floor_log2(unsigned int arg)
{
    // arg must be > 0
    unsigned int ret = 0;
    while(arg > 1)
    {
        arg /= 2;
        ret++;
    }
    return ret;
}

#endif /* ISOMATH_HPP */
