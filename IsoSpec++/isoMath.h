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
#include <fenv.h>

#define ISOSPEC_G_FACT_TABLE_SIZE 1024*1024*10

namespace IsoSpec
{

extern double* g_lfact_table;

//! Calculate \f$-\log(n!)\f$.
/*!
    \param n An integer.
    \return \f$-\log(n!)\f$.
*/
static inline double minuslogFactorial(int n) 
{ 
    if (n < 2) 
        return 0.0;
    if (g_lfact_table[n] == 0.0)
        g_lfact_table[n] = -lgamma(n+1);

    return g_lfact_table[n];
}

//! Calculate a quantile of the standard normal distibution.
/*!
    \param p A quantile to calculate.
*/
double NormalCDFInverse(double p);

//! Calculate a quantile of the standard normal distibution.
/*!
    \param p A quantile to calculate.
    \param mean The mean of the normal distribution.
    \param stdev The standard deviation of the normal distribution.
*/
double NormalCDFInverse(double p, double mean, double stdev);

//! Calculate the distribuant of the standard normal distibution.
/*!
    \param x The point at which to calculate the distribuant.
    \param mean The mean of the normal distribution.
    \param stdev The standard deviation of the normal distribution.
*/
double NormalCDF(double x, double mean, double stdev);

//! Calculate the density of the standard normal distibution.
/*!
    \param x The point at which to calculate the density.
    \param mean The mean of the normal distribution.
    \param stdev The standard deviation of the normal distribution.
*/
double NormalPDF(double x, double mean = 0.0, double stdev = 1.0);

} // namespace IsoSpec

