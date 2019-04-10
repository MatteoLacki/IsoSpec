/*
 *   Copyright (C) 2015-2019 Mateusz Łącki and Michał Startek.
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

#if !defined(ISOSPEC_G_FACT_TABLE_SIZE)
// 10M should be enough for anyone, right?
// Actually, yes. If anyone tries to input a molecule that has more than 10M atoms,
// he deserves to get an exception thrown in his face.
#define ISOSPEC_G_FACT_TABLE_SIZE 1024*1024*10
#endif

namespace IsoSpec
{

extern double* g_lfact_table;

static inline double minuslogFactorial(int n)
{
    if (n < 2)
        return 0.0;
    if (g_lfact_table[n] == 0.0)
        g_lfact_table[n] = -lgamma(n+1);

    return g_lfact_table[n];
}

const double pi = 3.14159265358979323846264338328;
const double log2pluslogpi = log(2.0) + log(pi);

double NormalCDFInverse(double p);
double NormalCDFInverse(double p, double mean, double stdev);
double NormalCDF(double x, double mean, double stdev);
double NormalPDF(double x, double mean = 0.0, double stdev = 1.0);

// Returns lower incomplete gamma function of a/2, x, where a is int and > 0.
double LowerIncompleteGamma2(int a, double x);

// Returns y such that LowerIncompleteGamma2(a, y) == x. Approximately.
double InverseLowerIncompleteGamma2(int a, double x);

// Computes the inverse Cumulative Distribution Funcion of the Chi-Square distribution with k degrees of freedom
inline double InverseChiSquareCDF2(int k, double x)
{
    return InverseLowerIncompleteGamma2(k, x*tgamma(static_cast<double>(k)/2.0)) * 2.0;
}



} // namespace IsoSpec

