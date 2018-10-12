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

#include <iostream>
#include <tuple>
#include <vector>
#include <fenv.h>
#include "isoMath.h"

namespace IsoSpec
{
//! Sum only those values contained in tables in valuesContainer that are selected by indices contained in conf.
/*!
    \param conf             An array of indices of elements in tables in valuesContainer.
    \param valuesContainer  A vector of tables with double floats, e.g. with masses or log-probabilities of subsequent subisotopologues.
    \param dimNumber        The total number of elements (i.e. the total number of subisotopologues that make up an isotopologue).
    \return The sum of selected elements in valuesContainer.
*/
inline double combinedSum(
    const int* conf, const std::vector<double>** valuesContainer, int dimNumber
){
    double res = 0.0;
    for(int i=0; i<dimNumber;i++)
        res += (*(valuesContainer[i]))[conf[i]];
    return res;
}

//! Get the counts of isotopes from the customized configuration-memory-chunk.
/*!
The configuration comprises as many bytes as to represent the log-probability and the counts of isotopologues.
This function selects those bytes that correspond to the counts of isotopologues and casts it into a usual array of ints.
*/
inline int* getConf(void* conf)
{
    return reinterpret_cast<int*>(
        reinterpret_cast<char*>(conf) + sizeof(double)
    );
}

//! Get the log-probability from the customized configuration-memory-chunk.
/*!
The configuration comprises as many bytes as to represent the log-probability and the counts of isotopologues.
This function selects the log-probability and casts it into a usual double precision float.
*/
inline double getLProb(void* conf)
{
    double ret = *reinterpret_cast<double*>(conf);
    return ret;
}

//! Calculate the isotope-count dependent part of the log-probability of an isotopologue.
/*!
    The probability of an isotopologue equals 
    \f$ \prod_{e\in\mathcal{E}} \frac{n_e!}{n_{e0},\dots,n_{e,i_e-1}} p_{e,0}^{n_{e,0}} \dots p_{e,i_e-1}^{n_{e,i_e-1}} \f$,
    where \f$n_{ej}\f$ is the count of element \f$e\f$'s \f$j^\text{th}\f$ isotope, and \f$p_{ej}\f$ are its abundances reported by IUPAC.
    In each multinomial distribution, \f$n_e!\f$ is constant under changes of isotope counts and so, we calculate it only once.

    \param conf The configuration (counts of isotopes).
    \param logProbs Table with the natural frequencies of individual isotopes.
    \param dim Total number of isotopes an isotopologue comprises.
    \return The log-probability of an isotopologue without the constant part.
*/
inline double unnormalized_logProb(const int* conf, const double* logProbs, int dim)
{
    double  res = 0.0;

    int curr_method = fegetround();

    fesetround(FE_TOWARDZERO);

    for(int i=0; i < dim; i++)
        res += minuslogFactorial(conf[i]);

    fesetround(FE_UPWARD);

    for(int i=0; i < dim; i++)
        res += conf[i] * logProbs[i];

    fesetround(curr_method);

    return res;
}

//! Calculate the mass of an isotopologue.
/*!
    The mass of an isotopologue equals \f$\sum_{e\in\mathcal{E}} \sum_{i=0}^{i_e-1} m_{e,i}n_{e,i}\f$,
    where \f$n_{e,j}\f$ is the count of element \f$e\f$'s \f$j^\text{th}\f$ isotope and \f$m_{e,j}\f$
    is its mass in daltons reported by IUPAC.

    \param conf The configuration (counts of isotopes).
    \param masses Table with the masses of individual isotopes.
    \param dim Total number of isotopes an isotopologue comprises.
    \return The mass of the isotopologue.
*/
inline double mass(const int* conf, const double* masses, int dim)
{
    double res = 0.0;

    for(int i=0; i < dim; i++)
    {
        res += conf[i] * masses[i];
    }

    return res;
}

//! Is the second element of the first tuple larger than that of the second tuple?
inline bool tupleCmp(
    std::tuple<double,double,int*> t1,
    std::tuple<double,double,int*> t2
){
    return std::get<1>(t1) > std::get<1>(t2);
}

//! Print an array of objects.
template<typename T> void printArray(const T* array, int size)
{
    for (int i=0; i<size; i++)
        std::cout << array[i] << " ";
    std::cout << std::endl;
}

//! Print a vector of tables of objects T.
template<typename T> void printVector(const std::vector<T>& vec)
{
    printArray<T>(vec.data(), vec.size());
}

//! Print a nester array.
template<typename T> void printNestedArray(const T** array, const int* shape, int size)
{
    for (int i=0; i<size; i++)
        printArray(array[i], shape[i]);
    std::cout << std::endl;
}

#define mswap(x, y) swapspace = x; x = y; y=swapspace;

//! Quickly select the n'th positional statistic, including the weights.
void* quickselect(void** array, int n, int start, int end);

//! Copy an array of T objects.
template <typename T> inline static T* array_copy(const T* A, int size)
{
    T* ret = new T[size];
    memcpy(ret, A, size*sizeof(T));
    return ret;
}

//! Deallocate table of objects T.
template<typename T> void dealloc_table(T* tbl, int dim)
{
    for(int i=0; i<dim; i++)
    {
        delete tbl[i];
    }
    delete[] tbl;
}

} // namespace IsoSpec

