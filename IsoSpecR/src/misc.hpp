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

#ifndef MISC_HPP
#define MISC_HPP

#include <iostream>
#include <tuple>
#include <vector>
#include "logFactorial.hpp"

inline double combinedSum(
    const int* conf, const std::vector<double>** valuesContainer, int dimNumber
){
    double res = 0.0;
    for(int i=0; i<dimNumber;i++)
        res += (*(valuesContainer[i]))[conf[i]];
    return res;
}

inline int* getConf(void* conf)
{
    return reinterpret_cast<int*>(
        reinterpret_cast<char*>(conf) + sizeof(double)
    );
}

inline double getLProb(void* conf)
{
    double ret = *reinterpret_cast<double*>(conf);
    return ret;
}


inline double logProb(const int* conf, const double* logProbs, int dim)
{
    int     N = 0;
    double  res = 0.0;

    for(int i=0; i < dim; i++)
    {
        N   += conf[i];
        res -= logFactorial(conf[i]);
        res += conf[i] * logProbs[i];
    }
    return res + logFactorial(N);
};

inline double mass(const int* conf, const double* masses, int dim)
{
    double res = 0.0;

    for(int i=0; i < dim; i++)
    {
        res += conf[i] * masses[i];
    }

    return res;
};


inline bool tupleCmp(
    std::tuple<double,double,int*> t1,
    std::tuple<double,double,int*> t2
){
    return std::get<1>(t1) > std::get<1>(t2);
};

template<typename T> void printArray(const T* array, int size)
{
    for (int i=0; i<size; i++)
        std::cout << array[i] << " ";
    std::cout << std::endl;
};

template<typename T> void printVector(const std::vector<T>& vec)
{
    printArray<T>(vec.data(), vec.size());
};


template<typename T> void printNestedArray(const T** array, const int* shape, int size)
{
    for (int i=0; i<size; i++)
        printArray(array[i], shape[i]);
    std::cout << std::endl;
};

#define mswap(x, y) swapspace = x; x = y; y=swapspace;

void* quickselect(void** array, int n, int start, int end);
#endif
