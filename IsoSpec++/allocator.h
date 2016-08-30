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


#ifndef ALLOCATOR_HPP
#define ALLOCATOR_HPP

#include <vector>
#include <iostream>
#include "conf.h"


template <typename T> inline void copyConf(
    const T* source, T* destination,
    int dim
){
    for(int i = 0; i < dim; i++) destination[i] = source[i];
}

template <typename T> class Allocator{
private:
    T*      currentTab;
    int currentId;
    const int       dim, tabSize;
    std::vector<T*>  prevTabs;
public:
    Allocator(const int dim, const int tabSize = 10000);
    ~Allocator();

    void shiftTables();

    inline T* newConf()
    {
        unsigned int idx = (currentId++) * dim;

        if (currentId >= tabSize)
        {
            shiftTables();
        }

        return &currentTab[ idx ];
    }

    inline T* makeCopy(const T* conf)
    {
        T* currentPlace = newConf();
        copyConf<T>( conf, currentPlace, dim );

        return currentPlace;
    }

    inline T* makeExternalCopy(const T* conf)
    {
        T* res = new T[dim];
        copyConf( conf, res, dim );

        return res;
    }
};

#endif
