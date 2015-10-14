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

#ifndef ALLOCATOR_HPP
#define ALLOCATOR_HPP

#include <vector>
#include <iostream>
#include "conf.hpp"


template <typename T> inline void copyConf(
    const T* source, T* destination,
    int dim
){
    for(int i = 0; i < dim; i++) destination[i] = source[i];
};

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
