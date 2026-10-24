/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
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

#include <cstring>
#include "conf.h"
#include "pod_vector.h"

namespace IsoSpec
{

// Configurations per allocator slab.  Small on purpose: with glibc, a slab of
// up to 1032 bytes (4 isotopes x 64 ints = 1024) is recycled through the
// per-thread tcache when the next molecule is built, so many small molecules
// in a row run 16-20% faster than with 1000 (2000 tryptic peptides, threshold
// generator, x86-64-v3 and v4).  Single large calculations do not notice.
constexpr int default_tab_size = 64;

template <typename T> inline void copyConf(
    const T* source, T* destination,
    int dim
){
    memcpy(destination, source, dim*sizeof(T));
}

template <typename T> class Allocator
{
 private:
    T*      currentTab;
    int currentId;
    const int       dim, tabSize;
    pod_vector<T*>  prevTabs;

 public:
    explicit Allocator(const int dim, const int tabSize = default_tab_size);
    ~Allocator();

    Allocator(const Allocator& other) = delete;
    Allocator& operator=(const Allocator& other) = delete;

    void shiftTables();

    inline T* newConf()
    {
        currentId++;

        if (currentId >= tabSize)
            shiftTables();

        return &(currentTab[ currentId * dim ]);
    }

    inline T* makeCopy(const T* conf)
    {
        T* currentPlace = newConf();
        copyConf<T>( conf, currentPlace, dim );

        return currentPlace;
    }
};

}  // namespace IsoSpec
