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

#include <vector>
#include <iostream>
#include <string.h>
#include "conf.h"

namespace IsoSpec
{


//! Copy a configuration.
/*!
    \param source       The memory from which to copy.
    \param destination  The memoru unto which to copy.
    \param dim          The number of copies of object T in a configuration.
*/
template <typename T> inline void copyConf(
    const T* source, T* destination,
    int dim
){
    memcpy(destination, source, dim*sizeof(T));
}


//! A template class for allocating memory for the calculated configurations.
/*!
    We store isotopes in cells.
    Each cell comprizes enough bytes for the log-probability (needed as a key in the Priority Queue) and
    the representation of the counts of individual isotopes.
    Cells are hold in arrays.
    If there is not enough space in an array, a new one is generated.
    Pointers to tables are collected in prevTabs vector.
*/
template <typename T> class Allocator{
private:
    T*              currentTab;
    int             currentId;
    const int       dim, tabSize;
    std::vector<T*> prevTabs;
public:

    //! Constructor.
    /*!
        \param dim The number of objects T in a cell.
        \param How many 'configurations' are to be stored in one array.
    */
    Allocator(const int dim, const int tabSize = 10000);
    ~Allocator();


    //! Switch to a new table.
    /*!
        Pushes the previous table on prevTabs and allocates new table.
    */
    void shiftTables();

    //! Return the pointer to the new configuration in the current array. If no place, create a new array.
    inline T* newConf()
    {
        currentId++;

        if (currentId >= tabSize)
            shiftTables();

        return &(currentTab[ currentId * dim ]);
    }

     //! Write the conf into the table and return pointer to that chunk.
    /*!
        \param conf A configuration (chunk of memory with log-probability and isotope counts).
    */
    inline T* makeCopy(const T* conf)
    {
        T* currentPlace = newConf();
        copyConf<T>( conf, currentPlace, dim );

        return currentPlace;
    }

    //! Like @ref makeCopy, but you have to remember to deallocate the configuration yourself, externally.
    inline T* makeExternalCopy(const T* conf)
    {
        T* res = new T[dim];
        copyConf( conf, res, dim );

        return res;
    }
};

}

