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

namespace IsoSpec
{

//! A class for allocating memory for the calculated configurations.
/*!
    We store isotopes in cells.
    Each cell comprizes enough bytes for the log-probability (needed as a key in the Priority Queue) and
    the representation of the counts of individual isotopes.
    Cells are hold in arrays.
    If there is not enough space in an array, a new one is generated.
    Pointers to tables are collected in prevTabs vector.
*/
class DirtyAllocator{
private:
    void*               currentTab;
    void*               currentConf;
    void*               endOfTablePtr;
    const int           tabSize;
    int                 cellSize;
    std::vector<void*>  prevTabs;
public:

    //! Constructor.
    /*!
        \param dim The total number of isotopes (the number of integers to be stored in one chunk of memory, 
        beside the space for log-probability and mass).
        \tabSize How many 'configurations' are to be stored in one array.
    */
    DirtyAllocator(const int dim, const int tabSize = 10000);
    
    //! Destructor.
    ~DirtyAllocator();

    //! Switch to a new table.
    /*!
        Pushes the previous table on prevTabs and allocates new table.
    */
    void shiftTables();

    //! Return the pointer to the new configuration in the current array. If no place, create a new array.
    inline void* newConf()
    {
        if (currentConf >= endOfTablePtr)
        {
            shiftTables();
        }

        void* ret = currentConf;
        currentConf = reinterpret_cast<char*>(currentConf) + cellSize;

        return ret;
    }

    //! Write the conf into the table and return pointer to that chunk.
    /*!
        \param conf A configuration (chunk of memory with log-probability and isotope counts).
    */
    inline void* makeCopy(const void* conf)
    {
        void* currentPlace = newConf();

        memcpy(currentPlace, conf, cellSize);

        return currentPlace;
    }

    //! Like @ref makeCopy, but you have to remember to deallocate the configuration yourself, externally.
    inline void* makeExternalCopy(const void* conf)
    {
        void* res = malloc(cellSize);

        memcpy(res, conf, cellSize);

        return res;
    }
};

} // namespace IsoSpec

