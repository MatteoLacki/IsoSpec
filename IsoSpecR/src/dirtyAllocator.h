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

#ifndef DIRTY_ALLOCATOR_HPP
#define DIRTY_ALLOCATOR_HPP

#include <vector>
#include <iostream>
#include <string.h>

class DirtyAllocator{
private:
    void*   currentTab;
    void*   currentConf;
    void*   endOfTablePtr;
    const int       tabSize, cellSize;
    std::vector<void*>  prevTabs;
public:
    DirtyAllocator(const int dim, const int tabSize = 10000);
    ~DirtyAllocator();

    void shiftTables();

    inline void* newConf()
    {
        if (currentConf >= endOfTablePtr)
        {
            shiftTables();
        }

        void*  ret = currentConf;
        currentConf = reinterpret_cast<char*>(currentConf) + cellSize;

        return ret;
    }

    inline void* makeCopy(const void* conf)
    {
        void* currentPlace = newConf();

        memcpy(currentPlace, conf, cellSize);

        return currentPlace;
    }

    inline void* makeExternalCopy(const void* conf)
    {
        void* res = malloc(cellSize);

        memcpy(res, conf, cellSize);

        return res;
    }
};

#endif
