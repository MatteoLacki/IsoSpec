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

#include <iostream>
#include <stdlib.h>
#include "dirtyAllocator.hpp"


DirtyAllocator::DirtyAllocator(
    const int dim, const int tabSize
): tabSize(tabSize), cellSize( sizeof(double) + sizeof(int) * dim )
{
    currentTab      = malloc( cellSize * tabSize );
    currentConf     = currentTab;
    endOfTablePtr = reinterpret_cast<char*>(currentTab) + cellSize*tabSize;
};


DirtyAllocator::~DirtyAllocator()
{
    for(unsigned int i = 0; i < prevTabs.size(); ++i) free(prevTabs[i]);
    free(currentTab);
};

void DirtyAllocator::shiftTables()
{
    prevTabs.push_back(currentTab);

    currentTab              = malloc( cellSize * tabSize );
    currentConf             = currentTab;
    endOfTablePtr   = reinterpret_cast<char*>(currentTab) + cellSize*tabSize;
};
