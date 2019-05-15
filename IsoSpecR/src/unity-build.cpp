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

#include "platform.h"

#if !ISOSPEC_BUILDING_R

#if !ISOSPEC_GOT_SYSTEM_MMAN && ISOSPEC_GOT_MMAN
    #include "mman.cpp"
#endif

// A poor-man's replacement for LTO. We're small enough that we can do that.

#include "allocator.cpp"
#include "dirtyAllocator.cpp"
#include "isoSpec++.cpp"
#include "isoMath.cpp"
#include "marginalTrek++.cpp"
#include "operators.cpp"
#include "element_tables.cpp"
#include "cwrapper.cpp"
#include "fixedEnvelopes.cpp"
#include "misc.cpp"

#endif
