#include "platform.h"

#if !ISOSPEC_BUILDING_R

#if !ISOSPEC_GOT_SYSTEM_MMAN && ISOSPEC_GOT_MMAN
	#include "mman.cpp"
#endif

#include "allocator.cpp"
#include "dirtyAllocator.cpp"
#include "isoSpec++.cpp"
#include "isoMath.cpp"
#include "marginalTrek++.cpp"
#include "operators.cpp"
#include "element_tables.cpp"
#include "cwrapper.cpp"
#include "tabulator.cpp"
#include "misc.cpp"

#endif
