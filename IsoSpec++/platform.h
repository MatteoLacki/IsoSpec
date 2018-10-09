
#pragma once

#if !defined(ISOSPEC_BUILDING_R)
#define ISOSPEC_BUILDING_R false
#endif

#if !defined(ISOSPEC_BUILDING_CPP)
#define ISOSPEC_BUILDING_CPP true
#endif

#if !defined(ISOSPEC_BUILDING_PYTHON)
#define ISOSPEC_BUILDING_PYTHON false
#endif


#if defined(__unix__) || defined(__unix) || \
        (defined(__APPLE__) && defined(__MACH__))
#include <sys/mman.h>
#define ISOSPEC_WE_ARE_ON_UNIX_YAY true
#define ISOSPEC_WE_ARE_ON_WINDOWS false
#define ISOSPEC_GOT_SYSTEM_MMAN true
#define ISOSPEC_GOT_MMAN true
#elif defined(__MINGW32__) || defined(_WIN32)
#include "mman.h"
#define ISOSPEC_WE_ARE_ON_UNIX_YAY false
#define ISOSPEC_WE_ARE_ON_WINDOWS true
#define ISOSPEC_GOT_SYSTEM_MMAN false
#define ISOSPEC_GOT_MMAN true
#else
#include <stdlib.h>     /* malloc, free, rand */
#define ISOSPEC_WE_ARE_ON_UNIX_YAY false /* Well, probably... */
#define ISOSPEC_WE_ARE_ON_WINDOWS false
#define ISOSPEC_GOT_SYSTEM_MMAN false
#define ISOSPEC_GOT_MMAN false
#endif


