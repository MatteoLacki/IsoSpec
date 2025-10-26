

#define ISOSPEC_BUILDING_PYTHON

#include "platform.h"


#if ISOSPEC_TEST_WE_ARE_ON_WINDOWS

#define ISOSPEC_C_API __declspec(dllexport)

#include <Python.h>

// Provide a dumy PyInit function on Windows/MSVC.
// We're not using it, as we'll load using CFFI - but it's easier
// than fighting with the build system.
extern "C" {
    __declspec(dllexport) PyObject* PyInit_IsoSpecCppPy(void) { return nullptr; };
}

#endif

#include "unity-build.cpp"
