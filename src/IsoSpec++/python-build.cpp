/*
 *   Copyright (C) 2025 Mateusz Łącki and Michał Startek.
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

#define ISOSPEC_BUILDING_PYTHON

#include "platform.h"


#if ISOSPEC_TEST_WE_ARE_ON_WINDOWS

#define ISOSPEC_C_API __declspec(dllexport)

#include <Python.h>

// Provide a dummy PyInit function on Windows/MSVC.
// We're not using it, as we'll load using CFFI - but it's easier
// than fighting with the build system.
extern "C" {
    __declspec(dllexport) PyObject* PyInit_IsoSpecCppPy(void) { return nullptr; }
}

#endif

#include "unity-build.cpp"  // NOLINT(build/include)
