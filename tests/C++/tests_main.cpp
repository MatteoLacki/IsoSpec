// doctest entry point for the IsoSpec C++ test suite.
//
// This is the ONLY translation unit that defines doctest's main().  Every other
// test_*.cpp includes doctest.h without the IMPLEMENT macro and just registers
// TEST_CASEs.  All of them, plus unity-build.cpp, are linked into one binary
// per build configuration (see the Makefile).

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
