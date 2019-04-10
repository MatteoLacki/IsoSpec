project("IsoSpec")
cmake_minimum_required(VERSION 2.6)

set(my_sources
cwrapper.cpp
allocator.cpp
dirtyAllocator.cpp
isoSpec++.cpp
isoMath.cpp
marginalTrek++.cpp
operators.cpp
element_tables.cpp
misc.cpp
tabulator.cpp
)

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
  CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
  if (WITH_FPIC)
    add_definitions(-fPIC)
  endif()

  # Only add those definitions on non-windows systems
  add_definitions(-std=c++17 -Wall -pedantic -Wextra)
else()
  ## On MSVS we need this for mmap
  set(my_sources ${my_sources} mman.cpp)
  add_definitions(-DMMAN_LIBRARY)
endif()

add_library(IsoSpec SHARED ${my_sources})

configure_file(IsoSpecConfig.cmake.in "${PROJECT_BINARY_DIR}/IsoSpecConfig.cmake" @ONLY)

export(TARGETS IsoSpec FILE IsoSpecLibrary.cmake)

