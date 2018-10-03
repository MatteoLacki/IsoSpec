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
)
## #include "allocator.cpp"
## #include "dirtyAllocator.cpp"
## #include "isoSpec++.cpp"
## #include "isoMath.cpp"
## #include "marginalTrek++.cpp"
## #include "operators.cpp"
## #include "element_tables.cpp"
## #include "misc.cpp"
## #include "spectrum2.cpp"
## #include "cwrapper.cpp"
## #include "tabulator.cpp"

# CXXFLAGS=

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
  CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
  if (WITH_FPIC)
    add_definitions(-fPIC)
  endif()

  # Only add those definitions on non-windows systems
  add_definitions(-std=c++11 -Wall -pedantic -Wextra)

endif()

add_library(IsoSpec ${my_sources})

add_executable(mytest test.c)
target_link_libraries(mytest IsoSpec)

configure_file(IsoSpecConfig.cmake.in "${PROJECT_BINARY_DIR}/IsoSpecConfig.cmake" @ONLY)

#### # Create imported target Crawdad
#### # adapted from http://www.cmake.org/Wiki/CMake/Tutorials/How_to_create_a_ProjectConfig.cmake_file#The_FooBar.2FFooBarBuildTreeSettings.cmake.in_file
#### 
#### # Compute paths
#### get_filename_component(CRAWDAD_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
#### if(EXISTS "${CRAWDAD_CMAKE_DIR}/CMakeCache.txt")
####    # In build tree
####    set(CRAWDAD_INCLUDE_DIRS
####      "/home/hr/projects/crawdad/"
####      "/home/hr/projects/crawdad/msmat")
#### else()
####   set(CRAWDAD_INCLUDE_DIRS "${CRAWDAD_CMAKE_DIR}/")
#### endif()
#### 
#### # Our library dependencies (contains definitions for IMPORTED targets)
#### include("${CRAWDAD_CMAKE_DIR}/CrawdadLibrary.cmake")
#### 
#### # These are IMPORTED targets created by CrawdadLibrary.cmake
#### set(CRAWDAD_LIBRARIES Crawdad)


export(TARGETS IsoSpec FILE CrawdadLibrary.cmake)

