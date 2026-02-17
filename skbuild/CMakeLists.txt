cmake_minimum_required(VERSION 3.15)

project(IsoSpecPy CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find Python and nanobind
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Try to find nanobind, if not found use FetchContent
find_package(nanobind CONFIG)
if(NOT nanobind_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        nanobind
        GIT_REPOSITORY https://github.com/wjakob/nanobind.git
        GIT_TAG v2.4.0
    )
    FetchContent_MakeAvailable(nanobind)
endif()

# Build the C++ library (for compatibility, though not strictly needed with nanobind)
add_library(IsoSpecCppPy SHARED ../src/IsoSpec++/unity-build.cpp)
target_compile_features(IsoSpecCppPy PRIVATE cxx_std_17)
install(TARGETS IsoSpecCppPy LIBRARY DESTINATION IsoSpecPy)

# Build the nanobind extension module
nanobind_add_module(_isospec_nb
    STABLE_ABI
    NB_STATIC
    ../src/IsoSpecPy/isospec_nb.cpp
    ../src/IsoSpec++/unity-build.cpp
)

target_compile_features(_isospec_nb PRIVATE cxx_std_17)
install(TARGETS _isospec_nb LIBRARY DESTINATION IsoSpecPy)

# Install header files
install(DIRECTORY ../src/IsoSpec++/ DESTINATION IsoSpecPy/IsoSpec++
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.hxx"
)
