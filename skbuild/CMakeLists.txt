cmake_minimum_required(VERSION 3.15)

project(IsoSpecPy CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(PythonExtensions REQUIRED)
#find_package(Python COMPONENTS Interpreter Development REQUIRED)

add_library(IsoSpecCppPy SHARED ../IsoSpec++/unity-build.cpp)
target_compile_features(IsoSpecCppPy PRIVATE cxx_std_17)
if(MSVC)
        set(IsoSpecCppPy_WARNINGS /W3)
else()
        set(IsoSpecCppPy_WARNINGS -Wall -Wextra -Wconversion -Wno-sign-conversion)# -Wsign-conversion -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Woverloaded-virtual -Wpedantic)
endif()

target_compile_options(IsoSpecCppPy PRIVATE ${IsoSpecCppPy_WARNINGS})
#python_extension_module(IsoSpecCppPy)
install(TARGETS IsoSpecCppPy LIBRARY DESTINATION IsoSpecPy)
install(DIRECTORY ../IsoSpec++/ DESTINATION IsoSpecPy/IsoSpec++
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.hxx"
)
