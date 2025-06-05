cmake_minimum_required(VERSION 3.15)

project(IsoSpecPy CXX)

#find_package(PythonExtensions REQUIRED)
#find_package(Python COMPONENTS Interpreter Development REQUIRED)

add_library(IsoSpecCppPy SHARED ../IsoSpec++/unity-build.cpp)
target_compile_features(IsoSpecCppPy PRIVATE cxx_std_17)
#python_extension_module(IsoSpecCppPy)
install(TARGETS IsoSpecCppPy LIBRARY DESTINATION IsoSpecPy)
install(DIRECTORY ../IsoSpec++/ DESTINATION IsoSpecPy/IsoSpec++
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.hxx"
)
