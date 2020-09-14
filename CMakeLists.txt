#############################################################
#############################################################
# CMake configuration

cmake_minimum_required(VERSION 3.12)

# CMake script for the IsoSpec library
# Author: Anonymous and Filippo Rusconi
# This script was present in the source tree but was
# highly rudimentary. I have implemented most of what is
# here (Filippo Rusconi).
# Last Modifs: July 2020

############################################################
############################################################
# Basic information about project

project(IsoSpec
    DESCRIPTION "A program to calculate isotopic clusters"
    HOMEPAGE_URL "https://github.com/MatteoLacki/IsoSpec")

set(ISOSPEC_VERSION_MAJOR "2")
set(ISOSPEC_VERSION_MINOR "1")
set(ISOSPEC_VERSION_PATCH "3")
set(ISOSPEC_LIB_NAME "libIsoSpec++")
set(ISOSPEC_VERSION "${ISOSPEC_VERSION_MAJOR}.${ISOSPEC_VERSION_MINOR}.${ISOSPEC_VERSION_PATCH}")
set(VERSION 2.1.3)
set(ISOSPEC_LIB_VERSION ${ISOSPEC_VERSION})
set(ISOSPEC_LIB_SOVERSION ${ISOSPEC_VERSION_MAJOR})

# Command to enable debug mode
# cmake -DCMAKE_BUILD_TYPE=Debug 


# This setting is also typically overridden by 
# configuring the build system like this:
# cmake -DCMAKE_INSTALL_PREFIX=/usr
# It needs to be set before include(GNUInstallDirs)
# below, because it roots all the other dirs' configs.

if(NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX /usr)
endif()

# This module sets a number of installation directories 
# according to the GNU coding standards. It takes the
# CMAKE_INSTALL_PREFIX and then add the directories to it.
include(GNUInstallDirs)

# Add folder where are supportive functions
# that contain toolchains: for UNIX/WIN/MAC
# that contain modules useful on any system.
set(CMAKE_UTILS_PATH ${CMAKE_SOURCE_DIR}/CMakeStuff)
set(CMAKE_TOOLCHAINS_PATH ${CMAKE_UTILS_PATH}/toolchains)
set(CMAKE_MODULE_PATH ${CMAKE_UTILS_PATH}/modules)
#message("CMAKE_MODULE_PATH:" ${CMAKE_MODULE_PATH})

# This says that when compiling, .h files must also be
# searche in the current directory where cpp files are 
# compiled.
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Include the system's uname that fills in SYSTEM_UNAME_S.
# Sets WIN64 if SYSTEM_UNAME_S is "^.*MING64.*"
include(${CMAKE_UTILS_PATH}/systemUname.cmake)

# Include the various colors we want to use in the output
include(${CMAKE_UTILS_PATH}/outputColors.cmake)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

message("")
message("\n${BoldRed}Configuring build for project ${PROJECT_NAME}${ColourReset}\n")
message("")

# This export will allow using the flags to be used by
# youcompleteme (vim plugin).
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
    execute_process( COMMAND cmake -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
        )
endif()


# We want C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "${BoldGreen}CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}${ColourReset}")

#############################################################
# We do not want warnings for unknown pragmas:
message(STATUS "${BoldGreen}Setting definition -Wno-unknown-pragmas.${ColourReset}")
add_definitions(-Wno-unknown-pragmas)

# Enable warnings and possibly treat them as errors
message(STATUS "${BoldGreen}Setting definition -Wall.${ColourReset}")
add_definitions(-Wall)
message(STATUS "${BoldGreen}Setting definition -Wextra.${ColourReset}")
add_definitions(-Wextra)
message(STATUS "${BoldGreen}Setting definition -pedantic.${ColourReset}")
add_definitions(-pedantic)

if(WARN_AS_ERROR)
    message(STATUS "${BoldYellow}Setting definition -Werror.${ColourReset}")
    add_definitions(-Werror)
endif()

message(STATUS "${BoldGreen}CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}${ColourReset}")
message(STATUS "${BoldGreen}CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}${ColourReset}")


#############################################################
#############################################################
# Platform-specific CMake configuration

# MXE is a pretty cool cross-compilation environment that allows one
# to build Win software right into GNU/Linux.
if(MXE)

    include(${CMAKE_TOOLCHAINS_PATH}/mxe-toolchain.cmake)
    #include(${CMAKE_TOOLCHAINS_PATH}/mxe-toolchain-olivier.cmake)

    # Set the name to the systemUname variable because in this situation that name
    # is not found, it it passed as a flag in the command line.
    set(SYSTEM_UNAME_S "mxe")

elseif(UNIX AND NOT APPLE)

    include(${CMAKE_TOOLCHAINS_PATH}/unix-toolchain.cmake)

elseif(WIN32 OR _WIN32)

    if(WIN10MINGW64)
        include(${CMAKE_TOOLCHAINS_PATH}/win10-mingw64-toolchain.cmake)
    endif()

elseif(APPLE)

    include(${CMAKE_TOOLCHAINS_PATH}/apple-macport-toolchain.cmake)

endif()

message("")
message("${BoldRed}The build toolchain is: ${SYSTEM_UNAME_S}${ColourReset}")
message("")


#############################################################
#############################################################
# Essential software configuration

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Type of build, options are: None, Debug, Release, RelWithDebInfo, MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES "Release")
    #message(STATUS "Compiling in release mode.")
    add_definitions("-DQT_NO_DEBUG_OUTPUT")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    #message(STATUS "Compiling in debug mode.")
    message(STATUS "Add definition -ggdb3 to format debug output for GDB.")
    add_definitions(-ggdb3)
endif()

if(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    #message(STATUS "Compiling in release with debug info mode.")
endif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")

message(STATUS "${BoldYellow}CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}.${ColourReset}")
message(STATUS "${BoldYellow}CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}.${ColourReset}")


# message(STATUS CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH})

# Now actually deal with the C++ source code.
add_subdirectory(IsoSpec++)
