# File:///home/langella/developpement/git/pappsomspp/CMakeStuff/toolchains/mxe-toolchain.cmake# 
# This file should be included if the command line reads like this:
# x86_64-w64-mingw32.shared-cmake -DCMAKE_BUILD_TYPE=Release -DMXE=1 ..

MESSAGE("MXE (M cross environment) https://mxe.cc/")
message("Please run the configuration like this:")
message("x86_64-w64-mingw32.shared-cmake -DMXE=1 -G \"Unix Makefiles\" -DCMAKE_BUILD_TYPE=Release ../../development")


set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES ${HOME_DEVEL_DIR}/mxe/usr/x86_64-w64-mingw32.shared/include)
set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES ${HOME_DEVEL_DIR}/mxe/usr/x86_64-w64-mingw32.shared/include)


if(WIN32 OR _WIN32)
	message(STATUS "Building with WIN32 defined.")
endif()


find_package(ZLIB REQUIRED)


set(QUAZIP_FOUND 1)
set(QUAZIP_INCLUDE_DIR "${HOME_DEVEL_DIR}/quazip5/development")
set(QUAZIP_LIBRARIES "${HOME_DEVEL_DIR}/quazip5/build-area/mxe/libquazip5.dll")
set(QUAZIP_ZLIB_INCLUDE_DIR ${ZLIB_INCLUDE_DIRS})
set(QUAZIP_INCLUDE_DIRS ${QUAZIP_INCLUDE_DIR} ${QUAZIP_ZLIB_INCLUDE_DIR})


message(STATUS "QUAZIP_INCLUDE_DIR :${QUAZIP_INCLUDE_DIR}")
