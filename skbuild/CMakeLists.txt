cmake_minimum_required(VERSION 3.15)

project(IsoSpecPy CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(PythonExtensions REQUIRED)
#find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Runtime ISA dispatch. This matters most here: the wheel is built with no
# -march at all, so without it every user gets the x86-64 baseline -- 2-double
# vectors -- no matter what their CPU can do.
include(${CMAKE_CURRENT_LIST_DIR}/../src/IsoSpec++/IsaDispatch.cmake)

add_library(IsoSpecCppPy SHARED ../src/IsoSpec++/unity-build.cpp ${isospec_ISA_SRCS})
target_compile_features(IsoSpecCppPy PRIVATE cxx_std_20)

# macOS wheels are built with Homebrew GCC (Apple's libc++ has no
# <experimental/simd>), and Homebrew bottles its libstdc++.6.dylib for the
# *current* macOS -- minos 26.0 today. delocate refuses to bundle a dylib whose
# minimum target is newer than the wheel's (11.0 on arm64), and hard-fails the
# build. Linking the C++ runtime statically leaves libSystem as the only
# dependency, so there is nothing for delocate to reject. The
# "built for newer macOS version" warnings this draws from ld are the archive
# members' own load commands and are harmless.
#
# Plain -static-libstdc++ is not enough, though: the archive members keep
# default visibility, so the dylib would re-export ~2000 libstdc++ symbols.
# Mach-O coalesces weak definitions (templates, inline functions, vtables)
# across *all* loaded images regardless of two-level namespace, so a process
# that also loads anything built against Homebrew's libstdc++.6.dylib ends up
# with two copies of that runtime spliced together and aborts inside dlopen --
# silently, with no dyld message and nothing Python can catch. A seven-line
# extension whose static initialiser touches std::locale is enough to trigger
# it; OldIsoSpecPy, which the wheel test suite installs, is one such consumer.
# -load_hidden links the archive with all of its symbols hidden, so our private
# copy of the runtime stays private. It replaces -static-libstdc++ (which would
# otherwise pull the same archive in again, visibly).
if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.a
                  OUTPUT_VARIABLE ISOSPEC_LIBSTDCXX_A
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT IS_ABSOLUTE "${ISOSPEC_LIBSTDCXX_A}" OR NOT EXISTS "${ISOSPEC_LIBSTDCXX_A}")
    message(FATAL_ERROR
      "${CMAKE_CXX_COMPILER} could not locate a static libstdc++.a "
      "(-print-file-name returned '${ISOSPEC_LIBSTDCXX_A}'). The macOS wheel "
      "needs one: linking libstdc++ dynamically makes delocate reject the "
      "Homebrew bottle, and linking it visibly breaks other libstdc++ users "
      "in the same process.")
  endif()
  # -load_hidden needs a reasonably modern ld64; an older one rejects the
  # option outright, which is a loud failure rather than a silent regression.
  #
  # The archive is only half of it: our own translation units instantiate std::
  # templates (std::to_string, sorts, the Mersenne twister), and those weak
  # definitions coalesce across images just as readily. std:: is not part of
  # this library's ABI -- callers reach it through cwrapper's extern "C"
  # interface -- so none of it has any business being exported.
  # tests/Python/test_cxx_runtime.py fails if any of this comes undone.
  target_link_options(IsoSpecCppPy PRIVATE
    -static-libgcc "-Wl,-load_hidden,${ISOSPEC_LIBSTDCXX_A}"
    "-Wl,-unexported_symbol,__ZNSt*" "-Wl,-unexported_symbol,__ZNKSt*"
    "-Wl,-unexported_symbol,__ZSt*")
endif()
#python_extension_module(IsoSpecCppPy)
install(TARGETS IsoSpecCppPy LIBRARY DESTINATION IsoSpecPy)
install(DIRECTORY ../src/IsoSpec++/ DESTINATION IsoSpecPy/IsoSpec++
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.hxx"
)

install(FILES ../data/unimod.csv DESTINATION IsoSpecPy/data)
