message("")
message("${BoldRed}APPLE environment${ColourReset}")
message("")
message("~~~~~~ Instructions ~~~~~~")
message("cmake -DCMAKE_BUILD_TYPE=<Debug | Release> ../../development")
message("")

set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/opt/local/include")
set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/opt/local/lib")


set(HOME_DEVEL_DIR "/Users/rusconi/devel")

set(CMAKE_MACOSX_RPATH 0)

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)

if(WITH_FPIC)
	add_definitions(-fPIC)
endif()


