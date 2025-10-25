# CMake script for IsoSpec++ library
# Author: Filippo Rusconi
# Created: June 2020

message("\n${BoldRed}Now configuring IsoSpec++/ for ${CMAKE_PROJECT_NAME}${ColourReset}\n")


set(isospec_SRCS
	unity-build.cpp
	)

######### The shared library #########

add_library(IsoSpec++-shared SHARED ${isospec_SRCS})

message("Library version: ${${PROJECT_NAME}_VERSION}")
message("Library soname: ${${PROJECT_NAME}_LIB_SOVERSION}")

add_definitions(-fPIC)

set_target_properties(IsoSpec++-shared
	PROPERTIES OUTPUT_NAME IsoSpec++
	CLEAN_DIRECT_OUTPUT 1
	VERSION  ${${PROJECT_NAME}_VERSION}
	SOVERSION ${${PROJECT_NAME}_LIB_SOVERSION}
	)


######### The static library #########

add_library(IsoSpec++-static STATIC ${isospec_SRCS})

set_target_properties(IsoSpec++-static
	PROPERTIES OUTPUT_NAME IsoSpec++
	CLEAN_DIRECT_OUTPUT 1
	)


# Install both shared and static libs

install(TARGETS IsoSpec++-shared IsoSpec++-static
	LIBRARY NAMELINK_SKIP
	DESTINATION ${CMAKE_INSTALL_LIBDIR}
	# This one is for WIN32
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	)

# Install the header files
install(
    DIRECTORY IsoSpec++
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
)
