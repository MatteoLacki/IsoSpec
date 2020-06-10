# CMake script for IsoSpec++ library
# Author: Filippo Rusconi
# Created: June 2020

message(\n${BoldRed}"Now configuring IsoSpec++/ for ${CMAKE_PROJECT_NAME}"${ColourReset}\n)


set(isospec_SRCS
	unity-build.cpp
	)


# Not sure I understand this. Is this the way to test for MSVS in CMake land ?
if(WIN32 OR _WIN32)

	if(MSVC)
		set(my_sources ${my_sources})
		add_definitions(-DMMAN_LIBRARY)
	endif()

endif()


#message(STATUS "ISOSPEC_LIB_VERSION: ${ISOSPEC_LIB_VERSION}")
#message(STATUS "ISOSPEC_LIB_SOVERSION: ${ISOSPEC_LIB_SOVERSION}")


######### The shared library #########

add_library(IsoSpec++-shared SHARED ${isospec_SRCS})

set_target_properties(IsoSpec++-shared
	PROPERTIES OUTPUT_NAME IsoSpec++
	CLEAN_DIRECT_OUTPUT 1
	VERSION ${ISOSPEC_LIB_VERSION}
	SOVERSION ${ISOSPEC_LIB_SOVERSION}
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


# Do not know how this would still be useful.
# export(TARGETS IsoSpec++ FILE IsoSpecLibrary.cmake)

