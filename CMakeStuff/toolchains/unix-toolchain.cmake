message("")
message("${BoldRed}UNIX non APPLE environment${ColourReset}")
message("")
message("~~~~~~ Instructions ~~~~~~")
message("On UNIX/GNU-Linux, please run the configuration like this:")
message("cmake -G \"Unix Makefiles\" -DCMAKE_BUILD_TYPE=<Debug | Release> -DCMAKE_INSTALL_PREFIX=</usr | your_dir> ../development")
message("")

set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES /usr/include)
set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES /usr/include)

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)

if (WITH_FPIC)
	add_definitions(-fPIC)
endif()

# Install cmake module
install(FILES ${CMAKE_MODULE_PATH}/FindIsoSpec++.cmake 
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/isospec++)

# Install cmake config
configure_file (${CMAKE_MODULE_PATH}/IsoSpec++Config.cmake.in
	${CMAKE_BINARY_DIR}/IsoSpec++Config.cmake)
install(FILES ${CMAKE_BINARY_DIR}/IsoSpec++Config.cmake 
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/isospec++)

# Install the PkgConfig config file
configure_file (${CMAKE_MODULE_PATH}/pkgconfig/libisospec++.pc.in
	${CMAKE_BINARY_DIR}/libisospec++.pc)
install(FILES ${CMAKE_BINARY_DIR}/libisospec++.pc 
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Prepare the logo picture files with the right version (configure_file)

configure_file(${CMAKE_SOURCE_DIR}/CMakeStuff/isospec_logo2_high.svg.in
	${CMAKE_SOURCE_DIR}/man/images/isospec_logo2_high_versioned.svg @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/CMakeStuff/isospec_logo2_long.svg.in
	${CMAKE_SOURCE_DIR}/man/images/isospec_logo2_long_versioned.svg @ONLY)

# Make the conversion of the svg file into a png, but only on GNU/Linux
# Produce a file with respected aspect ratio, 200 pixels wide.

if(UNIX AND NOT APPLE)
	execute_process(COMMAND gm convert -geometry 200x
		${CMAKE_SOURCE_DIR}/man/images/isospec_logo2_high_versioned.svg
		${CMAKE_SOURCE_DIR}/man/images/isospec_logo2_high_versioned.png)
endif()


if(UNIX AND NOT APPLE)
	execute_process(COMMAND gm convert -geometry 200x
		${CMAKE_SOURCE_DIR}/man/images/isospec_logo2_long_versioned.svg
		${CMAKE_SOURCE_DIR}/man/images/isospec_logo2_long_versioned.png)
endif()

# Ensure that the doxyfile configuration file for Doxygen has always
# the proper version number!

configure_file(${CMAKE_SOURCE_DIR}/CMakeStuff/doxyfile.in
	${CMAKE_SOURCE_DIR}/man/doxyfile @ONLY)

# Command:
# make doc
add_custom_target(doc
	COMMAND doxygen man/doxyfile 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Doxygen-based developer documentation generation")


