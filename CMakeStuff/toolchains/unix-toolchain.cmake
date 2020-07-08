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

if(WITH_FPIC)
	add_definitions(-fPIC)
endif()

# Install cmake module
install(FILES ${CMAKE_MODULE_PATH}/FindIsoSpec++.cmake 
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/isospec++)

# Install cmake config
configure_file(${CMAKE_MODULE_PATH}/IsoSpec++Config.cmake.in
	${CMAKE_BINARY_DIR}/IsoSpec++Config.cmake)
install(FILES ${CMAKE_BINARY_DIR}/IsoSpec++Config.cmake 
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/isospec++)

# Install the PkgConfig config file (only substitute the @VAR@
# because we need to preserve all the ${prefix} strings like
# they are.
configure_file(${CMAKE_MODULE_PATH}/pkgconfig/libisospec++.pc.in
	${CMAKE_BINARY_DIR}/libisospec++.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/libisospec++.pc 
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Now deal with the manual doc stuff
add_subdirectory(man)

# Now deal with the example  stuff
add_subdirectory(Examples)

