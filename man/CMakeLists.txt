
message("")
message(STATUS "${BoldGreen}Starting configuration for the doc material for ${PROJECT_NAME} ${ColourReset}")
message("")

# We want to build the documentation in the binary directory, not in the
# source directory (mind the ${CMAKE_CURRENT_BINARY_DIR}).

# Prepare the logo picture files with the right version (configure_file)
configure_file(${CMAKE_SOURCE_DIR}/CMakeStuff/isospec_logo2_high.svg.in
	${CMAKE_CURRENT_BINARY_DIR}/images/isospec_logo2_high_versioned.svg @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/CMakeStuff/isospec_logo2_long.svg.in
	${CMAKE_CURRENT_BINARY_DIR}/images/isospec_logo2_long_versioned.svg @ONLY)

# Make the conversion of the svg file into a png, but only on GNU/Linux
# Produce a file with respected aspect ratio, 200 pixels wide.

if(UNIX AND NOT APPLE)

	execute_process(COMMAND gm convert -geometry 200x
		${CMAKE_CURRENT_BINARY_DIR}/images/isospec_logo2_high_versioned.svg
		${CMAKE_CURRENT_BINARY_DIR}/images/isospec_logo2_high_versioned.png)

		execute_process(COMMAND gm convert -geometry 200x
		${CMAKE_CURRENT_BINARY_DIR}/images/isospec_logo2_long_versioned.svg
		${CMAKE_CURRENT_BINARY_DIR}/images/isospec_logo2_long_versioned.png)

endif()

# Ensure that the doxyfile configuration file for Doxygen has always
# the proper version number!

configure_file(${CMAKE_SOURCE_DIR}/CMakeStuff/doxyfile.in
	${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)

add_custom_target(html_doc
	COMMAND doxygen doxyfile
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Doxygen-based developer html documentation generation")

add_custom_target(pdf_doc
	COMMAND make
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex
	COMMENT "Doxygen-based developer pdf documentation generation")


###############
# install stuff

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/images
	DESTINATION ${CMAKE_INSTALL_FULL_DATAROOTDIR}/doc/isospec)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
	DESTINATION ${CMAKE_INSTALL_FULL_DATAROOTDIR}/doc/isospec)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/latex/refman.pdf
	DESTINATION ${CMAKE_INSTALL_FULL_DATAROOTDIR}/doc/isospec/pdf)

file(GLOB man3_pages "${CMAKE_CURRENT_BINARY_DIR}/man/man3/*")

install(FILES ${man3_pages}
	DESTINATION ${CMAKE_INSTALL_FULL_MANDIR}/man3)

message("")
message(STATUS "${BoldGreen}Finished configuration of the doc material.${ColourReset}")
message("")

