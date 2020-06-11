message("")
message(STATUS "${BoldGreen}Starting configuration for the example material for ${PROJECT_NAME} ${ColourReset}")
message("")


###############
# install stuff

install(DIRECTORY C++  Python  R
	DESTINATION ${CMAKE_INSTALL_FULL_DATAROOTDIR}/doc/isospec/examples)

message("")
message(STATUS "${BoldGreen}Finished configuration for the example material.${ColourReset}")
message("")


