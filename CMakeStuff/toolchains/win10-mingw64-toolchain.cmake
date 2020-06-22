message("\n${BoldRed}WIN10-MINGW64 environment${ColourReset}\n")
message("Please run the configuration like this:")
message("cmake -G \"Unix Makefiles\" -DCMAKE_BUILD_TYPE=Release ../development")


if(WIN32 OR _WIN32)
	message(STATUS "Building with WIN32 defined.")
endif()


message(STATUS "${BoldGreen}Setting definition -DISOSPEC_MAKE_DLL for symbol DLL export.${ColourReset}")
add_definitions(-DISOSPEC_MAKE_DLL)


# On Win10 all the code is relocatable.
remove_definitions(-fPIC -Wall -pedantic -Wextra)

