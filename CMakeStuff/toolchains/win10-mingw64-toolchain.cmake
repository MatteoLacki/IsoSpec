message("\n${BoldRed}WIN10-MINGW64 environment${ColourReset}\n")
message("Please run the configuration like this:")
message("cmake -G \"Unix Makefiles\" -DCMAKE_BUILD_TYPE=Release ../development")


if(WIN32 OR _WIN32)
	message(STATUS "Building with WIN32 defined.")
endif()


# On Win10 all the code is relocatable.
message(STATUS "Removing definitions -fPIC.")
remove_definitions(-fPIC)

