# This file should be included if the command line reads like this:
# x86_64-w64-mingw32.shared-cmake -DCMAKE_BUILD_TYPE=Release -DMXE=1 ..

message("MXE (M cross environment) https://mxe.cc/")
message("Please run the configuration like this:")
message("x86_64-w64-mingw32.shared-cmake -DMXE=1 -G \"Unix Makefiles\" -DCMAKE_BUILD_TYPE=Release ../../development")

set(HOME_DEVEL_DIR "/home/rusconi/devel")

set(MXE_ROOT_DIR "${HOME_DEVEL_DIR}/mxe")
set(MXE_SHIPPED_DLLS_DIR "${MXE_ROOT_DIR}/dll-set-for-packages")

set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES ${HOME_DEVEL_DIR}/mxe/usr/x86_64-w64-mingw32.shared/include)
set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES ${HOME_DEVEL_DIR}/mxe/usr/x86_64-w64-mingw32.shared/include)

if(WIN32 OR _WIN32)
    message(STATUS "Building with WIN32 defined.")
endif()

## We can build the package setup executable with this specific command.
add_custom_target(
    dllinstall
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/IsoSpec++/libIsoSpec++.dll ${MXE_SHIPPED_DLLS_DIR}
    COMMENT "Build and copy the dll file to its MXE destination."
    DEPENDS ${CMAKE_BINARY_DIR}/IsoSpec++/libIsoSpec++.dll
    VERBATIM
)
