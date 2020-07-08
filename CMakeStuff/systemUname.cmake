# Ask that uname -s be performed and store the value in SYSTEM_UNAME_S for
# later reference.

macro(get_uname_string)

execute_process(COMMAND uname -s OUTPUT_VARIABLE SYSTEM_UNAME_S) 

if(${SYSTEM_UNAME_S} MATCHES "MINGW64_NT-10.*")
    message(STATUS "System detected as Windows10 with MINGW64, setting WIN32 AND WIN10MINGW64")
	# Note that WIN32 is set even on 64 bits systems.
	set(WIN32 1)
	set(WIN10MINGW64 1)
#else()
	#message(STATUS "System is not Windows.")
endif()

endmacro()

get_uname_string()

