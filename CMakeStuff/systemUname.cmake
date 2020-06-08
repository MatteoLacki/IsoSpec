# Ask that uname -s be performed and store the value in SYSTEM_UNAME_S for
# later reference.

macro(get_uname_string)

execute_process(COMMAND uname -s OUTPUT_VARIABLE SYSTEM_UNAME_S) 

if(${SYSTEM_UNAME_S} MATCHES "^.*MINGW64.*")
	message(STATUS "System detected as Windows, setting WIN64")
	set(WIN64 1)
#else()
	#message(STATUS "System is not Windows.")
endif()

endmacro()

get_uname_string()

