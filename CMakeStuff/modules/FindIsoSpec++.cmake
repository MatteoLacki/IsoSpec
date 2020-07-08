# Copyright : Filippo Rusconi
# License : GPL-3.0+
# Authors : Filippo Rusconi

find_path(IsoSpec++_INCLUDE_DIRS IsoSpec++/isoSpec++.h
	PATHS /usr/local/include /usr/include
	PATH_SUFFIXES isospec++ libisospec++ ENV PATH)

find_library(IsoSpec++_LIBRARY NAMES isospec libisospec isospec++ libisospec++)

if(IsoSpec++_INCLUDE_DIRS AND IsoSpec++_LIBRARY)

	message(STATUS "Found Found IsoSpec++_LIBRARY: ${IsoSpec++_LIBRARY}")

	set(IsoSpec++_FOUND TRUE)

endif()

if(IsoSpec++_FOUND)

	if (NOT IsoSpec++_FIND_QUIETLY)

		message(STATUS "Found IsoSpec++_LIBRARY: ${IsoSpec++_LIBRARY}")

	endif()

	if(NOT TARGET IsoSpec++::IsoSpec++)

		add_library(IsoSpec++::IsoSpec++ UNKNOWN IMPORTED)

		set_target_properties(IsoSpec++::IsoSpec++ PROPERTIES
			IMPORTED_LOCATION             "${IsoSpec++_LIBRARY}"
			INTERFACE_INCLUDE_DIRECTORIES "${IsoSpec++_INCLUDE_DIRS}")

	endif()

endif()

