# IsoSpec C-ABI test suite.
#
# The point of this directory is that everything here is compiled by a *C*
# compiler: cwrapper.h is the header non-C++ callers include, and nothing else
# in the test suite would notice if a C++-only construct leaked into it.  The
# library itself is C++, so the link goes through the C++ driver.
#
#   make        # build the test binary
#   make test   # build and run it (non-zero exit status on any failed check)
#   make clean
#
# test_IsoThresholdGenerator.c / test_IsoOrderedGenerator.c / test_arrays.c are
# older hand-run demo programs (compiled as C++, printing rather than checking);
# they are not part of this build.

CC  ?= cc
CXX ?= c++

CFLAGS   := -std=c11 -Wall -Wextra -Werror -O2
CXXFLAGS := -std=c++20 -O2
INC      := -I../../src/IsoSpec++
UNITY    := ../../src/IsoSpec++/unity-build.cpp
BUILD    := build

.PHONY: all test clean

all: $(BUILD)/test_c_api

$(BUILD)/unity.o: $(UNITY) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $(UNITY) -o $@

$(BUILD)/test_c_api.o: test_c_api.c ../../src/IsoSpec++/cwrapper.h | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c test_c_api.c -o $@

$(BUILD)/test_c_api: $(BUILD)/test_c_api.o $(BUILD)/unity.o
	$(CXX) $^ -o $@ -lm

test: $(BUILD)/test_c_api
	$(BUILD)/test_c_api

$(BUILD):
	mkdir -p $@

clean:
	rm -rf $(BUILD)
