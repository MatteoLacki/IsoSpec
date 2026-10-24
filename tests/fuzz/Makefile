# libFuzzer targets for everything in the library that parses text:
#
#   fuzz_formula   a chemical formula, as Iso(const char*) takes it
#   fuzz_sequence  a peptide sequence, plain and with [UNIMOD:<id>] modifications
#   fuzz_unimod    an override Unimod CSV, then a sequence parsed against it
#
# clang only (libFuzzer), built with ASan + UBSan, every report fatal.
#
# Usage:
#   make                      # build the three targets under build/
#   make run TIME=600 JOBS=8  # fuzz each for TIME seconds with JOBS workers,
#                             # starting from seeds/<target>/; the growing
#                             # corpus and any crash files go under build/
#   make clean
#
# Not part of tests/C++'s `make test`: a fuzzing run has no natural end, and
# its value is in running long, which is a job for a person, not a check.

# make's built-in default CXX is g++, which has no libFuzzer.
ifeq ($(origin CXX),default)
CXX := clang++
endif
FLAGS := -std=c++20 -O1 -g -fsanitize=fuzzer,address,undefined -fno-sanitize-recover=all
INC   := -I../../src/IsoSpec++
UNITY := ../../src/IsoSpec++/unity-build.cpp
BUILD := build

TARGETS := fuzz_formula fuzz_sequence fuzz_unimod
LIB_SRCS := $(wildcard ../../src/IsoSpec++/*.cpp ../../src/IsoSpec++/*.h)

TIME ?= 60
JOBS ?= 4

.PHONY: all run clean

all: $(addprefix $(BUILD)/,$(TARGETS))

$(BUILD):
	mkdir -p $@

# The library is compiled once, instrumented for coverage but without
# libFuzzer's main() (-fsanitize=fuzzer-no-link); each target links it.
$(BUILD)/unity.o: $(UNITY) $(LIB_SRCS) | $(BUILD)
	$(CXX) $(subst fuzzer,fuzzer-no-link,$(FLAGS)) $(INC) -c $< -o $@

$(BUILD)/%: %.cpp $(BUILD)/unity.o
	$(CXX) $(FLAGS) $(INC) $< $(BUILD)/unity.o -o $@

# -jobs rather than -fork: libFuzzer's fork mode crashed in its own
# secondsSinceProcessStartUp() under clang 23. Each job writes fuzz-<n>.log into
# the directory it runs in, hence the cd into build/logs/<target>.
run: all
	for t in $(TARGETS); do \
	    mkdir -p $(BUILD)/corpus/$$t $(BUILD)/crashes/$$t $(BUILD)/logs/$$t; \
	    (cd $(BUILD)/logs/$$t && $(CURDIR)/$(BUILD)/$$t $(CURDIR)/$(BUILD)/corpus/$$t $(CURDIR)/seeds/$$t \
	        -max_total_time=$(TIME) -jobs=$(JOBS) -workers=$(JOBS) \
	        -artifact_prefix=$(CURDIR)/$(BUILD)/crashes/$$t/ -print_final_stats=1) || exit 1; \
	done

clean:
	rm -rf $(BUILD)
