# IsoSpec C++ test suite.
#
# The tests are doctest-based (doctest.h vendored here).  All test_*.cpp plus
# tests_main.cpp and unity-build.cpp link into one binary per configuration:
#
#   run_tests_gcc    release, g++
#   run_tests_clang  release, clang++
#   run_tests_dbg    -O0 -g, _GLIBCXX_DEBUG
#   run_tests_asan   AddressSanitizer + UBSan (clang)
#   run_tests_memsan MemorySanitizer + UBSan (opt-in; needs an MSan-instrumented
#                    libc++, so it is NOT part of `make test`)
#
# Usage:
#   make            # build the default test configs + CLIs
#   make test       # build the default configs and RUN them (fails on any error)
#   make asan       # just the ASan/UBSan config
#   make coverage   # gcov-instrumented run + line/branch coverage report
#   make cli        # build the command-line investigation tools under cli/
#   make clean
#
# Everything is written under build/ (out of tree; git-ignored).

CXX_GCC   ?= g++
CXX_CLANG ?= clang++

STD  := -std=c++20
WARN := -Wall -Wextra -pedantic
INC  := -I../../src/IsoSpec++
UNITY := ../../src/IsoSpec++/unity-build.cpp
BUILD := build

OPTFLAGS   := -O3 -march=native -mtune=native -g
DEBUGFLAGS := -O0 -g -DISOSPEC_DEBUG -DDEBUG -D_GLIBCXX_DEBUG
ASANFLAGS  := -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all
MSANFLAGS  := -O1 -g -fsanitize=memory,undefined -DISOSPEC_TESTS_MEMSAN

# --- runtime ISA dispatch ---------------------------------------------------
# The tests exercise the dispatching build, so each configuration also compiles
# the per-instruction-set-level kernels into their own objects.  Their flags are
# the configuration's minus -march: the baseline level has to stay baseline for
# the cross-level tests to mean anything.  Everything else about the flags must
# match the rest of the configuration -- _GLIBCXX_DEBUG and the sanitizers change
# layouts and interpositions, so mixing them across objects is an ODR violation.
ISA_DIR      := ../../src/IsoSpec++
ISA_DEFS     := -DISOSPEC_ISA_DISPATCH=1
ISA_LEVELS   := baseline avx v3 v4
ISA_OPTFLAGS := -O3 -g

isa-accepts = $(shell echo 'int main(){}' | $(CXX_GCC) $(1) -x c++ - -o /dev/null >/dev/null 2>&1 && echo yes)
AVX_FLAGS := -mavx
V3_FLAGS := -mavx2 -mfma
V4_FLAGS := -mavx512f -mavx512bw -mavx512dq -mavx512vl
ISA_FLAGS_baseline :=
ISA_FLAGS_avx := $(if $(call isa-accepts,$(AVX_FLAGS)),$(AVX_FLAGS))
ISA_FLAGS_v3 := $(if $(call isa-accepts,$(V3_FLAGS)),$(V3_FLAGS))
ISA_FLAGS_v4 := $(if $(call isa-accepts,$(V4_FLAGS)),$(V4_FLAGS))

# doctest suite: tests_main.cpp defines main(); the rest register TEST_CASEs.
# Keep this list alphabetical after tests_main.cpp.
TEST_SRCS := tests_main.cpp \
             test_aligned_ptr.cpp \
             test_analytic.cpp \
             test_api_iso.cpp \
             test_boundary_walk.cpp \
             test_containers.cpp \
             test_cross_generator.cpp \
             test_cwrapper.cpp \
             test_edge_cases.cpp \
             test_fasta.cpp \
             test_fixed_envelope.cpp \
             test_formula.cpp \
             test_generators.cpp \
             test_isa_dispatch.cpp \
             test_marginals.cpp \
             test_math.cpp \
             test_oracle.cpp \
             test_regression.cpp \
             test_simd_equivalence.cpp \
             test_single_atom.cpp \
             test_unimod.cpp \
             test_wasserstein.cpp
# Rebuild whenever anything in the library changes, not just the test sources:
# the unity build pulls in every .cpp, and the headers are full of inline code.
LIB_SRCS  := $(wildcard ../../src/IsoSpec++/*.cpp ../../src/IsoSpec++/*.h)
TEST_DEPS := $(TEST_SRCS) test_helpers.h doctest.h $(UNITY) $(LIB_SRCS)

CLI_SRCS := $(wildcard cli/*.cpp)
CLI_BINS := $(patsubst cli/%.cpp,$(BUILD)/cli/%,$(CLI_SRCS))

.PHONY: all tests test asan memsan coverage cli fasta_header_isolation clean

all: tests cli

# --- test-binary matrix -----------------------------------------------------
# $(1)=config name  $(2)=compiler  $(3)=extra flags
# $(1)=config  $(2)=compiler  $(3)=flags  $(4)=kernel flags (no -march)  $(5)=level
define ISA_rule
$(BUILD)/isa_$(5)_$(1).o: $(ISA_DIR)/isa_kernels_$(5).cpp $(LIB_SRCS) | $(BUILD)
	$(2) $(STD) $(WARN) $(INC) $(4) $(ISA_FLAGS_$(5)) $(ISA_DEFS) -c $$< -o $$@
endef

define TEST_rule
$(foreach lvl,$(ISA_LEVELS),$(eval $(call ISA_rule,$(1),$(2),$(3),$(4),$(lvl))))
ISA_OBJS_$(1) := $(foreach lvl,$(ISA_LEVELS),$(BUILD)/isa_$(lvl)_$(1).o)
$(BUILD)/run_tests_$(1): $(TEST_DEPS) $$(ISA_OBJS_$(1)) | $(BUILD)
	$(2) $(STD) $(WARN) $(INC) $(3) $(ISA_DEFS) $(TEST_SRCS) $(UNITY) $$(ISA_OBJS_$(1)) -o $$@
endef

$(eval $(call TEST_rule,gcc,$(CXX_GCC),$(OPTFLAGS),$(ISA_OPTFLAGS)))
$(eval $(call TEST_rule,clang,$(CXX_CLANG),$(OPTFLAGS),$(ISA_OPTFLAGS)))
$(eval $(call TEST_rule,dbg,$(CXX_GCC),$(DEBUGFLAGS),$(DEBUGFLAGS)))
$(eval $(call TEST_rule,asan,$(CXX_CLANG),$(ASANFLAGS),$(ASANFLAGS)))
$(eval $(call TEST_rule,memsan,$(CXX_CLANG),$(MSANFLAGS),$(MSANFLAGS)))

DEFAULT_TEST_BINS := $(BUILD)/run_tests_gcc $(BUILD)/run_tests_clang \
                     $(BUILD)/run_tests_dbg $(BUILD)/run_tests_asan

tests: $(DEFAULT_TEST_BINS)
asan: $(BUILD)/run_tests_asan
memsan: $(BUILD)/run_tests_memsan

# Build every default config and run it; stop at the first failure.
test: $(DEFAULT_TEST_BINS) fasta_header_isolation
	@set -e; for b in $(DEFAULT_TEST_BINS); do \
		echo "=================== running $$b ==================="; \
		$$b; \
	done; \
	echo "ALL C++ TESTS PASSED"

# --- coverage ---------------------------------------------------------------
# Line/branch coverage of src/IsoSpec++ under the doctest suite.  Built at -O0
# with gcov instrumentation into its own directory, so it never collides with
# the ordinary build.  Needs gcovr (pip install gcovr); without it the raw
# .gcda files are still produced and any gcov front-end can read them.
COVFLAGS := -O0 -g --coverage -fprofile-abs-path
COVDIR   := $(BUILD)/coverage

$(foreach lvl,$(ISA_LEVELS),$(eval $(call ISA_rule,cov,$(CXX_GCC),$(COVFLAGS),$(COVFLAGS),$(lvl))))
COV_ISA_OBJS := $(foreach lvl,$(ISA_LEVELS),$(BUILD)/isa_$(lvl)_cov.o)

$(COVDIR)/run_tests_cov: $(TEST_DEPS) $(COV_ISA_OBJS) | $(COVDIR)
	$(CXX_GCC) $(STD) $(WARN) $(INC) $(COVFLAGS) $(ISA_DEFS) $(TEST_SRCS) $(UNITY) $(COV_ISA_OBJS) -o $@

coverage: $(COVDIR)/run_tests_cov
	@rm -f $(COVDIR)/*.gcda
	$(COVDIR)/run_tests_cov
	@echo
	gcovr -r ../../src/IsoSpec++ --exclude-unreachable-branches \
	      --txt-metric line --txt $(COVDIR)/coverage.txt \
	      --html-details $(COVDIR)/coverage.html $(COVDIR)
	@cat $(COVDIR)/coverage.txt
	@echo "HTML report: $(COVDIR)/coverage.html"

# --- compile-only header-isolation check ------------------------------------
# fasta.h must be self-contained (it uses memset/size_t); compiling it in a TU
# that includes nothing else must succeed.
fasta_header_isolation: fasta_header_isolation.cpp | $(BUILD)
	$(CXX_CLANG) $(STD) $(INC) -c fasta_header_isolation.cpp -o $(BUILD)/fasta_header_isolation.o

# --- command-line investigation tools ---------------------------------------
cli: $(CLI_BINS)

$(BUILD)/cli/%: cli/%.cpp $(UNITY) | $(BUILD)/cli
	$(CXX_CLANG) $(STD) $(WARN) $(INC) $(OPTFLAGS) $(UNITY) $< -o $@

$(BUILD) $(BUILD)/cli $(COVDIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD)
