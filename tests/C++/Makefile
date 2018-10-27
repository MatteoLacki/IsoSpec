CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
TESTFLAGS=-fsanitize=leak,undefined
TESTMEMFLAGS= $(TESTFLAGS) -fsanitize=memory
TESTADDRFLAGS= $(TESTFLAGS) -fsanitize=address
LLVMTESTFLAGS= -fsanitize=dataflow,cfi,safe-stack
CXXFLAGS=-std=c++11 -Wall -I../../IsoSpec++ -Wextra -pedantic -lpthread
SRCFILES=../../IsoSpec++/unity-build.cpp

all: formula_layered formula_ordered formula_threshold

formula_ordered:
	clang++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_ordered.cpp -o ./from_formula_ordered_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_ordered.cpp -o ./from_formula_ordered_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_ordered.cpp -o ./from_formula_ordered_debug
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_ordered.cpp -fsanitize=address -o ./from_formula_ordered_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_ordered.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_ordered_memsan


formula_threshold:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold.cpp -o ./from_formula_threshold_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold.cpp -o ./from_formula_threshold_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold.cpp -o ./from_formula_threshold_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold.cpp -fsanitize=address -o ./from_formula_threshold_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_threshold_memsan

formula_layered:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -fsanitize=address -o ./from_formula_layered_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_layered_memsan


test:
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) $(SRCFILES) test.cpp -o ./test-clang
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) $(SRCFILES) test.cpp -o ./test-gcc
	g++-5.4.0 $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) $(SRCFILES) test.cpp -o ./test-gcc-5.4.0
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) test.cpp -o ./test-valgrind

ti:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) ../../IsoSpec++/unity-build.cpp titin-test.cpp -o ./titin

mr:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) ../../IsoSpec++/unity-build.cpp marginal-test.cpp -o marginal -g

tabulator:
	clang++ -std=c++11 ../../IsoSpec++/unity-build.cpp tabulator_test.cpp -o tabulator

la:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) ../../IsoSpec++/unity-build.cpp layered-test.cpp -o layered

IsoThresholdGenerator:
	clang++ -std=c++11 IsoThresholdGenerator.cpp -o IsoThresholdGenerator -lpthread

nr_conf:
	clang++ -std=c++11 nr_conf.cpp -o nr_conf -lpthread
