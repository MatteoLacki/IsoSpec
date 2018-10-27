CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
TESTFLAGS=-fsanitize=leak,undefined
TESTMEMFLAGS= $(TESTFLAGS) -fsanitize=memory
TESTADDRFLAGS= $(TESTFLAGS) -fsanitize=address
LLVMTESTFLAGS= -fsanitize=dataflow,cfi,safe-stack
CXXFLAGS=-std=c++11 -Wall -I../../IsoSpec++ -Wextra -pedantic
SRCFILES=../../IsoSpec++/unity-build.cpp

all: main_test

main_test:
	clang++ $(CXXFLAGS) $(OPTFLAGS) main_test.cpp -o main_test_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) main_test.cpp -o main_test_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -o main_test_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -fsanitize=address,undefined -o main_test_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -fsanitize=cfi -flto -fvisibility=hidden -o main_test_cfi
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -fsanitize=safe-stack -o main_test_ss
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o main_test_memsan

cmdlines: formula_layered formula_ordered formula_threshold formula_threshold_simple

formula_ordered:
	clang++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_ordered.cpp -o ./from_formula_ordered_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_ordered.cpp -o ./from_formula_ordered_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_ordered.cpp -o ./from_formula_ordered_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_ordered.cpp -fsanitize=address,undefined -o ./from_formula_ordered_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_ordered.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_ordered_memsan

formula_threshold_simple:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold_simple.cpp -o ./from_formula_threshold_simple_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold_simple.cpp -o ./from_formula_threshold_simple_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold_simple.cpp -o ./from_formula_threshold_simple_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold_simple.cpp -fsanitize=address,undefined -o ./from_formula_threshold_simple_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold_simple.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_threshold_simple_memsan

formula_threshold:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold.cpp -o ./from_formula_threshold_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold.cpp -o ./from_formula_threshold_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold.cpp -o ./from_formula_threshold_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold.cpp -fsanitize=address,undefined -o ./from_formula_threshold_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_threshold.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_threshold_memsan

formula_layered:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -fsanitize=address,undefined -o ./from_formula_layered_asan
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

clean:
	rm -f *_gcc *_clang *_dbg *_memsan *_asan
