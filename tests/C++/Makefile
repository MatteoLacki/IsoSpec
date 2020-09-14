CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native -g
DEBUGFLAGS=-O0 -g -DISOSPEC_DEBUG -DDEBUG -D_GLIBCXX_DEBUG
TESTFLAGS=-fsanitize=leak,undefined
TESTMEMFLAGS= $(TESTFLAGS) -fsanitize=memory
TESTADDRFLAGS= $(TESTFLAGS) -fsanitize=address
LLVMTESTFLAGS= -fsanitize=dataflow,cfi,safe-stack
CXXFLAGS=-std=c++17 -Wall -I../../IsoSpec++ -Wextra -pedantic
SRCFILES=../../IsoSpec++/unity-build.cpp

all: main_test cmdlines

main_test:
	clang++ $(CXXFLAGS) $(OPTFLAGS) main_test.cpp -o main_test_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) main_test.cpp -o main_test_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -o main_test_dbg
	g++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -DISOSPEC_SKIP_SLOW_TESTS -o main_test_dbg_fast
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -DISOSPEC_SKIP_SLOW_TESTS -fsanitize=address,undefined -o main_test_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -fsanitize=cfi -flto -fvisibility=hidden -o main_test_cfi
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -fsanitize=safe-stack -o main_test_ss
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) main_test.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o main_test_memsan

cmdlines: formula_layered formula_layered_generator formula_ordered formula_threshold formula_threshold_simple formula_threshold_profile mass_range formula_stochastic

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

formula_threshold_profile:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold_profile.cpp -g -pg -o ./from_formula_threshold_gprof_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold_profile.cpp -g -pg -o ./from_formula_threshold_gprof_gcc
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold_profile.cpp -g -o ./from_formula_threshold_profile_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_threshold_profile.cpp -g -o ./from_formula_threshold_profile_gcc

formula_stochastic:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_stochastic.cpp -o ./from_formula_stochastic_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_stochastic.cpp -o ./from_formula_stochastic_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_stochastic.cpp -o ./from_formula_stochastic_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_stochastic.cpp -fsanitize=address,undefined -o ./from_formula_stochastic_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_stochastic.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_stochastic_memsan

formula_layered:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -o ./from_formula_layered_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -fsanitize=address,undefined -o ./from_formula_layered_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_layered_memsan

formula_layered_generator:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered_generator.cpp -o ./from_formula_layered_generator_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula_layered_generator.cpp -o ./from_formula_layered_generator_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered_generator.cpp -o ./from_formula_layered_generator_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered_generator.cpp -fsanitize=address,undefined -o ./from_formula_layered_generator_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) from_formula_layered_generator.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./from_formula_layered_generator_memsan

mass_range:
	clang++  $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) mass_range.cpp -o ./mass_range_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) mass_range.cpp -o ./mass_range_gcc
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) mass_range.cpp -o ./mass_range_dbg
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) mass_range.cpp -fsanitize=address,undefined -o ./mass_range_asan
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) mass_range.cpp -DISOSPEC_TESTS_MEMSAN -fsanitize=memory,undefined -o ./mass_range_memsan


ti:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) ../../IsoSpec++/unity-build.cpp titin-test.cpp -o ./titin

mr:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) marginal-test.cpp -o marginal -g

tabulator:
	clang++ -std=c++11 tabulator_test.cpp -o tabulator

la:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) ../../IsoSpec++/unity-build.cpp layered-test.cpp -o layered

IsoThresholdGenerator:
	clang++ -std=c++11 IsoThresholdGenerator.cpp -o IsoThresholdGenerator -lpthread

nr_conf:
	clang++ -std=c++11 nr_conf.cpp -o nr_conf -lpthread

clean:
	rm -rf *_gcc *_clang *_dbg *_memsan *_asan IsoThresholdGenerator layered main_test_cfi main_test_ss marginal nr_conf tabulator titin *.dSYM a.out from_formula_layered_clang_more main_test_dbg_fast
