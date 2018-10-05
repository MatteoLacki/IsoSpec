CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
TESTFLAGS=-fsanitize=leak,undefined
TESTMEMFLAGS= $(TESTFLAGS) -fsanitize=memory
TESTADDRFLAGS= $(TESTFLAGS) -fsanitize=address
LLVMTESTFLAGS= -fsanitize=dataflow,cfi,safe-stack
CXXFLAGS=-std=c++11 -Wall -I../../IsoSpec++ -Wextra -pedantic -lpthread
SHARED_FLAGS=-L../../IsoSpec++
LINK_FLAGS=-Wl,-rpath,../../IsoSpec++ ../../IsoSpec++/libIsoSpec++.so

SRCFILES=../../IsoSpec++/unity-build.cpp

all: testcc


testcc:
	clang++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula.cpp -o ./from_formula_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula.cpp -o ./from_formula_gcc
	g++-5 $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula.cpp -o ./from_formula_gcc-5

simpletest-debug:
	clang++ $(SHARED_FLAGS) $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) -lIsoSpec++ simpletest.cpp -o ./simpletest-clang-debug $(LINK_FLAGS)
	g++ $(SHARED_FLAGS) $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) -lIsoSpec++ simpletest.cpp -o ./simpletest-gcc-debug $(LINK_FLAGS)

simpletest:
	clang++ $(SHARED_FLAGS) $(CXXFLAGS) -lIsoSpec++ simpletest.cpp -o ./simpletest-clang $(LINK_FLAGS)
	g++ $(SHARED_FLAGS) $(CXXFLAGS) -pthread -lIsoSpec++ simpletest.cpp -o ./simpletest-gcc $(LINK_FLAGS)

segfault:
	clang++ $(SHARED_FLAGS) $(CXXFLAGS) -lIsoSpec++ segfault.cpp -o ./segfault-clang $(LINK_FLAGS)
	clang++ $(SHARED_FLAGS) $(CXXFLAGS) -DPREVENT_SEGFAULT -lIsoSpec++ segfault.cpp -o ./segfault-ok-clang $(LINK_FLAGS)
	clang++ $(SHARED_FLAGS) $(CXXFLAGS) -fsanitize=address -lIsoSpec++ segfault.cpp -o ./segfault-ok2-clang $(LINK_FLAGS)

test:
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) $(SRCFILES) test.cpp -o ./test-clang
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) $(SRCFILES) test.cpp -o ./test-gcc
	g++-5.4.0 $(CXXFLAGS) $(DEBUGFLAGS) $(TESTADDRFLAGS) $(SRCFILES) test.cpp -o ./test-gcc-5.4.0
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) test.cpp -o ./test-valgrind

ti:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) ../../IsoSpec++/unity-build.cpp titin-test.cpp -o ./titin

mt:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) ../../IsoSpec++/unity-build.cpp titin-multithreaded.cpp -o ./multithreaded -lpthread

mr:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) ../../IsoSpec++/unity-build.cpp marginal-test.cpp -o marginal -g

tabulator:
	clang++ -std=c++11 ../../IsoSpec++/unity-build.cpp tabulator_test.cpp -o tabulator

la:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) ../../IsoSpec++/unity-build.cpp layered-test.cpp -o layered

IsoThresholdGenerator:
	clang++ -std=c++11 IsoThresholdGenerator.cpp -o IsoThresholdGenerator -lpthread
