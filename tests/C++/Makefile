CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
CXXFLAGS=-std=c++11 -Wall -I../../IsoSpec++
SRCFILES=`ls ../../IsoSpec++/*.cpp`

all: testcc


testcc:
	clang++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula.cpp -o ./from_formula_clang
	g++ $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula.cpp -o ./from_formula_gcc
	g++-5.4.0 $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) from_formula.cpp -o ./from_formula_gcc-5.4.0
