CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
CXXFLAGS=-std=c++11 -Wall
SRCFILES=cwrapper.cpp allocator.cpp  dirtyAllocator.cpp  isoSpec++.cpp  isoMath.cpp  marginalTrek++.cpp  operators.cpp element_tables.cpp misc.cpp

all: lib

lib:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) -fPIC -shared -o libIsoSpec++.so 

debug:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) -DDEBUG -fPIC -shared -o libIsoSpec++.so

perftest:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) -DDEBUG -fPIC -shared -o libIsoSpec++.so 

ctest: lib
	$(CC) -g -L. -lIsoSpec++ test.c -o ctest

tests: lib
	$(CXX) -std=c++11 -I. -O3 -L. -Wl,-rpath,. -lIsoSpec++ test2.cpp -o test2
	$(CXX) -std=c++11 -I. -O3 -L. -Wl,-rpath,. -lIsoSpec++ test3.cpp -o test3
clean:
	rm -f libIsoSpec++.so
