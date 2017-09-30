CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
CXXFLAGS=-std=c++11 -Wall -pedantic -Wextra
SRCFILES=cwrapper.cpp allocator.cpp  dirtyAllocator.cpp  isoSpec++.cpp  isoMath.cpp  marginalTrek++.cpp  operators.cpp element_tables.cpp misc.cpp

all: unitylib

lib:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) -fPIC -shared -o libIsoSpec++.so 

unitylib:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) unity-build.cpp -fPIC -shared -o libIsoSpec++.so

debug:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) unity-build.cpp -DDEBUG -fPIC -shared -o libIsoSpec++.so

perftest:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) -DDEBUG -fPIC -shared -o libIsoSpec++.so 

ctest: lib
	$(CC) -g -L. -lIsoSpec++ test.c -o ctest

tests: lib
	$(CXX) -std=c++11 -I. -O3 -L. -Wl,-rpath,. -lIsoSpec++ test2.cpp -o test2
	$(CXX) -std=c++11 -I. -O3 -L. -Wl,-rpath,. -lIsoSpec++ test3.cpp -o test3
clean:
	rm -f libIsoSpec++.so

windows:
	g++ -O3 -std=gnu++11 -O3 -shared -static -static-libstdc++ -static-libgcc *.cpp -o lib.dll

unity:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) unity-build.cpp -fPIC -shared -o libIsoSpec++.so

