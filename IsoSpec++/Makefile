OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g -Werror
CXXFLAGS=-std=c++17 -Wall -pedantic -Wextra
SRCFILES=cwrapper.cpp allocator.cpp  dirtyAllocator.cpp  isoSpec++.cpp  isoMath.cpp  marginalTrek++.cpp  operators.cpp element_tables.cpp misc.cpp mman.cpp precalculatedEnvelopes.cpp

all: unitylib

unitylib:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) unity-build.cpp -fPIC -shared -o libIsoSpec++.so

debug:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) unity-build.cpp -DDEBUG -fPIC -shared -o libIsoSpec++.so

debug-gcc:
	g++ $(CXXFLAGS) $(DEBUGFLAGS) unity-build.cpp -DDEBUG -fPIC -shared -o libIsoSpec++.so

nonunity:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) -DDEBUG -fPIC -shared -o libIsoSpec++.so

clean:
	rm -f libIsoSpec++.so

windows:
	g++ -O3 -std=gnu++11 -O3 -shared -static -static-libstdc++ -static-libgcc unity-build.cpp -o ../IsoSpecPy/IsoSpec++.dll
	x86_64-w64-mingw32-g++.exe -std=gnu++11 -O3 -shared -static -static-libstdc++ unity-build.cpp -o ../IsoSpecPy/IsoSpecPy/prebuilt-libIsoSpec++1.9.3-x64.dll
	i686-w64-mingw32-g++.exe -std=gnu++11 -O3 -shared -static -static-libstdc++ unity-build.cpp -o ../IsoSpecPy/IsoSpecPy/prebuilt-libIsoSpec++1.9.3-x32.dll
