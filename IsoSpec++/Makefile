OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g -Werror -DISOSPEC_DEBUG -DDEBUG -D_GLIBCXX_DEBUG
CXXFLAGS=-std=c++11 -Wall -pedantic -Wextra -Wshadow -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wno-old-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wno-sign-conversion -Wsign-promo -Wswitch-default -Wundef
SRCFILES=cwrapper.cpp allocator.cpp  dirtyAllocator.cpp  isoSpec++.cpp  isoMath.cpp  marginalTrek++.cpp  operators.cpp element_tables.cpp misc.cpp mman.cpp fixedEnvelopes.cpp fasta.cpp

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
	#g++ -O3 -std=gnu++1z -O3 -shared -static -static-libstdc++ -static-libgcc unity-build.cpp -o ../IsoSpecPy/IsoSpec++.dll
	x86_64-w64-mingw32-g++.exe -std=c++11 -O3 -shared -static -static-libstdc++ unity-build.cpp -o ../IsoSpecPy/IsoSpecPy/prebuilt-libIsoSpec++-2.1.0-x64.dll
	i686-w64-mingw32-g++.exe -std=c++11 -O3 -shared -static -static-libstdc++ unity-build.cpp -o ../IsoSpecPy/IsoSpecPy/prebuilt-libIsoSpec++-2.1.0-x32.dll
check:
	cpplint *.cpp *.h
	#cppcheck --std=c++11 --language=c++ --enable=style,portability,missingInclude,performance -q --force *.cpp *.h --suppress=memleakOnRealloc:*:* --inline-suppr --template={file}:{line}\ \({severity}\)\ \<{id}\>:\ {message}

check11:
	$(CXX) -std=c++11 -Wall -pedantic -Wextra -Wshadow -fsyntax-only unity-build.cpp

asm:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) unity-build.cpp -S -o out.S
