OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g -Werror -DISOSPEC_DEBUG -DDEBUG -D_GLIBCXX_DEBUG
CXXFLAGS=-std=c++20 -Wall -pedantic -Wextra -Wshadow -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wno-old-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wno-sign-conversion -Wsign-promo -Wswitch-default -Wundef
SRCFILES=cwrapper.cpp allocator.cpp  dirtyAllocator.cpp  isoSpec++.cpp  isoMath.cpp  marginalTrek++.cpp  operators.cpp element_tables.cpp misc.cpp mman.cpp fixedEnvelopes.cpp fasta.cpp element_lookup.cpp unimod.cpp fasta_mods.cpp

# ---------------------------------------------------------------------------
# Runtime ISA dispatch
#
# The two batched fill kernels are compiled once per instruction-set level, each
# into its own translation unit with its own -march, and the one to run is chosen
# from the CPU at startup. They must never share a compilation with the rest of
# the library: a raised -march would otherwise reach every inline function in the
# headers, and the linker could then pick an AVX-512 copy of one for a code path
# reached on a CPU without it. See isa_kernels.h for the full reasoning.
#
# Note the kernels deliberately do NOT get $(OPTFLAGS): -march=native there would
# defeat the point of having a baseline level at all.
# ---------------------------------------------------------------------------
ISA_DEFS=-DISOSPEC_ISA_DISPATCH=1
ISA_OPTFLAGS=-O3
ISA_OBJS=isa_kernels_baseline.o isa_kernels_avx.o isa_kernels_v3.o isa_kernels_v4.o

# Offer each level its flags only where the compiler will take them (they are
# meaningless off x86). The .cpp files then self-check the predefined macros the
# flags produce, so a level whose flags did not arrive is simply absent rather
# than quietly built at the wrong width.
cxx-accepts = $(shell echo 'int main(){}' | $(CXX) $(1) -x c++ - -o /dev/null >/dev/null 2>&1 && echo yes)
AVX_FLAGS=-mavx
V3_FLAGS=-mavx2 -mfma
V4_FLAGS=-mavx512f -mavx512bw -mavx512dq -mavx512vl
ISA_AVX_FLAGS:=$(if $(call cxx-accepts,$(AVX_FLAGS)),$(AVX_FLAGS))
ISA_V3_FLAGS:=$(if $(call cxx-accepts,$(V3_FLAGS)),$(V3_FLAGS))
ISA_V4_FLAGS:=$(if $(call cxx-accepts,$(V4_FLAGS)),$(V4_FLAGS))

all: unitylib

isa_kernels_baseline.o: isa_kernels_baseline.cpp
	$(CXX) $(CXXFLAGS) $(ISA_OPTFLAGS) $(ISA_DEFS) -fPIC -c $< -o $@

isa_kernels_avx.o: isa_kernels_avx.cpp
	$(CXX) $(CXXFLAGS) $(ISA_OPTFLAGS) $(ISA_AVX_FLAGS) $(ISA_DEFS) -fPIC -c $< -o $@

isa_kernels_v3.o: isa_kernels_v3.cpp
	$(CXX) $(CXXFLAGS) $(ISA_OPTFLAGS) $(ISA_V3_FLAGS) $(ISA_DEFS) -fPIC -c $< -o $@

isa_kernels_v4.o: isa_kernels_v4.cpp
	$(CXX) $(CXXFLAGS) $(ISA_OPTFLAGS) $(ISA_V4_FLAGS) $(ISA_DEFS) -fPIC -c $< -o $@

unitylib: $(ISA_OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(ISA_DEFS) unity-build.cpp $(ISA_OBJS) -fPIC -shared -o libIsoSpec++.so

debug: $(ISA_OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) $(ISA_DEFS) unity-build.cpp $(ISA_OBJS) -DDEBUG -fPIC -shared -o libIsoSpec++.so

debug-gcc: $(ISA_OBJS)
	g++ $(CXXFLAGS) $(DEBUGFLAGS) $(ISA_DEFS) unity-build.cpp $(ISA_OBJS) -DDEBUG -fPIC -shared -o libIsoSpec++.so

nonunity: $(ISA_OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(ISA_DEFS) $(SRCFILES) isa_dispatch.cpp $(ISA_OBJS) -DDEBUG -fPIC -shared -o libIsoSpec++.so

# The library without runtime dispatch: one level, whatever $(OPTFLAGS) targets.
# This is what the R package and a plain `c++ unity-build.cpp` get.
singlelevel:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) unity-build.cpp -fPIC -shared -o libIsoSpec++.so

clean:
	rm -f libIsoSpec++.so $(ISA_OBJS)

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

.PHONY: all unitylib debug debug-gcc nonunity singlelevel clean windows check check11 asm
