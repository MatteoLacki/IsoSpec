CXX=clang++
OPTFLAGS=-O2 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
CXXFLAGS=-std=c++11 
SRCFILES=cwrapper.cpp allocator.cpp  dirtyAllocator.cpp  isoSpec++.cpp   marginalTrek++.cpp  operators.cpp element_tables.cpp misc.cpp Rinterface.cpp RcppExports.cpp isoMath.cpp
CXX_STD = CXX11


RCPPFLAGS = `R CMD config CPPFLAGS` `R CMD config --cppflags`
RLDFLAGS =  `R CMD config LDFLAGS`  `R CMD config --ldflags`


RCPPINCL = `echo 'Rcpp:::CxxFlags()' | R --vanilla --slave`
RCPPLIBS = `echo 'Rcpp:::LdFlags()'  | R --vanilla --slave`

all: lib

lib: rcppstuff
	$(CXX) $(RCPPFLAGS) $(RLDFLAGS) $(RCPPINCL) $(RCPPLIBS) $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) -DBUILDING_R -fPIC -shared -o IsoSpecR.so 

debug: rcppstuff
	$(CXX) $(RCPPFLAGS) $(RLDFLAGS) $(RCPPINCL) $(RCPPLIBS) $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) -DBUILDING_R -fPIC -shared -o IsoSpecR.so

rcppstuff:
	echo "Rcpp:::compileAttributes('..')" | R --vanilla --slave

clean:
	rm -f libIsoSpec++.so
