CXX=clang++
OPTFLAGS=-O3 -march=native -mtune=native
DEBUGFLAGS=-O0 -g
CXXFLAGS=-std=c++11 
SRCFILES=cwrapper.cpp allocator.cpp  dirtyAllocator.cpp  isoSpec++.cpp   marginalTrek++.cpp  operators.cpp element_tables.cpp misc.cpp Rinterface.cpp RcppExports.cpp isoMath.cpp
CXX_STD = CXX11

R_HOME = $(shell R RHOME)

#RCPPFLAGS = $(shell $(R_HOME)/bin/R CMD config --cppflags)
#RLDFLAGS = $(shell $(R_HOME)/bin/R CMD config --ldflags)
RCPPFLAGS = $(shell $(R_HOME)/bin/R CMD config CPPFLAGS) $(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS = $(shell $(R_HOME)/bin/R CMD config LDFLAGS) $(shell $(R_HOME)/bin/R CMD config --ldflags) 


RCPPINCL = $(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS = $(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

all: lib

lib: rcppstuff
	$(CXX) $(RCPPFLAGS) $(RLDFLAGS) $(RCPPINCL) $(RCPPLIBS) $(CXXFLAGS) $(OPTFLAGS) $(SRCFILES) -DBUILDING_R -fPIC -shared -o IsoSpecR.so 

debug: rcppstuff
	$(CXX) $(RCPPFLAGS) $(RLDFLAGS) $(RCPPINCL) $(RCPPLIBS) $(CXXFLAGS) $(DEBUGFLAGS) $(SRCFILES) -DBUILDING_R -fPIC -shared -o IsoSpecR.so

rcppstuff:
	echo "Rcpp:::compileAttributes('..')" | $(R_HOME)/bin/R --vanilla --slave
clean:
	rm -f libIsoSpec++.so
