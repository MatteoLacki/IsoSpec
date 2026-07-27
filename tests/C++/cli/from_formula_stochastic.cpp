#include <iostream>
#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_stochastic(const char* formula, size_t molecules, double precision, bool print_confs);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_stochastic C10000H1000O1000N1000 10000000" << std::endl;
		std::cout << "...will the configurations with ionic current of 10000000 for the above molecule" << std::endl;
                std::cout << "Proper usage (for example): ./from_formula_stochastic C10000H1000O1000N1000 10000000 false" << std::endl;
                std::cout << "...will just count them" << std::endl;
		return -1;
	}

        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);

	size_t no_visited = test_stochastic(argv[1], atoi(argv[2]), 0.999, print_confs);
	
	std::cout << "The number of visited configurations is:" << no_visited << std::endl;

}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_stochastic(const char* formula, size_t molecules, double precision, bool print_confs)
{
    FixedEnvelope fe = FixedEnvelope::FromStochastic(formula, molecules, precision, 5.0, print_confs);
    if(!print_confs)
        return fe.confs_no();

    const double* probs = fe.probs();
    double* masses = fe.release_masses();
    const int* confs = fe.confs();

    for(size_t ii=0; ii<fe.confs_no(); ii++)
    {
        std::cout << "PROB: " << probs[ii] << "  \tMASS: " << masses[ii] << "\tCONF: ";
        const int* space = confs + ii*fe.getAllDim();
        for(int ii=0; ii<fe.getAllDim(); ii++)
            std::cout << space[ii] << " ";
        std::cout << std::endl;

    }
    free(masses);
    return 0;
}
