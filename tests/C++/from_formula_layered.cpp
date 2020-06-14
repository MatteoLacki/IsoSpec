#include <iostream>
#include "isoSpec++.h"
#include "fixedEnvelopes.h"

using namespace IsoSpec;

size_t test_layered_tabulator(const char* formula, double total_prob, bool print_confs = false);

#ifndef ISOSPEC_TESTS_SKIP_MAIN
int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_layered C10000H1000O1000N1000 0.9999" << std::endl;
		std::cout << "...will print the configurations necessary to cover 0.9999 probability of the above molecule" << std::endl;
                std::cout << "./from_formula_layered C10000H1000O1000N1000 0.9999 false" << std::endl;
                std::cout << "will just count them" << std::endl;

		return -1;
	}

        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);

	size_t no_confs = test_layered_tabulator(argv[1], atof(argv[2]), print_confs);

	std::cout << "The number of visited configurations is:" << no_confs << std::endl;
}
#endif /* #ifndef ISOSPEC_TESTS_SKIP_MAIN */

size_t test_layered_tabulator(const char* formula, double total_prob, bool print_confs)
{
//	IsoLayeredGenerator i(formula, 1000, 1000);
        FixedEnvelope t = FixedEnvelope::FromTotalProb(formula, total_prob, true, print_confs);
        const double* probs = t.probs();
        double* masses = t.release_masses();
        const int* confs = t.confs();

        if(print_confs)
            for(size_t ii = 0; ii<t.confs_no(); ii++)
            {
                std::cout << "PROB: " << probs[ii] << "  \tMASS: " << masses[ii] << "\tCONF: ";
                const int* space = confs + ii*t.getAllDim();
                for(int ii=0; ii<t.getAllDim(); ii++)
                    std::cout << space[ii] << " ";
                std::cout << std::endl;

	    }

        free(masses);

	return t.confs_no();
}
