#include <iostream>
#include "isoSpec++.h"
#include "tabulator.h"

using namespace IsoSpec;

size_t test_layered_tabulator(const char* formula, double total_prob, bool print_confs = false);

#ifndef ISOSPEC_TESTS_SKIP_MAIN
int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_layered C10000H1000O1000N1000 0.9999" << std::endl;
		std::cout << "...will print the minimal number of configurations necessary to cover 0.9999 probability of the above molecule" << std::endl;
		return -1;
	}

	size_t no_confs = test_layered_tabulator(argv[1], atof(argv[2]), true);

	std::cout << "The number of visited configurations is:" << no_confs << std::endl;
}
#endif /* #ifndef ISOSPEC_TESTS_SKIP_MAIN */

size_t test_layered_tabulator(const char* formula, double total_prob, bool print_confs)
{
	IsoLayeredGenerator i(formula, 1000, 1000);
        LayeredTabulator t(&i, true, true, true, true, total_prob, true);
        double* probs = t.probs(false);
        double* masses = t.masses(true);
        int* confs = t.confs();

        for(size_t ii = 0; ii<t.confs_no(); ii++)
	{
		if(print_confs)
		{
			std::cout << "PROB: " << probs[ii] << "  \tMASS: " << masses[ii] << "\tCONF: ";
                        int* space = confs + ii*i.getAllDim();
			for(int ii=0; ii<i.getAllDim(); ii++)
			    std::cout << space[ii] << " ";
			std::cout << std::endl;
		}

	}

        free(masses);

	return t.confs_no();
}
