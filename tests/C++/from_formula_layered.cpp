#include <iostream>
#include "isoSpec++.h"

using namespace IsoSpec;

size_t test_layered(const char* formula, double total_prob, bool print_confs = false);

#ifndef ISOSPEC_TESTS_SKIP_MAIN
int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_layered C10000H1000O1000N1000 0.9999" << std::endl;
		std::cout << "...will print the minimal number of configurations necessary to cover 0.9999 probability of the above molecule" << std::endl;
		return -1;
	}

	size_t no_confs = test_layered(argv[1], atof(argv[2]), true);

	std::cout << "The number of visited configurations is:" << no_confs << std::endl;
}
#endif /* #ifndef ISOSPEC_TESTS_SKIP_MAIN */

size_t test_layered(const char* formula, double total_prob, bool print_confs)
{
	IsoLayeredGenerator i(formula, total_prob, 0.3, 1000, 1000, true);
	size_t no_visited = 0;
        int* space = new int[i.getAllDim()];
	while(i.advanceToNextConfiguration())
	{
		no_visited += 1;
		if(print_confs)
		{
			std::cout << "EPROB: " << i.eprob() << "  \tMASS: " << i.mass() << "\tCONF: ";
			i.get_conf_signature(space);
			for(int ii=0; ii<i.getAllDim(); ii++)
			    std::cout << space[ii] << " ";
			std::cout << std::endl;
		}

	}
        delete[] space;

	return no_visited;
}
