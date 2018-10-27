#include <iostream>
#include "isoSpec++.h"

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN
size_t test_ordered(const char* formula, double total_prob, bool print_confs = false);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): " << argv[0] << " C10000H1000O1000N1000 0.9999" << std::endl;
		std::cout << "...will print the minimal number of configurations necessary to cover 0.9999 probability of the above molecule" << std::endl;
		return -1;
	}
	
	#ifndef ISOSPEC_TESTS_MEMSAN
	size_t no_confs = 
	#endif 
			  test_ordered(argv[1], atof(argv[2]), true);

	#ifndef ISOSPEC_TESTS_MEMSAN
	std::cout << "The number of visited configurations is: " << no_confs << std::endl;
	#endif
}

#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_ordered(const char* formula, double total_prob, bool print_confs)
{
	IsoOrderedGenerator i(formula);
	double target_prob = total_prob;

	size_t no_visited = 0;
        int* space = new int[i.getAllDim()];
	while(target_prob > 0.0 && i.advanceToNextConfiguration())
	{
		target_prob -= i.eprob();
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
