#include <iostream>
#include <cassert>
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
                std::cout << argv[0] << " C10000H1000O1000N1000 0.9999 false" << std::endl;
                std::cout << "will just count them" << std::endl;
		return -1;
	}
        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);
	
	#ifndef ISOSPEC_TESTS_MEMSAN
	size_t no_confs = 
	#endif 
			  test_ordered(argv[1], atof(argv[2]), print_confs);

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
        double last_p = 1.0;
	while(target_prob > 0.0 && i.advanceToNextConfiguration())
	{
                double curr_p = i.prob();
                assert(last_p >= curr_p);
                last_p = curr_p;
		target_prob -= curr_p;
		no_visited += 1;
		if(print_confs)
		{
			std::cout << "PROB: " << i.prob() << "  \tMASS: " << i.mass() << "\tCONF: ";
			i.get_conf_signature(space);
			for(int ii=0; ii<i.getAllDim(); ii++)
			    std::cout << space[ii] << " ";
			std::cout << std::endl;
		}

	}
        delete[] space;
	return no_visited;

}
