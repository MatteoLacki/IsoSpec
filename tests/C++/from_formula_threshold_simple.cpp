#include <iostream>
#include "isoSpec++.h"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_threshold_simple(const char* formula, double threshold, bool print_confs);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_threshold C10000H1000O1000N1000 0.01" << std::endl;
		std::cout << "...will the configurations with probability above 0.01 for the above molecule" << std::endl;
		return -1;
	}

        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);

	size_t no_visited = test_threshold_simple(argv[1], atof(argv[2]), print_confs);
	
	std::cout << "The number of visited configurations is:" << no_visited << std::endl;

}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_threshold_simple(const char* formula, double threshold, bool print_confs)
{

	IsoThresholdGenerator i(formula, threshold, true, 10, 10, true);
#ifdef TEST_SIZE
        size_t confs_no = i.count_confs();
	if(print_confs)
        	std::cout << "No. confs is: " << confs_no << std::endl;
        i.reset();
#endif
        int* confspace = print_confs ? new int[i.getAllDim()] : nullptr;
        size_t no_visited = 0;
        double total_prob = 0.0;
	while(i.advanceToNextConfiguration())
	{
            if(print_confs)
            {
                i.get_conf_signature(confspace);
                std::cout << "lprob: " << i.lprob() << " prob: " << i.prob() << " log(prob): " << log(i.prob()) << " mass: " << i.mass() << " conf: ";
                printArray<int>(confspace, i.getAllDim());
            }
            no_visited += 1;
            total_prob += i.prob();
	}
	delete[] confspace;
	return no_visited;
}
