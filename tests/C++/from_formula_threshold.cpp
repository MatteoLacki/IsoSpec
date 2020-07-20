#include <iostream>
#include "isoSpec++.h"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_threshold(const char* formula, double threshold, bool print_confs);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_threshold C10000H1000O1000N1000 0.01" << std::endl;
		std::cout << "...will the configurations with probability above 0.01 for the above molecule" << std::endl;
                std::cout << "Proper usage (for example): ./from_formula_threshold C10000H1000O1000N1000 0.01 false" << std::endl;
                std::cout << "...will just count them" << std::endl;
		return -1;
	}

        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);

	size_t no_visited = test_threshold(argv[1], atof(argv[2]), print_confs);
	
	std::cout << "The number of visited configurations is:" << no_visited << std::endl;

}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_threshold(const char* formula, double threshold, bool print_confs)
{

	IsoThresholdGenerator i(formula, threshold, true, 100, 100, true);
        size_t confs_no = i.count_confs();
	//if(print_confs)
        //	std::cout << "No. confs is: " << confs_no << std::endl;
        i.reset();
        IsoThresholdGenerator i2(formula, threshold, true, 100, 100, true);
        IsoThresholdGenerator i3(formula, threshold, true, 100, 100, false);
        int* confspace = new int[i.getAllDim()];
        int* confspace2 = new int[i2.getAllDim()];
        //int* confspace3 = new int[i3.getAllDim()]; // these will be in different order...
        size_t no_visited = 0;
        double total_prob = 0.0;
	while(i.advanceToNextConfiguration())
	{
                assert(i2.advanceToNextConfiguration());
                assert(i3.advanceToNextConfiguration());
		if(print_confs)
                	std::cout << "lprob: " << i.lprob() << " prob: " << i.prob() << " log(prob): " << log(i.prob()) << " mass: " << i.mass() << " conf: ";
                assert(i.lprob() == i2.lprob());
                assert(i.mass() == i2.mass());
                assert(i.prob() == i2.prob());
                i.get_conf_signature(confspace);
                i2.get_conf_signature(confspace2);
                assert(memcmp(confspace, confspace2, i.getAllDim()*sizeof(int)) == 0);
		if(print_confs)
                	printArray<int>(confspace, i.getAllDim());
		no_visited += 1;
                total_prob += i.prob();
	}
        delete[] confspace;
        delete[] confspace2;
        assert(!i2.advanceToNextConfiguration());
        assert(!i3.advanceToNextConfiguration());
	assert(!i.advanceToNextConfiguration());
	assert(!i2.advanceToNextConfiguration());
	assert(!i3.advanceToNextConfiguration());
        assert(confs_no == no_visited);
	return no_visited;

}
