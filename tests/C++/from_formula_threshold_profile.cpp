#include <iostream>
#include "isoSpec++.h"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_threshold_simple(const char* formula, double threshold);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_threshold C10000H1000O1000N1000 0.01" << std::endl;
		std::cout << "...will the configurations with probability above 0.01 for the above molecule" << std::endl;
		return -1;
	}
        size_t v = 0;
        for (int ii=0; ii<10000; ii++)
	    v += test_threshold_simple(argv[1], atof(argv[2]));
	
	std::cout << "The number of visited configurations is:" << v << std::endl;

}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_threshold_simple(const char* formula, double threshold)
{

	IsoThresholdGenerator i(formula, threshold, true, 100, 100, true);
        int no_visited = 0;
        double total_prob = 0.0;
	while(i.advanceToNextConfiguration())
	{
		no_visited += 1;
                total_prob += i.prob();
	}
	return no_visited;

}
