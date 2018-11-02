#include <iostream>
#include "isoSpec++.h"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_zero(double, bool);

int main(int argc, char** argv)
{
	if(argc < 1) // yeah... never.
	{
		std::cout << "Proper usage (for example): " << argv[1] << std::endl;
		std::cout << "...will hopefully throw an exception instead of segfaulting." << std::endl;
		return -1;
	}
	size_t nc = test_zero(0.1, true);
        std::cout << "No confs visited: " << nc << std::endl;
	

}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_zero(double threshold, bool print_confs)
{
        int dimNumber = 1;
        int isoNumbers[] = { 3 };
        int atomCounts[] = { 100 };
        double isoMasses[] = { 1.0, 2.0, 3.0 };
        double isoProbs[] = { 0.0, 0.4, 0.6 };

        double* IM = { isoMasses };
        double* IP = { isoProbs };

        Iso ii(dimNumber, isoNumbers, atomCounts, &IM, &IP);
	IsoThresholdGenerator i(std::move(ii), false, threshold);
        size_t confs_no = i.count_confs();
	if(print_confs)
        	std::cout << "No. confs is: " << confs_no << std::endl;
        i.reset();
        int* confspace = new int[i.getAllDim()];
        int no_visited = 0;
        double total_prob = 0.0;
	while(i.advanceToNextConfiguration())
	{
                i.get_conf_signature(confspace);
		if(print_confs)
                	std::cout << "lprob: " << i.lprob() << " prob: " << i.prob() << " log(prob): " << log(i.prob()) << " mass: " << i.mass() << " conf: ";
		if(print_confs)
                	printArray<int>(confspace, i.getAllDim());
		no_visited += 1;
                total_prob += i.prob();
	}
        delete[] confspace;
	return no_visited;

}
