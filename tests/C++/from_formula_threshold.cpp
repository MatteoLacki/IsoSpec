#include <iostream>
#include "isoSpec++.h"

using namespace IsoSpec;

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_threshold C10000H1000O1000N1000 0.01" << std::endl;
		std::cout << "...will the configurations with probability above 0.01 for the above molecule" << std::endl;
		return -1;
	}
	double target_prob = atof(argv[2]);
	IsoThresholdGenerator i(argv[1], target_prob, true, 100, 100, true);
        int confspace[1000];
        int no_visited = 0;
        double total_prob = 0.0;
	while(i.advanceToNextConfiguration())
	{
                std::cout << "lprob: " << i.lprob() << " prob: " << i.eprob() << " log(prob): " << log(i.eprob()) << " mass: " << i.mass() << " conf: ";
                i.get_conf_signature(confspace);
                printArray<int>(confspace, i.getAllDim());
		no_visited += 1;
                total_prob += i.eprob();
	}
	std::cout << "The number of visited configurations is:" << no_visited << " for total prob of: " << total_prob << std::endl;

}
