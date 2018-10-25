#include <iostream>
#include "isoSpec++.h"
#include <cassert>

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
        size_t confs_no = i.count_confs();
        std::cout << "No. confs is: " << confs_no << std::endl;
        i.reset();
        IsoThresholdGenerator i2(argv[1], target_prob, true, 100, 100, true);
        IsoThresholdGenerator i3(argv[1], target_prob, true, 100, 100, false);
        int* confspace = new int[i.getAllDim()];
        int* confspace2 = new int[i2.getAllDim()];
        //int* confspace3 = new int[i3.getAllDim()]; // these will be in different order...
        int no_visited = 0;
        double total_prob = 0.0;
	while(i.advanceToNextConfiguration())
	{
                assert(i2.advanceToNextConfiguration());
                assert(i3.advanceToNextConfiguration());
                std::cout << "lprob: " << i.lprob() << " prob: " << i.eprob() << " log(prob): " << log(i.eprob()) << " mass: " << i.mass() << " conf: ";
                assert(i.lprob() == i2.lprob());
                assert(i.mass() == i2.mass());
                assert(i.eprob() == i2.eprob());
                i.get_conf_signature(confspace);
                i2.get_conf_signature(confspace2);
                assert(memcmp(confspace, confspace2, i.getAllDim()*sizeof(int)) == 0);
                printArray<int>(confspace, i.getAllDim());
		no_visited += 1;
                total_prob += i.eprob();
	}
        delete[] confspace;
        delete[] confspace2;
        assert(!i2.advanceToNextConfiguration());
        assert(!i3.advanceToNextConfiguration());
	std::cout << "The number of visited configurations is:" << no_visited << " for total prob of: " << total_prob << std::endl;

}
