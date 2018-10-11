#include <iostream>
#include "isoSpec++.h"

using namespace IsoSpec;

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula C10000H1000O1000N1000 0.9999" << std::endl;
		std::cout << "...will print the minimal number of configurations necessary to cover 0.9999 probability of the above molecule" << std::endl;
		return -1;
	}
	IsoOrderedGenerator i(argv[1]);
	double target_prob = atof(argv[2]);
	int no_visited = 0;
	while(target_prob > 0.0 && i.advanceToNextConfiguration())
	{
		target_prob -= i.eprob();
		no_visited += 1;
	}
	std::cout << "The number of visited configurations is:" << no_visited << std::endl;

}
