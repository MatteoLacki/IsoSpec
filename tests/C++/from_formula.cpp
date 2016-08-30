#include <iostream>
#include "isoSpec++.h"

int main(int argc, char** argv)
{
	auto i = IsoSpec::IsoFromFormula(argv[1], atof(argv[2]));
	i->processConfigurationsUntilCutoff();
	std::cout << i->getNoVisitedConfs() << std::endl;
	delete i;

}
