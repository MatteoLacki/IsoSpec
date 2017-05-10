#include <iostream>
#include "isoSpec++.h"
#include "summator.h"


int main()
{

    IsoThresholdGenerator* iso = IsoSpec::IsoFromFormula("C169719H270464N45688O52237S911", 0.03);

    SSummator s;
    unsigned int cnt = 0;
    while(iso->advanceToNextConfiguration())
    {
    	cnt++;
	s.add(exp(iso->lprob()));
	if(cnt % 10000000 == 0)
		std::cout << cnt << "	" << s.get() << '\n';
    }

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << cnt << " element(s)" << std::endl;

    delete iso;
}
