#include <iostream>
#include "marginalTrek++.h"
#include "summator.h"
#include "spectrum2.h"


int main()
{

    SSummator s;
    unsigned int cnt_tot = 0;
    int total_t = 10;
    double threshold = 0.1;
    double masses[] = {1.0, 1000.0};
    double probs[] = {0.9, 0.1};
    Marginal m(masses, probs, 2, 100);
    LayeredMarginal mr(std::move(m));

    unsigned int cnt = 0;
    mr.extend(-10.0);
    

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << mr.get_no_confs() << " element(s)" << std::endl;
    std::cout << "prob: " << s.get() << std::endl;
#if 0
    IsoThresholdGeneratorBoundMass* isob = new IsoThresholdGeneratorBoundMass("C169719H270464N45688O52237S911", threshold, mmin, mmax, false);
    std::cout << isob->getModeLProb() << std::endl;

    unsigned int confsig[5];
    double cnt = 1.0;
    cnt_tot = 0;
    double lc = isob->getModeLProb() + log(threshold);
    while(isob->advanceToNextConfiguration())
    {
    	cnt_tot++;
    }
    delete isob;

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << cnt_tot << " element(s)" << std::endl;
#endif
}
