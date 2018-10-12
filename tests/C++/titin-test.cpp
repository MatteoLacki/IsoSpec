#include <iostream>
#include "isoSpec++.h"
#include "summator.h"

using namespace IsoSpec;

int main()
{

    SSummator s;
//    unsigned int cnt_tot = 0;
//    int total_t = 10;
    double threshold = 0.1;
//    double mmin = 3815900.0;
//    double mmax = 3816000.0;
    double mmin = -100000000000.0;
    double mmax = 100000000000.0;
    #if 1
    IsoThresholdGenerator* iso = new IsoThresholdGenerator("C169719H270464N45688O52237S911", threshold, false);
    unsigned int cnt = 0;
    while(iso->advanceToNextConfiguration())
    {
        if(iso->mass() >= mmin and mmax >= iso->mass())
            cnt++;
        s.add(exp(iso->lprob()));
    }
    delete iso;

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << cnt << " element(s)" << std::endl;
    std::cout << "prob: " << s.get() << std::endl;
    #endif
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
