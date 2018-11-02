#include <iostream>
#include "isoSpec++.h"
#include "summator.h"

using namespace IsoSpec;

int main()
{

    SSummator s;
//    unsigned int cnt_tot = 0;
//    int total_t = 10;
    double threshold = 0.01;
//    double mmin = 3815900.0;
//    double mmax = 3816000.0;
    double mmin = -100000000000.0;
    double mmax = 100000000000.0;
    const char formula[] = "C169719H270464N45688O52237S911";
//    const char formula[] = "H2O1";
//    const char formula[] = "C63H98N18O13S1"; // substance P
//    const char formula[] = "C520H817N139O147S8"; // Human insulin
    IsoThresholdGenerator* iso = new IsoThresholdGenerator(formula, threshold, false);
    unsigned int cnt = 0;
    while(iso->advanceToNextConfiguration())
    {
        if(iso->mass() >= mmin and mmax >= iso->mass())
            cnt++;
        s.add(iso->prob());
    }
    delete iso;

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << cnt << " element(s)" << std::endl;
    std::cout << "prob: " << s.get() << std::endl;
}
