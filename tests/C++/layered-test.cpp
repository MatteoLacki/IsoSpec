#include <iostream>
#include "isoSpec++.h"
#include "summator.h"

using namespace IsoSpec;


int main()
{

    SSummator s;
//    unsigned int cnt_tot = 0;
//    int total_t = 10;
    double threshold = 0.999;
//    double mmin = 3815900.0;
//    double mmax = 3816000.0;
    double mmin = -100000000000.0;
    double mmax = 100000000000.0;

    IsoLayeredGenerator* iso = new IsoLayeredGenerator("C1600H2700N1000", threshold, 0.3);
//    IsoLayeredGenerator* iso = new IsoLayeredGenerator("H2C2", 0.5, 0.3);
    unsigned int cnt = 0;
    int cntr[10];
    while(iso->advanceToNextConfiguration())
    {
        if(iso->mass() >= mmin and mmax >= iso->mass())
            cnt++;
        iso->get_conf_signature(cntr);

//        std::cout << cntr[0] << " " << cntr[1] << " " << cntr[2] << " " << " " << cntr[3] << " " << cntr[4] << " " << cntr[5] << " " << iso->lprob() << std::endl;
        if(cnt%10000 == 0)
            std::cout << cnt << " " << s.get() << std::endl;

        s.add(exp(iso->lprob()));
    }
    delete iso;

    std::cout <<  "The isotopologue set containing at least " << threshold << " probability has " << cnt << " element(s)" << std::endl;
    std::cout << "prob: " << s.get() << std::endl;
}
