#include <iostream>
#include "isoSpec++.h"
#include "summator.h"
#include "spectrum2.h"


int main()
{

    SSummator s;
    unsigned int cnt_tot = 0;
    int total_t = 10;
    double threshold = 0.5;
/*    for (int ii=0; ii<total_t; ii++)
    {
        IsoThresholdGeneratorMultithreaded* iso = new IsoThresholdGeneratorMultithreaded(total_t, ii, "C169719H270464N45688O52237S911", threshold, false);
	std::cout << exp(iso->lprob()) << std::endl;
	double last = 0.0;
        unsigned int cnt = 0;
        do
        {
    	    cnt++;
	    s.add(exp(iso->lprob()));
	    if(cnt % 10000000 == 0)
		std::cout << cnt << "	" << s.get() << "\t" << exp(iso->lprob()) << '\n';
	    last = exp(iso->lprob());
        } while (iso->advanceToNextConfiguration());
	delete iso;
	cnt_tot += cnt;
	std::cout <<  "Slice: " << cnt << " element(s), last: " << last << std::endl;
    }

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << cnt_tot << " element(s)" << std::endl;
*/
    IsoThresholdGenerator* iso = new IsoThresholdGenerator("C169719H270464N45688O52237S911", threshold, false);

    unsigned int confsig[5];
    double cnt = 1.0;
    cnt_tot = 0;
    while(iso->advanceToNextConfiguration())
    {
    //    cnt *= 1.000000001;
    	cnt_tot++;
//        s.add(exp(iso->lprob()));
//        std::cout << "Lprob: " << iso->lprob() << std::endl;
//        iso->get_conf_signature(confsig);
//        printArray(confsig, 5);
    }
    delete iso;

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << cnt_tot << " element(s)" << std::endl;
    std::cout << s.get() << cnt << std::endl;

    //IsoThresholdGenerator* isosp = new IsoThresholdGenerator("C169719H270464N45688O52237S911", threshold, false);
    //Spectrum sp(*isosp, 0.01);
    //delete isosp;
}
