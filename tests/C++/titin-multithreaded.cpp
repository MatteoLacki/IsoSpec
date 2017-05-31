#include <iostream>
#include "isoSpec++.h"
#include "summator.h"
#include "spectrum2.h"
#include <pthread.h>


#define n_threads 4

void* thread(void* nr);

double fin_probs[n_threads];
Spectrum* spectra[n_threads];
SyncMarginal* SM = nullptr;
const double threshold = 0.2;

int main()
{
    pthread_t threads[n_threads];
    int threadargs[n_threads];

    Iso I("C169719H270464N45688O52237S911");
    SM = I.get_last_marginal(1000, 1000, log(threshold)+I.getModeLProb());
    std::cout << SM << std::endl;
//    SM->reset();

    for(int  index = 0; index < n_threads; ++index )
    {
    	threadargs[index] = index;
        pthread_create( &threads[index], NULL, thread, (void*) &threadargs[index] );
    };
    for(int  index = 0; index < n_threads; ++index )
    	pthread_join( threads[ index ], NULL );

    double total = 0.0;
    for(int  index = 0; index < n_threads; ++index )
    	total += fin_probs[index];

    std::cout << "Final summary: prob: " << total << '\n';
/*
    for(int ii=1; ii<n_threads; ii++)
    {
    	spectra[0]->add_other(*spectra[ii]);
	delete spectra[ii];
    }
    spectra[0]->print();
    delete spectra[0];
*/}





void* thread(void* nr)
{
    int numer = *((int*) nr);

    char padding[64];
    SSummator s;
    unsigned int cnt_tot = 0;
        IsoThresholdGeneratorMT* iso = new IsoThresholdGeneratorMT("C169719H270464N45688O52237S911", threshold, SM, false);
        iso->advanceToNextConfiguration();
	std::cout << "Ready: " << exp(iso->lprob()) << "Range: " << iso->getLightestPeakMass() << " - " << iso->getHeaviestPeakMass() << std::endl;
        unsigned int cnt = 1;
        while (iso->advanceToNextConfiguration())
        {
    	    cnt++;
//	    s.add(exp(iso->lprob()));
//	    if(cnt % 10000000 == 0)
//		std::cout << cnt << "	" << s.get() << "\t" << exp(iso->lprob()) << '\n';
	//    last = exp(iso->lprob());
        };
        std::cout <<  "Slice: " << cnt << " element(s), last: " << "totalprob: " << s.get() << std::endl;
	/*
	spectra[numer] = new Spectrum (*iso, 1.0);
	Spectrum* spctr = spectra[numer];
	delete iso;
	cnt_tot += cnt;
	std::cout <<  "Slice: " << cnt << " element(s), last: " << "totalprob: " << spctr->sum.get() << std::endl;
	fin_probs[numer] = spctr->sum.get();*/
	return NULL;
}
