#include <iostream>
#include "isoSpec++.h"
#include "summator.h"
#include "spectrum2.h"
#include <pthread.h>


#define n_threads 64

void* thread(void* nr);

double fin_probs[n_threads];

int main()
{
    pthread_t threads[n_threads];
    int threadargs[n_threads];

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

}





void* thread(void* nr)
{
    int numer = *((int*) nr);


    SSummator s;
    unsigned int cnt_tot = 0;
    double threshold = 0.0001;
        IsoThresholdGeneratorMultithreaded* iso = new IsoThresholdGeneratorMultithreaded(n_threads, numer, "C169719H270464N45688O52237S911", threshold, false);
	std::cout << "Ready: " << exp(iso->lprob()) << "Range: " << iso->getLightestPeakMass() << " - " << iso->getHeaviestPeakMass() << std::endl;
        unsigned int cnt = 0;
/*        while (iso->advanceToNextConfiguration())
        {
    	    cnt++;
	    s.add(exp(iso->lprob()));
	    if(cnt % 10000000 == 0)
		std::cout << cnt << "	" << s.get() << "\t" << exp(iso->lprob()) << '\n';
	    last = exp(iso->lprob());
        };
*/	
	Spectrum spctr(*iso, 0.1);
	delete iso;
	cnt_tot += cnt;
	std::cout <<  "Slice: " << cnt << " element(s), last: " << "totalprob: " << spctr.sum.get() << std::endl;
	fin_probs[numer] = spctr.sum.get();
	return NULL;
}
