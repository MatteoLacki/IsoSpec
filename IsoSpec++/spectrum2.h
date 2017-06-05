#include "isoSpec++.h"



class Spectrum
{
private:
        Iso&& iso;
	double lowest_mass;
	const double bucket_width;
	unsigned long n_buckets;
	double* storage;
        double* ofset_store;
        pthread_t* threads;
        const double cutoff;
        SyncMarginal* SM;
        unsigned int n_threads;
        bool absolute;
        std::atomic<unsigned int> thread_idxes;
        double** thread_storages;
        double* thread_partials;
public:
	Spectrum(Iso&& I, double bucket_width, double cutoff, bool _absolute);
	~Spectrum();
	void add_other(Spectrum& other);
        void run(unsigned int threads = 0, bool sync = true);
        void worker_thread();
        void wait();
        void calc_sum();
	
	void print(std::ostream& o = std::cout);

};
