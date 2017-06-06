#include "isoSpec++.h"



class Spectrum
{
private:
        Iso* isoptr;
        Iso&& iso;
	double lowest_mass;
	const double bucket_width;
	unsigned long n_buckets;
	double* storage;
        double* ofset_store;
        pthread_t* threads;
        const double cutoff;
        PrecalculatedMarginal** PMs;
        unsigned int n_threads;
        bool absolute;
        std::atomic<unsigned int> thread_idxes;
        double** thread_storages;
        double* thread_partials;
        unsigned int* thread_numbers;
        unsigned int total_confs;
public:
	Spectrum(Iso&& I, double bucket_width, double cutoff, bool _absolute);
	~Spectrum();
	void add_other(Spectrum& other);
        void run(unsigned int threads = 0, bool sync = true);
        void worker_thread();
        void wait();
        void calc_sum();
	inline unsigned int get_total_confs() { return total_confs; };
	void print(std::ostream& o = std::cout);

};
