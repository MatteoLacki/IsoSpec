#include <cmath>
#include "spectrum2.h"
#include <assert.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/sysinfo.h>




Spectrum::Spectrum(Iso&& I, double _bucket_width, double _cutoff, bool _absolute) : 
iso(std::move(I)),
lowest_mass(I.getLightestPeakMass()),
bucket_width(_bucket_width),
n_buckets(static_cast<unsigned int>(ceil(I.getHeaviestPeakMass()-lowest_mass)/bucket_width)),
cutoff(_cutoff),
absolute(_absolute),
thread_idxes(0)
{
        PMs = I.get_MT_marginal_set(log(cutoff), absolute, 1024, 1024);
	long pagesize = sysconf(_SC_PAGESIZE);
	unsigned long mmap_len = n_buckets * sizeof(double);
	mmap_len += pagesize - mmap_len%pagesize;
	storage = reinterpret_cast<double*>(mmap(NULL, n_buckets*sizeof(double), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
	unsigned int diff = static_cast<unsigned int>(floor(lowest_mass/bucket_width));
	ofset_store = storage - diff;
}

void* wrapper_func_thr(void* spc)
{
    reinterpret_cast<Spectrum*>(spc)->worker_thread();
    return NULL;
}

void Spectrum::run(unsigned int nthreads, bool sync)
{
    if(nthreads == 0)
        nthreads = get_nprocs();
    n_threads = nthreads;

    threads = new pthread_t[n_threads];
    thread_storages = new double*[n_threads];
    thread_partials = new double[n_threads];
    thread_numbers = new unsigned int[n_threads];

    for(unsigned int ii = 0; ii < n_threads; ii++)
        pthread_create(&threads[ii], NULL, wrapper_func_thr, this);

    if(sync)
        wait();
}

void Spectrum::wait()
{
    for(unsigned int ii = 0; ii < n_threads; ii++)
        pthread_join(threads[ii], NULL);

    delete[] threads;

    calc_sum();
}

void Spectrum::calc_sum()
{
    total_confs = 0;
    total_prob = 0;

    for(unsigned int ii = 0; ii < n_threads; ii++)
    {
        total_confs += thread_numbers[ii];
        total_prob += thread_partials[ii];
    };
}


void Spectrum::worker_thread()
{
    unsigned int thread_id = thread_idxes.fetch_add(1);
    IsoThresholdGeneratorMT* isoMT = new IsoThresholdGeneratorMT(std::move(iso), cutoff, PMs, absolute);
    double* local_storage = reinterpret_cast<double*>(mmap(NULL, n_buckets*sizeof(double), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
    unsigned int diff = static_cast<unsigned int>(floor(lowest_mass/bucket_width));
    double* local_ofset_store = local_storage - diff;
    double prob;
    Summator sum;
    unsigned int cnt = 0;
    while(isoMT->advanceToNextConfiguration())
    {
        prob = isoMT->eprob();
        local_ofset_store[static_cast<unsigned int>(floor(isoMT->mass()/bucket_width))] += prob;
        sum.add(prob);
        cnt++;
    }
    thread_storages[thread_id] = local_storage;
    thread_partials[thread_id] = sum.get();
    thread_numbers[thread_id] = cnt;
    delete isoMT;
}

Spectrum::~Spectrum()
{
	munmap(storage, n_buckets*sizeof(double));
}

void Spectrum::add_other(Spectrum& other)
{
	assert(n_buckets == other.n_buckets);
	assert(bucket_width == other.bucket_width);
	assert(lowest_mass == other.lowest_mass);
	for(unsigned long ii=0; ii<n_buckets; ii++)
	    if(other.storage[ii] > 0.0) // This seemingly unnecesary if is here to avoid needless writes into mmaped memory
	        storage[ii] += other.storage[ii];
}

void Spectrum::print(std::ostream& o)
{
	for(unsigned long ii=0; ii<n_buckets; ii++)
	    o << lowest_mass + static_cast<double>(ii)*bucket_width << "\t" << storage[ii] << std::endl;
}
