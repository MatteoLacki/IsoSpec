#include <cmath>
#include "spectrum2.h"
#include <assert.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>



Spectrum::Spectrum(IsoGenerator& G, double _bucket_width) : 
lowest_mass(G.getLightestPeakMass()),
bucket_width(_bucket_width),
n_buckets(static_cast<unsigned int>(ceil(G.getHeaviestPeakMass()-lowest_mass)/bucket_width))
{	
	double prob;
	long pagesize = sysconf(_SC_PAGESIZE);
	unsigned long mmap_len = n_buckets * sizeof(double);
	mmap_len += pagesize - mmap_len%pagesize;
	storage = reinterpret_cast<double*>(mmap(NULL, n_buckets*sizeof(double), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
	unsigned int diff = static_cast<unsigned int>(floor(lowest_mass/bucket_width));
	double* ofset_store = storage - diff;
	while(G.advanceToNextConfiguration())
	{
		prob = exp(G.lprob());
		ofset_store[static_cast<unsigned int>(floor(G.mass()/bucket_width))] += exp(G.lprob());
		sum.add(prob);
	}
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
	sum.add(other.sum.get());
	for(unsigned long ii=0; ii<n_buckets; ii++)
	    if(other.storage[ii] > 0.0) // This seemingly unnecesary if is here to avoid needless writes into mmaped memory
	        storage[ii] += other.storage[ii];
}

void Spectrum::print(std::ostream& o)
{
	for(unsigned long ii=0; ii<n_buckets; ii++)
	    o << lowest_mass + static_cast<double>(ii)*bucket_width << "\t" << storage[ii] << std::endl;
}
