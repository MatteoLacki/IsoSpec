#include <cmath>
#include "spectrum2.h"




Spectrum::Spectrum(IsoGenerator& G, double _bucket_width) : 
lowest_mass(G.getLightestPeakMass()),
bucket_width(_bucket_width),
n_buckets(static_cast<unsigned int>(ceil(G.getHeaviestPeakMass()-lowest_mass)/bucket_width)),
storage(new double[n_buckets])
{
	memset(storage, 0, n_buckets*sizeof(double));
	unsigned int diff = static_cast<unsigned int>(floor(lowest_mass/bucket_width));
	double* ofset_store = storage - diff;
	std::cerr << diff << "n_buck: " << n_buckets << '\n';
	while(G.advanceToNextConfiguration())
	{
		std::cerr << lowest_mass << " " << G.mass() << " " << G.getHeaviestPeakMass() << '\n';
		std::cerr << static_cast<unsigned int>(floor(G.mass()/bucket_width));
		ofset_store[static_cast<unsigned int>(floor(G.mass()/bucket_width))] += exp(G.lprob());
	}
}
