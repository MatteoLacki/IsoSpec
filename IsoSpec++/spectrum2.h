#include "isoSpec++.h"



class Spectrum
{
private:
	double lowest_mass;
	double bucket_width;
	unsigned int n_buckets;
	double* storage;
public:
	Spectrum(IsoGenerator& G, double bucket_width);

};
