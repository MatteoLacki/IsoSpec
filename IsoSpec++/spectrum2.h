#include "isoSpec++.h"



class Spectrum
{
private:
	double lowest_mass;
	double bucket_width;
	unsigned long n_buckets;
	double* storage;
public:
	Spectrum(IsoGenerator& G, double bucket_width);
	~Spectrum();
	Summator sum;
	void add_other(Spectrum& other);
	
	void print(std::ostream& o = std::cout);

};
