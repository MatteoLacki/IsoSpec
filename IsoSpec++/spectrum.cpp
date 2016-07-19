#include <math.h>
#include <string.h>
#include <iostream>
#include "spectrum.hpp"
#include "isoMath.hpp"
#include "isoSpec++.hpp"

Kernel::Kernel(double _delta, double* _k, double _bucketsize, double _buckets) : 
delta(_delta),
k(_k),
bucketsize(_bucketsize),
buckets(_buckets)
{}

Kernel* Kernel::SinglePoint(double bucketsize)
{
    double* k = new double[1];
    k[0] = 1.0;
    return new Kernel(bucketsize/2.0, k, bucketsize, 1);
}

Kernel* Kernel::Gaussian(double stdev, double bucketsize, double prob)
{
	double rg_end = -NormalCDFInverse((1.0 - prob)/2.0, 0.0, stdev);
	double bucklen = ceil(rg_end - bucketsize/2.0);
	unsigned int buck_offset = static_cast<unsigned int>(bucklen);
	unsigned int buckets = 2 * buck_offset + 1;
	double* k = new double[buckets];
	for (int ii=0; ii<buckets; ii++)
		k[ii] = NormalPDF((ii - buck_offset) * bucketsize, 0.0, stdev);
	return new Kernel((static_cast<double>(buck_offset) + 0.5) * bucketsize, k, bucketsize, buckets);
}

void Kernel::print()
{
    for( unsigned int ii = 0; ii < buckets; ii++ )
    {
    	std::cout << k[ii] << std::endl;
    }
}


Spectrum::Spectrum(double _start, double _bucketsize, int _buckets, bool _clear)
: start(_start), end(_start + _bucketsize * static_cast<float>(_buckets)), bucketsize(_bucketsize), buckets(_buckets)
{
	spectrum = new double[buckets];
	memset(spectrum, buckets, sizeof(double));
}
	

Spectrum::Spectrum(IsoSpec& iso, Kernel& k)
{
	bucketsize = k.bucketsize;

	double lpm = iso.getLightestPeakMass();
	double hpm = iso.getHeaviestPeakMass();

	start = floor(lpm) - bucketsize * 1.5;
	buckets = floor((hpm - start)/bucketsize) + 2;
	end = start + buckets * bucketsize;

	spectrum = new double[buckets];
	memset(spectrum, buckets, sizeof(double));
}



