#include <math.h>
#include <string.h>


Kernel::Kernel(double _delta, double* _k, double _bucketsize, double buckets) : 
delta(_delta),
k(_k),
bucketsize(_bucketsize),
buckets(_buckets)
{}

Kernel::SinglePoint(double _bucketsize)
{
    double* k = new double[1];
    k[0] = 1.0;
    return new Kernel(_bucketsize/2.0, k, _bucketsize, 1);
}

Kernel::Gaussian


Spectrum::Spectrum(double _start, double _bucketsize, int _buckets, bool _clear)
: start(_start), bucketsize(_bucketsize), buckets(_buckets), end(_start + bucketsize * static_cast(float)(buckets))
{
	spectrum = new double[buckets];
	memset(spectrum, buckets, sizeof(double));
}
	

Spectrum::Spectrum(IsoSpec& iso, Kernel& k, double _bucketsize)
{
	bucketsize = _bucketsize;

	double lpm = iso.getLightestPeakMass();
	double hpm = iso.getHeaviestPeakMass();

	start = floor(lpm) - bucketsize * 1.5;
	buckets = floor((hpm - start)/bucketsize) + 2;
	end = start + buckets * bucketsize;

	spectrum = new double[buckets];
	memset(spectrum, buckets, sizeof(double));
}



