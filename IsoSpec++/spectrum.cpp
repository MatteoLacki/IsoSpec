#include <math.h>
#include <string.h>





Spectrum::Spectrum(double _start, double _bucketsize, int _buckets, bool _clear)
: start(_start), bucketsize(_bucketsize), buckets(_buckets), end(_start + bucketsize * static_cast(float)(buckets))
{
	spectrum = new double[buckets];
	if(clear)
		memset(spectrum, buckets, sizeof(double));
}
	

Spectrum::FromIso(


