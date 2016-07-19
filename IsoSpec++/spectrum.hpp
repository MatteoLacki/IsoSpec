#ifndef SPECTRUM_HPP
#define SPECTRUM_HPP
#include <math.h>
#include "isoSpec++.hpp"


class Kernel
{
public:
	const double delta;
	const double* k;
	const double bucketsize;
	const int buckets;

	Kernel(double _width, double* _k, double _bucketsize, double _buckets);

	static Kernel* SinglePoint(double _bucketsize);
	static Kernel* Gaussian(double stdev, double bucketsize, double prob);
	static Kernel* Rectangular(int width, double bucketsize);
	static Kernel* Triangular(int width, double bucketsize);

	void print();
};



class Spectrum
{
public:
//	double* kernel = nullptr;
	double* spectrum = nullptr;
	double  start = 0.0;
	double  end = 1.0;
	double  bucketsize = 1.0;
	int	buckets = 0;

	inline int position(double mass)
	{
		if(mass >= end)
			return buckets - 1;
		if(mass <= start)
			return 0;
		return static_cast<int>((mass - start) / bucketsize);
	};

	inline double value(double mass)
	{
		return spectrum[position(mass)];
	};

	Spectrum(double _start = -0.5, double _bucketsize = 1.0, int _buckets = 1, bool _clear = true);
	Spectrum(IsoSpec& iso, Kernel& k);

	

};


#endif /* SPECTRUM_HPP */
