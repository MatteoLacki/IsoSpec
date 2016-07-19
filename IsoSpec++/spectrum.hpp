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


class FunctionalKernel
{
public:
	virtual double getMass(double bucketStart, double bucketEnd) = 0;
	virtual double getSupportMin() = 0;
	virtual double getSupportMax() = 0;
};

class SinglePointFunctionalKernel : FunctionalKernel
{
public:
	SinglePointFunctionalKernel();
	virtual double getMass(double bucketStart, double bucketEnd);
	virtual double getSupportMin();
	virtual double getSupportMax();
};

class TruncatedGaussianFunctionalKernel : FunctionalKernel
{
	double stdev;
	double prob;
	double support_min;
	double support_max;
//	double support_len;
	double correction;
public:
	TruncatedGaussianFunctionalKernel(double _stdev, double _prob = 0.99);
	double getMass(double bucketStart, double bucketEnd);
	double getSupportMin();
	double getSupportMax();
};


class RectangularFunctionalKernel : FunctionalKernel
{
	double support_min;
	double support_max;
	double support_len;
public:
	RectangularFunctionalKernel(double start, double end);
	double getSupportMin();
	double getSupportMax();
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

	inline unsigned int position(double mass)
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

	inline double mass_at_index_start(unsigned int idx)
	{
		return start + bucketsize * idx;
	}

	Spectrum(double _start = -0.5, double _bucketsize = 1.0, int _buckets = 1, bool _clear = true);
	Spectrum(IsoSpec& iso, FunctionalKernel& k, double _bucketsize);

	

};


#endif /* SPECTRUM_HPP */
