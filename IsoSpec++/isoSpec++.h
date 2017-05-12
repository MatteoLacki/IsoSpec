/*
 *   Copyright (C) 2015-2016 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */


#ifndef ISOSPEC_PLUS_PLUS_HPP
#define ISOSPEC_PLUS_PLUS_HPP

#include <tuple>
#include <unordered_map>
#include <queue>
#include "lang.h"
#include "dirtyAllocator.h"
#include "summator.h"
#include "operators.h"
#include "marginalTrek++.h"


#ifdef BUILDING_R
 #include <Rcpp.h>
 using namespace Rcpp;
#endif /* BUILDING_R */


class IsoSpecLayered;
class IsoThresholdGenerator;

class Iso {
private:
	void setupMarginals(const double** _isotopeMasses, const double** _isotopeProbabilities);
protected:
	int 			dimNumber;
	int*			isotopeNumbers;
	int*			atomCounts;
	unsigned int		confSize;
	int			allDim;
	MarginalTrek**          marginalResults;
	const int             	tabSize;
        const int             	hashSize;

public:
	Iso(
	    int             _dimNumber,
	    const int*      _isotopeNumbers,
	    const int*      _atomCounts,
	    const double**  _isotopeMasses,
	    const double**  _isotopeProbabilities,
	    int		    _tabSize,
	    int             _hashSize
	);

	Iso(
	    const char* formula,
            int tabsize = 1000,
            int hashsize = 1000
        );


	virtual ~Iso();

	double getLightestPeakMass();
	double getHeaviestPeakMass();

};

 class IsoSpec : public Iso {
 protected:
     const double            cutOff;
     const std::vector<double>**     logProbs;
     const std::vector<double>**     masses;
     const std::vector<int*>**       marginalConfs;
     DirtyAllocator          allocator;
     std::vector<void*>      newaccepted;
     Summator                totalProb;
     unsigned int            cnt;
     int*                    candidate;
     void*                   topConf;
     int                     allDim;
     void*                   initialConf;

 public:
     IsoSpec(
         int             _dimNumber,
         const int*      _isotopeNumbers,
         const int*      _atomCounts,
         const double**  _isotopeMasses,
         const double**  _isotopeProbabilities,
         const double    _cutOff,
         int             tabSize = 1000,
         int             hashSize = 1000
     );

     static IsoThresholdGenerator* IsoFromFormula(
         const char* formula,
         double cutoff,
         int tabsize = 1000,
         int hashsize = 1000
     );

     virtual ~IsoSpec();

     virtual bool advanceToNextConfiguration() = 0;
     void processConfigurationsUntilCutoff();
     int getNoVisitedConfs();
     int getNoIsotopesTotal();


     void getCurrentProduct(double* res_mass, double* res_logProb, int* res_isoCounts);
     void getProduct(double* res_mass, double* res_logProb, int* res_isoCounts);
     std::tuple<double*,double*,int*,int> getCurrentProduct();
     std::tuple<double*,double*,int*,int> getProduct();

     #ifdef BUILDING_R
    // An R friend should be considered the worst enemy.
    //                              Sun Tzu.
    friend  NumericMatrix Rinterface(
         	const IntegerVector&  molecule,
         	const DataFrame&      isotopes,
         	double  stopCondition,
         	int		algo,
         	int 	tabSize,
         	int		hashSize,
         	double 	step,
         	bool 	showCounts,
            bool    trim
        );
     #endif

     friend class Spectrum;
 };

 class IsoSpecOrdered : public IsoSpec
 {
 protected:
     std::priority_queue<void*,std::vector<void*>,ConfOrder>  pq;

 public:
     IsoSpecOrdered(
         int             _dimNumber,
         const int*      _isotopeNumbers,
         const int*      _atomCounts,
         const double**  _isotopeMasses,
         const double**  _isotopeProbabilities,
         const double    _cutOff,
         int             tabSize = 1000,
         int             hashSize = 1000
     );

     virtual ~IsoSpecOrdered();

     bool advanceToNextConfiguration();


 };

 class IsoSpecLayered : public IsoSpec
 {
 protected:
     std::vector<void*>*         current;
     std::vector<void*>*         next;
     double                      lprobThr;
     double                      percentageToExpand;
     bool                        estimateThresholds;
     bool                        do_trim;
     int layers;
#ifdef DEBUG
     int moves = 0;
     int hits = 0;
#endif /* DEBUG */

 public:
     IsoSpecLayered(
         int             _dimNumber,
         const int*      _isotopeNumbers,
         const int*      _atomCounts,
         const double**  _isotopeMasses,
         const double**  _isotopeProbabilities,
         const double    _cutOff,
         int             tabSize = 1000,
         int             hashSize = 1000,
         double          layerStep = 0.3,
         bool            _estimateThresholds = false,
	 bool            trim = true
     );

     virtual ~IsoSpecLayered();

     bool advanceToNextConfiguration();
 };


 class IsoSpecThreshold : public IsoSpec
 {
 protected:
     std::vector<void*>          current;
     double                      lprobThr;
 public:
     IsoSpecThreshold(
         int             _dimNumber,
         const int*      _isotopeNumbers,
         const int*      _atomCounts,
         const double**  _isotopeMasses,
         const double**  _isotopeProbabilities,
         double          _threshold,
         bool            _absolute = true,
         int             tabSize = 1000,
         int             hashSize = 1000
     );

     virtual ~IsoSpecThreshold();

     bool advanceToNextConfiguration();
     void processConfigurationsAboveThreshold();
 };



class IsoGenerator : public Iso
{
public:
	virtual bool advanceToNextConfiguration() = 0;
	virtual const double& lprob() const = 0;
	virtual const double& mass() const = 0;
//	virtual const int* const & conf() const = 0;

	inline IsoGenerator(int _dimNumber,
            	const int*      _isotopeNumbers,
            	const int*      _atomCounts,
            	const double**  _isotopeMasses,
            	const double**  _isotopeProbabilities,
            	int             _tabSize,
            	int             _hashSize) : 
	    Iso(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses, _isotopeProbabilities, _tabSize, _hashSize) {};
	inline IsoGenerator(const char* formula, int _tabsize, int _hashsize) :
		Iso(formula, _tabsize, _hashsize) {}
	inline virtual ~IsoGenerator() {};

};

class IsoOrderedGenerator : public IsoGenerator
{
private:
	double cutOff;
	std::priority_queue<void*,std::vector<void*>,ConfOrder> pq;
	void IsoOrderedGenerator_init(double _cutoff);
	void* 				topConf;
	DirtyAllocator          	allocator;
        const std::vector<double>**     logProbs;
        const std::vector<double>**     masses;
        const std::vector<int*>**       marginalConfs;
	double 				currentLProb;
	double 				currentMass;
	int*				candidate;

public:
	virtual bool advanceToNextConfiguration();
	virtual const double& lprob() const;
	virtual const double& mass() const;

	inline IsoOrderedGenerator(int             _dimNumber,
                                   const int*      _isotopeNumbers,
                                   const int*      _atomCounts,
                                   const double**  _isotopeMasses,
                                   const double**  _isotopeProbabilities,
                                   const double    _cutOff   = std::numeric_limits<double>::infinity(),
                                   int             _tabSize  = 1000,
                                   int             _hashSize = 1000) : 
			IsoGenerator(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses, _isotopeProbabilities, _tabSize, _hashSize),
			allocator(dimNumber, _tabSize)
			{ IsoOrderedGenerator_init(_cutOff); };

	inline IsoOrderedGenerator(const char* formula,
                		   double  _cutOff   = std::numeric_limits<double>::infinity(),
		                   int     _tabSize  = 1000,
                		   int     _hashSize = 1000) : IsoGenerator(formula, _tabSize, _hashSize),
				   			       allocator(dimNumber, _tabSize)
				   			       { IsoOrderedGenerator_init(_cutOff); };

	virtual ~IsoOrderedGenerator();

};

class IsoThresholdGenerator : public IsoGenerator
{
private:
	unsigned int* counter;
	double* partialLProbs;
	double* partialMasses;
	double* maxConfsLPSum;
	double Lcutoff;

	void IsoThresholdGenerator_init(double _threshold, bool _absolute);

public:
	virtual bool advanceToNextConfiguration();
	virtual inline const double& lprob() const { return partialLProbs[0]; };
	virtual inline const double& mass() const { return partialMasses[0]; };
//	virtual const int* const & conf() const;

	inline IsoThresholdGenerator(int _dimNumber,
                const int*      _isotopeNumbers,
                const int*      _atomCounts,
                const double**  _isotopeMasses,
                const double**  _isotopeProbabilities,
		double          _threshold,
		bool            _absolute = true,
                int             _tabSize = 1000,
                int             _hashSize = 1000) : IsoGenerator(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses, _isotopeProbabilities, _tabSize, _hashSize)
						{ IsoThresholdGenerator_init(_threshold, _absolute); };

	inline IsoThresholdGenerator(const char* formula, 
	        double 	_threshold,
		bool 	_absolute = true,
		int 	_tabSize  = 1000,
		int 	_hashSize = 1000) : IsoGenerator(formula, _tabSize, _hashSize) { IsoThresholdGenerator_init(_threshold, _absolute); };

	inline virtual ~IsoThresholdGenerator() { delete[] counter; delete[] partialLProbs; delete[] partialMasses; delete[] maxConfsLPSum;};

private:
	inline void recalc(int idx)
	{
		for(; idx >=0; idx--)
		{
			partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->conf_probs()[counter[idx]]; 
			partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->conf_masses()[counter[idx]];
		}
	}

	

};


class IsoThresholdGeneratorMultithreaded : public IsoGenerator
{
private:
        unsigned int* counter;
	unsigned int total_threads;
	unsigned int thread_offset;
        double* partialLProbs;
        double* partialMasses;
        double* maxConfsLPSum;
        double Lcutoff;

        void IsoThresholdGeneratorMultithreaded_init(unsigned int _total_threads, unsigned int _thread_offset, double _threshold, bool _absolute);

public:
        virtual bool advanceToNextConfiguration();
        virtual inline const double& lprob() const { return partialLProbs[0]; };
        virtual inline const double& mass() const { return partialMasses[0]; };
//      virtual const int* const & conf() const;

        inline IsoThresholdGeneratorMultithreaded(
		unsigned int 	_total_threads,
		unsigned int 	_thread_offset,
		int 		_dimNumber,
                const int*      _isotopeNumbers,
                const int*      _atomCounts,
                const double**  _isotopeMasses,
                const double**  _isotopeProbabilities,
                double          _threshold,
                bool            _absolute = true,
                int             _tabSize = 1000,
                int             _hashSize = 1000) : IsoGenerator(_dimNumber, _isotopeNumbers, _atomCounts, _isotopeMasses, _isotopeProbabilities, _tabSize, _hashSize)
                                                { IsoThresholdGeneratorMultithreaded_init(_total_threads, _thread_offset, _threshold, _absolute); };

        inline IsoThresholdGeneratorMultithreaded(
		unsigned int    _total_threads,
                unsigned int    _thread_offset,
		const char* 	formula,
                double  	_threshold,
                bool    	_absolute = true,
                int     	_tabSize  = 1000,
                int     	_hashSize = 1000) : IsoGenerator(formula, _tabSize, _hashSize) 
						{ IsoThresholdGeneratorMultithreaded_init(_total_threads, _thread_offset, _threshold, _absolute); };

        inline virtual ~IsoThresholdGeneratorMultithreaded() { delete[] counter; delete[] partialLProbs; delete[] partialMasses; delete[] maxConfsLPSum;};

private:
        inline void recalc(int idx)
        {
                for(; idx >=0; idx--)
                {
                        partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->conf_probs()[counter[idx]];
                        partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->conf_masses()[counter[idx]];
                }
        }



};


 void printConfigurations(
     const   std::tuple<double*,double*,int*,int>& results,
     int     dimNumber,
     int*    isotopeNumbers
 );
#endif
