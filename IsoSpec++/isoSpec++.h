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
#include <limits>
#include "lang.h"
#include "dirtyAllocator.h"
#include "summator.h"
#include "operators.h"
#include "marginalTrek++.h"


#ifdef BUILDING_R
 #include <Rcpp.h>
 using namespace Rcpp;
#endif /* BUILDING_R */



unsigned int parse_formula(const char* formula, std::vector<const double*>& isotope_masses, std::vector<const double*>& isotope_probabilities, int** isotopeNumbers, int** atomCounts, unsigned int* confSize);

class IsoSpecLayered;
class IsoThresholdGenerator;

class Iso {
private:
	void setupMarginals(const double** _isotopeMasses, const double** _isotopeProbabilities);
public:
        bool disowned;
protected:
	int 			dimNumber;
	int*			isotopeNumbers;
	int*			atomCounts;
	unsigned int		confSize;
	int			allDim;
	Marginal**              marginals;
        double                  modeLProb;

public:
	Iso(
	    int             _dimNumber,
	    const int*      _isotopeNumbers,
	    const int*      _atomCounts,
	    const double**  _isotopeMasses,
	    const double**  _isotopeProbabilities
	);

	Iso(const char* formula);

        Iso(Iso&& other);

        Iso(const Iso& other, bool fullcopy);

	virtual ~Iso();

	double getLightestPeakMass() const;
	double getHeaviestPeakMass() const;
        inline double getModeLProb() const { return modeLProb; };
        inline int getDimNumber() const { return dimNumber; };
        PrecalculatedMarginal** get_MT_marginal_set(double Lcutoff, bool absolute, int tabSize, int hashSize);

};

 class IsoSpec: public Iso {
 protected:
     MarginalTrek** marginalResults;
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
         int             tabSize = 1000
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
         double          layerStep = 0.3,
         bool            _estimateThresholds = false,
	 bool            trim = true
     );

     virtual ~IsoSpecLayered();

     bool advanceToNextConfiguration();
 };


// Be very absolutely safe vs. false-sharing cache lines between threads...
#define PADDING 64

class IsoGenerator : public Iso
{
protected:
        double* partialLProbs;
        double* partialMasses;
        double* partialExpProbs;

public:
	virtual bool advanceToNextConfiguration() = 0;
        inline const double& lprob() const { return partialLProbs[0]; };
        inline const double& mass()  const { return partialMasses[0]; };
        inline const double& eprob() const { return partialExpProbs[0]; };
//	virtual const int* const & conf() const = 0;


    inline IsoGenerator(Iso&& iso);
	inline virtual ~IsoGenerator() { delete[] partialLProbs; delete[] partialMasses; delete[] partialExpProbs; };

};

class IsoOrderedGenerator: public IsoGenerator
{
private:
    MarginalTrek** marginalResults;
    std::priority_queue<void*,std::vector<void*>,ConfOrder> pq;
    void* topConf;
    DirtyAllocator allocator;
    const std::vector<double>**     logProbs;
    const std::vector<double>**     masses;
    const std::vector<int*>**       marginalConfs;
    double currentLProb;
    double currentMass;
    int*   candidate;

public:
    virtual bool advanceToNextConfiguration();
    virtual const double& lprob() const { return currentLProb; };
    virtual const double& mass() const { return currentMass; };
    virtual inline void get_conf_signature(int* space) const 
    { 
        const int* c = getConf(topConf);
        for(int ii=0; ii<dimNumber; ii++)
        {
            memcpy(space, marginalResults[ii]->confs()[c[ii]], isotopeNumbers[ii]*sizeof(int));
            space += isotopeNumbers[ii];
        }
    };

    IsoOrderedGenerator(Iso&& iso, int _tabSize  = 1000, int _hashSize = 1000);

    virtual ~IsoOrderedGenerator();

};

class IsoThresholdGenerator: public IsoGenerator
{
private:
    int* counter;
    double* maxConfsLPSum;
    const double Lcutoff;
    PrecalculatedMarginal** marginalResults;

public:
    virtual bool advanceToNextConfiguration();
    virtual inline void get_conf_signature(int* space) 
    {
        for(int ii=0; ii<dimNumber; ii++)
        {
            memcpy(space, marginalResults[ii]->get_conf(counter[ii]), isotopeNumbers[ii]*sizeof(int));
            space += isotopeNumbers[ii];
        }
    };

    IsoThresholdGenerator(Iso&& iso, double _threshold, bool _absolute=true,
                          int _tabSize=1000, int _hashSize=1000);

    inline virtual ~IsoThresholdGenerator() { delete[] counter;
                                              delete[] maxConfsLPSum;
                                              dealloc_table(marginalResults, dimNumber); };

    void terminate_search();

private:
    inline void recalc(int idx)
    {
	for(; idx >=0; idx--)
	{
	    partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
	    partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
            partialExpProbs[idx] = partialExpProbs[idx+1] * marginalResults[idx]->get_eProb(counter[idx]);
	}
    }



};


class IsoThresholdGeneratorBoundMass : public IsoGenerator
{
private:
	double* maxConfsLPSum, *minMassCSum, *maxMassCSum;
	double Lcutoff;
        double min_mass, max_mass;
        RGTMarginal** marginalResults;

public:
    //TODO work on virtuals
	virtual bool advanceToNextConfiguration();
//        virtual inline void get_conf_signature(unsigned int* target);
//	virtual const int* const & conf() const;

        IsoThresholdGeneratorBoundMass(Iso&& iso, double  _threshold, double min_mass, double max_mass, bool _absolute = true, int _tabSize  = 1000, int _hashSize = 1000);

	virtual ~IsoThresholdGeneratorBoundMass();

private:
	void setup_ith_marginal_range(unsigned int idx);
        inline void recalc(int idx)
        {
            partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->current_lProb();
            partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->current_mass();
            partialExpProbs[idx] = partialExpProbs[idx+1] * marginalResults[idx]->current_eProb();
        }




};



class IsoThresholdGeneratorMT : public IsoGenerator
{
private:
	unsigned int* counter;
	double* maxConfsLPSum;
	const double Lcutoff;
        SyncMarginal* last_marginal;
        PrecalculatedMarginal** marginalResults;

public:
	virtual bool advanceToNextConfiguration();
        virtual inline void get_conf_signature(unsigned int* target) { memcpy(target, counter, sizeof(unsigned int)*dimNumber); };
//	virtual const int* const & conf() const;

        IsoThresholdGeneratorMT(Iso&& iso, double  _threshold, PrecalculatedMarginal** marginals, bool _absolute = true);

	inline virtual ~IsoThresholdGeneratorMT() { delete[] counter; delete[] maxConfsLPSum;};
        void terminate_search();

private:
	inline void recalc(int idx)
	{
		for(; idx >=0; idx--)
		{
			partialLProbs[idx] = partialLProbs[idx+1] + marginalResults[idx]->get_lProb(counter[idx]);
			partialMasses[idx] = partialMasses[idx+1] + marginalResults[idx]->get_mass(counter[idx]);
                        partialExpProbs[idx] = partialExpProbs[idx+1] * marginalResults[idx]->get_eProb(counter[idx]);
		}
	}



};





#ifndef BUILDING_R

 void printConfigurations(
     const   std::tuple<double*,double*,int*,int>& results,
     int     dimNumber,
     int*    isotopeNumbers
 );
#endif
#endif /* ISOSPEC_PLUS_PLUS_HPP */
