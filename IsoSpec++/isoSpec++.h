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

 class IsoSpec{
 protected:
     const int               dimNumber;
     int*                    isotopeNumbers;
     int*                    atomCounts;
     const double            cutOff;
     MarginalTrek**          marginalResults;
     const std::vector<double>**     logProbs;
     const std::vector<double>**     masses;
     const std::vector<int*>**       marginalConfs;
     DirtyAllocator          allocator;
     std::vector<void*>      newaccepted;
     Summator                totalProb;
     unsigned int            cnt;
     const unsigned int      confSize;
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

     template<typename T = IsoSpecLayered> static T* IsoFromFormula(
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
     double getLightestPeakMass();
     double getHeaviestPeakMass();


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



 void printConfigurations(
     const   std::tuple<double*,double*,int*,int>& results,
     int     dimNumber,
     int*    isotopeNumbers
 );
#endif
