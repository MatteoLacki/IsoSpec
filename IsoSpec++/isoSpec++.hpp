/*
 *   Copyright (C) 2015 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License
 *   version 3, as published by the Free Software Foundation.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with IsoSpec.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ISOSPEC_PLUS_PLUS_HPP
#define ISOSPEC_PLUS_PLUS_HPP

#include <tuple>
#include <unordered_map>
#include <queue>
#include "dirtyAllocator.hpp"
#include "summator.hpp"
#include "operators.hpp"
#include "marginalTrek++.hpp"


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
     int                     cnt = 0;
     const int               confSize;
     int*                    candidate;
     void*                   topConf;
     int                     allDim = 0;
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


     void getCurrentProduct(double* res_mass, double* res_logProb, int* res_isoCounts);
     void getProduct(double* res_mass, double* res_logProb, int* res_isoCounts);
     std::tuple<double*,double*,int*,int> getCurrentProduct();
     std::tuple<double*,double*,int*,int> getProduct();

     #ifdef BUILDING_R
     friend List Rinterface(
         IntegerVector isotopeNumbers,
         IntegerVector atomCounts,
         NumericVector isotopeMasses,
         NumericVector isotopeProbabilities,
         double stopCondition, int algo, int tabSize, int hashSize, double step);
     #endif
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
     int layers = 0;
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
	     bool            _estimateThresholds = false
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
