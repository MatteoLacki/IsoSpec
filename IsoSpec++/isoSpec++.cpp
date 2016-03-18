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

#include <cmath>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <tuple>
#include <unordered_map>
#include <queue>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cctype>
#include "conf.hpp"
#include "dirtyAllocator.hpp"
#include "operators.hpp"
#include "summator.hpp"
#include "marginalTrek++.hpp"
#include "isoSpec++.hpp"
#include "misc.hpp"
#include "element_tables.h"



using namespace std;



IsoSpec::IsoSpec(
    int             _dimNumber,
    const int*      _isotopeNumbers,
    const int*      _atomCounts,
    const double**  isotopeMasses,
    const double**  isotopeProbabilities,
    const double    _cutOff,
    int             tabSize,
    int             hashSize
) :
dimNumber(_dimNumber),
cutOff(_cutOff),
allocator(_dimNumber, tabSize),
confSize(_dimNumber * sizeof(int)),
candidate(new int[dimNumber])

{
    isotopeNumbers  = new int[_dimNumber];
    memcpy(isotopeNumbers, _isotopeNumbers, _dimNumber*sizeof(int));
    atomCounts      = new int[_dimNumber];
    memcpy(atomCounts, _atomCounts, _dimNumber*sizeof(int));


    for(int i=0; i<dimNumber;i++) allDim += isotopeNumbers[i];

    logProbs        = new const vector<double>*[dimNumber];
    masses          = new const vector<double>*[dimNumber];
    marginalConfs   = new const vector<int*>*[dimNumber];

    marginalResults = new MarginalTrek*[dimNumber];
    for(int i = 0; i<dimNumber; i++)
    {
        marginalResults[i] = new MarginalTrek(
            isotopeMasses[i],
            isotopeProbabilities[i],
            isotopeNumbers[i],
            atomCounts[i],
            tabSize,
            hashSize
        );

        masses[i] = &marginalResults[i]->conf_masses();
        logProbs[i] = &marginalResults[i]->conf_probs();
        marginalConfs[i] = &marginalResults[i]->confs();
    }

    initialConf     = allocator.newConf();
    memset(
        reinterpret_cast<char*>(initialConf) + sizeof(double),
           0,
           sizeof(int)*dimNumber
    );

    *(reinterpret_cast<double*>(initialConf)) =
    combinedSum(
        getConf(initialConf),
                logProbs,
                dimNumber
    );

}

template<typename T> T* IsoSpec::IsoFromFormula(const char* formula, double cutoff, int tabsize, int hashsize)
{
    static_assert(std::is_base_of<IsoSpec, T>::value, "Template argument must be derived from IsoSpec");

    string cpp_formula(formula);
    int last_modeswitch = 0;
    int mode = 0;
    int pos = 0;
    std::vector<string> elements;
    std::vector<int> numbers;
    while(formula[pos] != '\0')
    {
        if(isdigit(formula[pos]) && mode == 0)
        {
            elements.push_back(cpp_formula.substr(last_modeswitch, pos-last_modeswitch));
            last_modeswitch = pos;
            mode = 1;
        }
        else if(isalpha(formula[pos]) && mode == 1)
        {
            numbers.push_back(stoi(cpp_formula.substr(last_modeswitch, pos-last_modeswitch)));
            last_modeswitch = pos;
            mode = 0;
        }
        pos++;
    }

    numbers.push_back(stoi(cpp_formula.substr(last_modeswitch, pos)));


    if(elements.size() != numbers.size())
        throw invalid_argument("Invalid formula");

    std::vector<int> element_indexes;

    for (unsigned int i=0; i<elements.size(); i++)
    {
        int idx = -1;
        for(int j=0; j<NUMBER_OF_ISOTOPIC_ENTRIES; j++)
        {
            if (elements[i].compare(elem_table_symbol[j]) == 0)
            {
                idx = j;
                break;
            }
        }
        if(idx < 0)
            throw invalid_argument("Invalid formula");
        element_indexes.push_back(idx);

    }

    vector<int> isotope_numbers;

    for(auto it = element_indexes.begin(); it != element_indexes.end(); ++it)
    {
        int num = 0;
        int at_idx = *it;
        int atomicNo = elem_table_atomicNo[at_idx];
        while(at_idx < NUMBER_OF_ISOTOPIC_ENTRIES && elem_table_atomicNo[at_idx] == atomicNo)
        {
            at_idx++;
            num++;
        }
        isotope_numbers.push_back(num);
    }

    vector<const double*> isotope_masses;
    vector<const double*> isotope_probabilities;
    for(auto it = element_indexes.begin(); it != element_indexes.end(); ++it)
    {
        isotope_masses.push_back(&elem_table_mass[*it]);
        isotope_probabilities.push_back(&elem_table_probability[*it]);
    }

    return new T(
        elements.size(),
                 isotope_numbers.data(),
                 numbers.data(),
                 isotope_masses.data(),
                 isotope_probabilities.data(),
                 cutoff,
                 tabsize,
                 hashsize
    );

}

template  IsoSpecOrdered*
IsoSpec::IsoFromFormula<IsoSpecOrdered>(const char* formula,
                                        double cutoff,
                                        int tabsize,
                                        int hashsize);
template  IsoSpecLayered*
IsoSpec::IsoFromFormula<IsoSpecLayered>(const char* formula,
                                        double cutoff,
                                        int tabsize,
                                        int hashsize);


void IsoSpec::processConfigurationsUntilCutoff()
{
    while( cutOff > totalProb.get() && advanceToNextConfiguration() ) {}
}


bool IsoSpecOrdered::advanceToNextConfiguration()
{
    if(pq.size() < 1)
        return false;


    topConf = pq.top();
    pq.pop();

    cnt++;

    int* topConfIsoCounts = getConf(topConf);

    //      visited[topConfIsoCounts] = cnt;
    newaccepted.push_back(topConf);

    totalProb.add( exp(*reinterpret_cast<double*>(topConf) ));

    for(int j = 0; j < dimNumber; ++j)
    {
        // candidate cannot refer to a position that is
        // out of range of the stored marginal distribution.
        if(marginalResults[j]->probeConfigurationIdx(topConfIsoCounts[j] + 1))
        {
            memcpy(candidate, topConfIsoCounts, confSize);
            candidate[j]++;

            void*       acceptedCandidate                       = allocator.newConf();
            int*        acceptedCandidateIsoCounts      = getConf(acceptedCandidate);
            memcpy(     acceptedCandidateIsoCounts, candidate, confSize);

            *(reinterpret_cast<double*>(acceptedCandidate)) =
            combinedSum(
                candidate,
                logProbs,
                dimNumber
            );
            pq.push(acceptedCandidate);
        }
        if(topConfIsoCounts[j] > 0)
            break;
    }


    return true;
}



std::tuple<double*,double*,int*,int> IsoSpec::getProduct()
{
    processConfigurationsUntilCutoff();
    return getCurrentProduct();
}


std::tuple<double*,double*,int*,int> IsoSpec::getCurrentProduct()
{

    double*         res_mass        = new double[cnt];
    double*         res_logProb     = new double[cnt];
    int*            res_isoCounts   = new int[cnt*allDim];

    getProduct(res_mass, res_logProb, res_isoCounts);

    return std::tuple<double*,double*,int*,int>(
        res_mass,
        res_logProb,
        res_isoCounts,
        cnt
    );
}


void IsoSpec::getProduct(double* res_mass, double* res_logProb, int* res_isoCounts)
{
    processConfigurationsUntilCutoff();
    getCurrentProduct(res_mass, res_logProb, res_isoCounts);
}

void IsoSpec::getCurrentProduct(double* res_mass, double* res_logProb, int* res_isoCounts)
{

    int i = 0;
    int j = 0;

    for(auto it = newaccepted.cbegin(); it != newaccepted.cend(); it++)
    {
        int* curr_conf  = getConf(*it);
        res_mass[i]     = combinedSum( curr_conf, masses, dimNumber );
        res_logProb[i]  = getLProb(*it);

        for(int isotopeNumber=0; isotopeNumber<dimNumber; isotopeNumber++)
        {
            int currentConfIndex = curr_conf[isotopeNumber];
            int locIsoNo = isotopeNumbers[isotopeNumber];
            memcpy(
                &res_isoCounts[j],
                (*marginalConfs[isotopeNumber])[currentConfIndex],
                   sizeof(int)*locIsoNo
            );
            j += locIsoNo;
        }
        i++;
    }
}
/*
std::unordered_map<double, double> IsoSpec::getPlot(double resolution)
{
    xret = new unordered_map<double, double>

    for(auto it = newaccepted.cbegin(); it != newaccepted.cend(); it++)
    {
        int* curr_conf  = getConf(*it);
    double mass     = nearbyint(combinedSum( curr_conf, masses, dimNumber )/resolution)*resolution;
    double logProb  = getLProb(*it);
        
    if(xret.count(mass) == 0)
        xret[mass] = new Summator();

    xret[mass].add(exp(logProb));
    }

    return xret;
}
*/
int IsoSpec::getNoVisitedConfs()
{
    return newaccepted.size();
}

int IsoSpec::getNoIsotopesTotal()
{
    int ret = 0;
    for(int i=0; i<dimNumber; i++)
        ret += isotopeNumbers[i];
    return ret;
}

template<typename T> void dealloc_table(T* tbl, int dim)
{
    for(int i=0; i<dim; i++)
    {
        delete tbl[i];
    }
    delete[] tbl;
}

IsoSpec::~IsoSpec()
{
    dealloc_table(marginalResults, dimNumber);
    delete[] candidate;
    delete[] logProbs;
    delete[] masses;
    delete[] marginalConfs;
    delete[] atomCounts;
    delete[] isotopeNumbers;
}

IsoSpecOrdered::IsoSpecOrdered( int             _dimNumber,
                                const int*      _isotopeNumbers,
                                const int*      _atomCounts,
                                const double**  _isotopeMasses,
                                const double**  _isotopeProbabilities,
                                const double    _cutOff,
                                int             tabSize,
                                int             hashSize
) : IsoSpec( _dimNumber,
             _isotopeNumbers,
             _atomCounts,
             _isotopeMasses,
             _isotopeProbabilities,
             _cutOff,
             tabSize = 1000,
             hashSize = 1000)
{pq.push(initialConf);};


IsoSpecOrdered::~IsoSpecOrdered(){};


IsoSpecLayered::IsoSpecLayered( int             _dimNumber,
                                const int*      _isotopeNumbers,
                                const int*      _atomCounts,
                                const double**  _isotopeMasses,
                                const double**  _isotopeProbabilities,
                                const double    _cutOff,
                                int             tabSize,
                                int             hashSize,
                                double          layerStep,
                bool            _estimateThresholds
) : IsoSpec( _dimNumber,
             _isotopeNumbers,
             _atomCounts,
             _isotopeMasses,
             _isotopeProbabilities,
             _cutOff,
             tabSize = 1000,
             hashSize = 1000
),
estimateThresholds(_estimateThresholds)
{
    current = new std::vector<void*>();
    next    = new std::vector<void*>();

    current->push_back(initialConf);

    percentageToExpand = layerStep;
    lprobThr = (*reinterpret_cast<double*>(initialConf));
};


IsoSpecLayered::~IsoSpecLayered()
{
    if(current != NULL)
        delete current;
    if(next != NULL)
        delete next;
};

bool IsoSpecLayered::advanceToNextConfiguration()
{
    layers += 1;
    double maxFringeLprob = -std::numeric_limits<double>::infinity();

    if(current == nullptr)
        return false;
    int accepted_in_this_layer = 0;
    Summator prob_in_this_layer(totalProb);


    while(current->size() > 0)
    {
        topConf = current->back();
        current->pop_back();

        cnt++;

        double top_lprob = getLProb(topConf);

        if(top_lprob >= lprobThr)
        {
#ifdef DEBUG
            hits += 1;
#endif /* DEBUG */
            newaccepted.push_back(topConf);
            accepted_in_this_layer++;
            prob_in_this_layer.add(exp(top_lprob));
        }
        else
        {
#ifdef DEBUG
            moves += 1;
#endif /* DEBUG */
            next->push_back(topConf);
            continue;
        }

    int* topConfIsoCounts = getConf(topConf);

        for(int j = 0; j < dimNumber; ++j)
        {
            // candidate cannot refer to a position that is
            // out of range of the stored marginal distribution.
            if(marginalResults[j]->probeConfigurationIdx(topConfIsoCounts[j] + 1))
            {
                memcpy(candidate, topConfIsoCounts, confSize);
                candidate[j]++;

                void*       acceptedCandidate                       = allocator.newConf();
                int*        acceptedCandidateIsoCounts      = getConf(acceptedCandidate);
                memcpy(     acceptedCandidateIsoCounts, candidate, confSize);

                double newConfProb = combinedSum(
                    candidate,
                    logProbs,
                    dimNumber
                );



                *(reinterpret_cast<double*>(acceptedCandidate)) = newConfProb;

                if(newConfProb >= lprobThr)
                    current->push_back(acceptedCandidate);
                else
        {
                    next->push_back(acceptedCandidate);
            if(newConfProb > maxFringeLprob)
                maxFringeLprob = top_lprob;
        }
            }
            if(topConfIsoCounts[j] > 0)
                break;
        }
    }

    if(next == nullptr || next->size() < 1)
        return false;
    else
    {
        if(prob_in_this_layer.get() < cutOff)
        {
#ifdef DEBUG
            Summator testDupa(prob_in_this_layer);
            for (auto it = next->begin(); it != next->end(); it++) {
                testDupa.add(exp(getLProb(*it)));
            }
            std::cout << "Prob(Layer) = " << prob_in_this_layer.get() << std::endl;
            std::cout << "Prob(Layer)+Prob(Fringe) = " << testDupa.get() << std::endl;
            std::cout << "Layers = " << layers << std::endl;
            std::cout << std::endl;
#endif /* DEBUG */

        // // This was an attempt to merge two methods: layered and layered_estimating 
        // // that does not work so good as predicted.
//             if( estimateThresholds and ( prob_in_this_layer.get() >= cutOff*.99 ) ){
//                 estimateThresholds = false;
//                 percentageToExpand = .25; // The ratio of one rectangle to the rectangle.
// #ifdef DEBUG
//                 std::cout << "We switch!" << std::endl;
// #endif /* DEBUG */
//             }

#ifdef DEBUG
                std::cout << "percentageToExpand = " << percentageToExpand << std::endl;
#endif /* DEBUG */

            std::vector<void*>* nnew = current;
            nnew->clear();
            current = next;
            next = nnew;
            int howmany = floor(current->size()*percentageToExpand);
            if(estimateThresholds){
                // Screw numeric correctness, ARRRRRRR!!! Well, this is an estimate anyway, doesn't have to be that precise
                // lprobThr += log((1.0-cutOff))+log1p((percentageToExpand-1.0)/layers) - log(1.0-prob_in_this_layer.get());
                lprobThr += log(1.0-cutOff) + log(1.0-(1.0-percentageToExpand)/pow(layers, 2.0)) - log(1.0 -prob_in_this_layer.get());
                if(lprobThr > maxFringeLprob){
                    lprobThr = maxFringeLprob;
                    estimateThresholds = false;
                    percentageToExpand = .3;
#ifdef DEBUG
                    std::cout << "We switch to other method because density estimates where higher than max on fringe." << std::endl;
#endif /* DEBUG */
                    lprobThr = getLProb(quickselect(current->data(), howmany, 0, current->size()));                    
                }
            } else
                lprobThr = getLProb(quickselect(current->data(), howmany, 0, current->size()));
            totalProb = prob_in_this_layer;
        }
        else
        {
#ifdef DEBUG
            std::cerr << "No. layers: " << layers << "  hits: " << hits << "    misses: " << moves << " miss ratio: " << static_cast<double>(moves) / static_cast<double>(hits) << std::endl;
#endif /* DEBUG */
            delete next;
            next = nullptr;
            delete current;
            current = nullptr;
            int start = 0;
            int end = accepted_in_this_layer - 1;
            void* swapspace;

            void** lastLayer = &(newaccepted.data()[newaccepted.size()-accepted_in_this_layer]);

            Summator qsprob(totalProb);
            while(totalProb.get() < cutOff)
            {
                if(start == end)
                    break;

                // Partition part

                int len = end - start;
                int pivot = rand() % len + start;
                void* pval = lastLayer[pivot];
                double pprob = getLProb(pval);
                mswap(lastLayer[pivot], lastLayer[end-1]);
                int loweridx = start;
                for(int i=start; i<end-1; i++)
                {
                    if(getLProb(lastLayer[i]) > pprob)
                    {
                        mswap(lastLayer[i], lastLayer[loweridx]);
                        loweridx++;
                    }
                }
                mswap(lastLayer[end-1], lastLayer[loweridx]);

                // Selection part

                Summator leftProb(qsprob);
                for(int i=start; i<=loweridx; i++)
                {
                    leftProb.add(exp(getLProb(lastLayer[i])));
                }
                if(leftProb.get() < cutOff)
                {
                    start = loweridx+1;
                    qsprob = leftProb;
                }
                else
                    end = loweridx;
            }
        int accend = newaccepted.size()-accepted_in_this_layer+start+1;
#ifdef DEBUG
            std::cerr << "Last layer size: " << accepted_in_this_layer << " Total size: " << newaccepted.size() << "    Total size after trimming: " << accend << " No. trimmed: " << -start-1+accepted_in_this_layer 
        << "    Trimmed to left ratio: " << static_cast<double>(-start-1+accepted_in_this_layer) / static_cast<double>(accend) << std::endl;
#endif /* DEBUG */

            totalProb = qsprob;
            newaccepted.resize(accend);
            return true;
        }
    }
    return true;

};




IsoSpecThreshold::IsoSpecThreshold( int             _dimNumber,
                                    const int*      _isotopeNumbers,
                                    const int*      _atomCounts,
                                    const double**  _isotopeMasses,
                                    const double**  _isotopeProbabilities,
                                    double          _threshold,
                                    bool            _absolute,
                                    int             tabSize,
                                    int             hashSize
) : IsoSpec( _dimNumber,
             _isotopeNumbers,
             _atomCounts,
             _isotopeMasses,
             _isotopeProbabilities,
             0.0,
             tabSize = 1000,
             hashSize = 1000
)
{
    current.push_back(initialConf);

    lprobThr = log(_threshold);

    if(not _absolute)
        lprobThr += *reinterpret_cast<double*>(initialConf);
};


IsoSpecThreshold::~IsoSpecThreshold(){};


bool IsoSpecThreshold::advanceToNextConfiguration()
{
    if(current.size() < 1)
        return false;


    if(current.size() > 0)
    {
        topConf = current.back();
        current.pop_back();

        cnt++;

        int* topConfIsoCounts = getConf(topConf);

        double top_lprob = getLProb(topConf);

        if(top_lprob >= lprobThr)
        {
            newaccepted.push_back(topConf);
            totalProb.add(exp(top_lprob));
        }
        else
        {
            return true;
        }

        memcpy(candidate, topConfIsoCounts, confSize);

        for(int j = 0; j < dimNumber; ++j)
        {
            // candidate cannot refer to a position that is
            // out of range of the stored marginal distribution.

            if(marginalResults[j]->probeConfigurationIdx(topConfIsoCounts[j] + 1))
            {
                candidate[j]++;

                double newConfProb = combinedSum(
                    candidate,
                    logProbs,
                    dimNumber
                );

                if(newConfProb >= lprobThr)
                {
                    void*       acceptedCandidate               = allocator.newConf();
                    int*        acceptedCandidateIsoCounts      = getConf(acceptedCandidate);

                    memcpy(acceptedCandidateIsoCounts, candidate, confSize);


                    *(reinterpret_cast<double*>(acceptedCandidate)) = newConfProb;

                    current.push_back(acceptedCandidate);
                }

                candidate[j]--;
            }
            if(topConfIsoCounts[j] > 0)
                break;
        }
    }

    return true;

};


void IsoSpecThreshold::processConfigurationsAboveThreshold()
{
    while (advanceToNextConfiguration()) {}
}



void printConfigurations(
    const   std::tuple<double*,double*,int*,int>& results,
    int     dimNumber,
    int*    isotopeNumbers
){
    int m = 0;

    for(int i=0; i<std::get<3>(results); i++){

        std::cout << "Mass = "  << std::get<0>(results)[i] <<
        "\tand log-prob = "     << std::get<1>(results)[i] <<
        "\tand prob = "                 << exp(std::get<1>(results)[i]) <<
        "\tand configuration =\t";


        for(int j=0; j<dimNumber; j++){
            for(int k=0; k<isotopeNumbers[j]; k++ )
            {
                std::cout << std::get<2>(results)[m] << " ";
                m++;
            }
            std::cout << '\t';
        }


        std::cout << std::endl;
    }
}






