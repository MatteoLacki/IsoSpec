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
#include <stdexcept>
#include <string>
#include <limits>
#include "platform.h"
#include "conf.h"
#include "dirtyAllocator.h"
#include "operators.h"
#include "summator.h"
#include "marginalTrek++.h"
#include "isoSpec++.h"
#include "misc.h"
#include "element_tables.h"


using namespace std;





IsoLayered::IsoLayered( Iso&&     iso,
                        double    _cutOff,
                        double    _delta,
                        int       _tabSize,
                        int       _hashSize,
                        bool      trim
) : Iso(std::move(iso)),
allocator(dimNumber, _tabSize),
candidate(new int[dimNumber]),
cutOff(_cutOff),
do_trim(trim),
layers(0)
{
    marginalResults = new MarginalTrek*[dimNumber];

    for(int i = 0; i<dimNumber; i++)
        marginalResults[i] = new MarginalTrek(std::move(*(marginals[i])), _tabSize, _hashSize);

    logProbs        = new const vector<double>*[dimNumber];
    masses          = new const vector<double>*[dimNumber];
    marginalConfs   = new const vector<int*>*[dimNumber];

    for(int i = 0; i<dimNumber; i++)
    {
        masses[i] = &marginalResults[i]->conf_masses();
        logProbs[i] = &marginalResults[i]->conf_lprobs();
        marginalConfs[i] = &marginalResults[i]->confs();
    }

    void* topConf = allocator.newConf();
    bzero(reinterpret_cast<char*>(topConf) + sizeof(double), sizeof(int)*dimNumber);
    *(reinterpret_cast<double*>(topConf)) = combinedSum(getConf(topConf), logProbs, dimNumber);

    current = new std::vector<void*>();
    next    = new std::vector<void*>();

    current->push_back(topConf);

    lprobThr = (*reinterpret_cast<double*>(topConf));
}


IsoLayered::~IsoLayered()
{
    if(current != nullptr)
        delete current;
    if(next != nullptr)
        delete next;
    delete[] logProbs;
    delete[] masses;
    delete[] marginalConfs;
    delete[] candidate;
    dealloc_table(marginalResults, dimNumber);
}

bool IsoLayered::advanceToNextLayer()
{
    layers += 1;
    double maxFringeLprob = -std::numeric_limits<double>::infinity();

    if(current == nullptr)
        return false;
    int accepted_in_this_layer = 0;
    Summator prob_in_this_layer(totalProb);

    void* topConf;

    while(current->size() > 0)
    {
        topConf = current->back();
        current->pop_back();

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

                void*       acceptedCandidate          = allocator.newConf();
                int*        acceptedCandidateIsoCounts = getConf(acceptedCandidate);
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
            for (std::vector<void*>::iterator it = next->begin(); it != next->end(); it++) {
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

            if(do_trim)
            {
                void** lastLayer = &(newaccepted.data()[newaccepted.size()-accepted_in_this_layer]);

                Summator qsprob(totalProb);
                while(totalProb.get() < cutOff)
                {
                    if(start == end)
                        break;

                    // Partition part

                    int len = end - start;
#ifdef BUILDING_R
            int pivot = len/2 + start;  // We're very definitely NOT switching to R to use a RNG, and if R sees us use C RNG it complains...
#else
            int pivot = rand() % len + start;
#endif
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
            else // No trimming
            {
                totalProb = prob_in_this_layer;
                return true;
            }
        }
    }
    return true;

}




