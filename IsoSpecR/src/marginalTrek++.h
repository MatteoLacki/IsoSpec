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


#ifndef MARGINALTREK_HPP
#define MARGINALTREK_HPP
#include <tuple>
#include <unordered_map>
#include <queue>
#include "conf.h"
#include "allocator.h"
#include "operators.h"
#include "summator.h"


Conf initialConfigure(const int atomCnt, const int isotopeNo, const double* probs);

std::tuple<double*,double*,int*,int> getMarginal(
    const double* masses,   // masses size = logProbs size = isotopeNo
    const double* probs,
    int isotopeNo,                  // No of isotope configurations.
    int atomCnt,
    const double cutOff = 0.99999,
    const int tabSize = 1000,
    const int hashSize= 1000
);

void printMarginal(const std::tuple<double*,double*,int*,int>& results, int dim);


class MarginalTrek
{
    int current_count;
    const int _tabSize;
    const int _hashSize;
public:
    const unsigned int _isotopeNo;
private:
    const unsigned int _atomCnt;
    const double* iso_masses;
    const double* logProbs;
    Allocator<int> allocator;
    KeyHasher keyHasher;
    const ConfEqual equalizer;
    const ConfOrderMarginal orderMarginal;
    std::unordered_map<Conf,int,KeyHasher,ConfEqual> visited;
    std::priority_queue<Conf,std::vector<Conf>,ConfOrderMarginal> pq;
    Summator totalProb;
    Conf candidate;
    std::vector<double> _conf_probs;
    std::vector<double> _conf_masses;
    std::vector<int*> _confs;

    bool add_next_conf();
    void sort_configs();

public:
    MarginalTrek(
        const double* masses,   // masses size = logProbs size = isotopeNo
        const double* probs,
        int isotopeNo,                  // No of isotope configurations.
        int atomCnt,
        int tabSize = 1000,
        int hashSize = 1000
    );

    inline bool probeConfigurationIdx(int idx)
    {
        if(current_count > idx)
            return true;
        while(current_count <= idx)
            if(not add_next_conf()) 
                return false;
        return true;
    }

    int processUntilCutoff(double cutoff);

    inline const std::vector<double>& conf_probs() const { return _conf_probs; };
    inline const std::vector<double>& conf_masses() const { return _conf_masses; };
    inline const std::vector<int*>& confs() const { return _confs; };
    inline int get_isotopeNo() const { return _isotopeNo; };
    ~MarginalTrek();

    double getLightestConfMass();
    double getHeaviestConfMass();
};






#endif
