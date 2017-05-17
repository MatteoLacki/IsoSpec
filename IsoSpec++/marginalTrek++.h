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
#include <atomic>
#include "conf.h"
#include "allocator.h"
#include "operators.h"
#include "summator.h"


Conf initialConfigure(const int atomCnt, const int isotopeNo, const double* probs);


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
    const KeyHasher keyHasher;
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
//    	std::cerr << "PC\n" << idx << '\n';
        while(current_count <= idx)
            if(not add_next_conf()) 
	    {	//std::cerr << "PC\n" << idx << "FALSE" << '\n';
                return false;
	    }
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



class PrecalculatedMarginal
{
protected:
    std::vector<Conf> configurations;
    Conf* confs;
    int no_confs;
    double* masses;
    double* lProbs;
    const unsigned int isotopeNo;
    const unsigned int atomCnt;
    const double* isoMasses;
    const double* isoLProbs;
    Allocator<int> allocator;
public: 
    PrecalculatedMarginal(
    	const double* masses,
	const double* probs,
	int isotopeNo,
	int atomCnt,
	double cutOff,
	bool sort = true,
	int tabSize = 1000,
	int hashSize = 1000
    );
    virtual ~PrecalculatedMarginal();
    inline bool inRange(int idx) { return idx < no_confs; };
};

class SyncMarginal : public PrecalculatedMarginal
{
private:
    std::atomic<int> counter;
public:
    inline SyncMarginal(
        const double* masses,
        const double* probs,
        int isotopeNo,
        int atomCnt,
        double cutOff,
        int tabSize = 1000,
        int hashSize = 1000
    ) : PrecalculatedMarginal(
        masses,
        probs,
        isotopeNo,
        atomCnt,
        cutOff,
        false,
        tabSize,
        hashSize
    ), counter(0) {};


    inline int getNextConfIdx() { return counter.fetch_add(1, std::memory_order_relaxed); };

};

class RGTMarginal : public PrecalculatedMarginal
{
private:
    const unsigned int tree_overhead;
    const unsigned int tree_size;
    const double* probs_tree;
    const unsigned int* subtree_sizes;
    const unsigned int mass_layer_size;
    const TableOrder<double> mass_order;
    const unsigned int* subtree_locations;
    const double* mass_table;
    const unsigned int* subintervals;
public:
    RGTMarginal(
	const double* masses,
        const double* probs,
        int isotopeNo,
        int atomCnt,
        double cutOff,
        int tabSize = 1000,
        int hashSize = 1000
    );
private:
    double* alloc_and_construct_ptree();
    void construct_ptree(double* new_tree, unsigned int where, double* pstart, unsigned int howmany);
    unsigned int compute_total_no_masses();
    unsigned int* alloc_and_setup_subtree_sizes();
    unsigned int setup_subtree_sizes(unsigned int* T, unsigned int idx);
    unsigned int* alloc_and_setup_subtree_locations();
    unsigned int setup_subtree_locations(unsigned int* T, unsigned int idx, unsigned int csum);
    unsigned int* alloc_and_setup_subintervals();
    unsigned int setup_subintervals(unsigned int* T, unsigned int idx, bool left);
    double* alloc_and_setup_mass_table();




};

#endif
