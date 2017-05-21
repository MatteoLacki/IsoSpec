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
    unsigned int no_confs;
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
	double lCutOff,
	bool sort = true,
	int tabSize = 1000,
	int hashSize = 1000
    );
    virtual ~PrecalculatedMarginal();
    inline bool inRange(unsigned int idx) const { return idx < no_confs; };
    inline const double& get_lProb(unsigned int idx) const { return lProbs[idx]; };
    inline const double& get_mass(unsigned int idx) const { return masses[idx]; };
    inline const Conf& get_conf(unsigned int idx) const { return confs[idx]; };
};

class SyncMarginal : public PrecalculatedMarginal
{
private:
    std::atomic<unsigned int> counter;
public:
    inline SyncMarginal(
        const double* masses,
        const double* probs,
        int isotopeNo,
        int atomCnt,
        double lCutOff,
        int tabSize = 1000,
        int hashSize = 1000
    ) : PrecalculatedMarginal(
        masses,
        probs,
        isotopeNo,
        atomCnt,
        lCutOff,
        false,
        tabSize,
        hashSize
    ), counter(0) {};


    inline unsigned int getNextConfIdx() { return counter.fetch_add(1, std::memory_order_relaxed); };
    inline unsigned int getNextConfIdxwMass(double mmin, double mmax)
    {
    	unsigned int local = counter.fetch_add(1, std::memory_order_relaxed);
	while(local < no_confs and (mmin > masses[local] or mmax < masses[local]))
	    local = counter.fetch_add(1, std::memory_order_relaxed);
	return local;
    }

};

class RGTMarginal : public PrecalculatedMarginal
{
private:
    const unsigned int TOV_NP2, TOV_NP2M1;
    const unsigned int mass_table_rows_no, mass_table_row_size, mass_table_size;
    unsigned int* subintervals;
    double* mass_table;
    double pmin, pmax, mmin, mmax;
    unsigned int lower, upper, arridx, arrend, mask, current_level, cidx;
    bool goingleft, going_up;
public:
    RGTMarginal(
	const double* masses,
        const double* probs,
        int isotopeNo,
        int atomCnt,
        double lCutOff,
        int tabSize = 1000,
        int hashSize = 1000
    );
    ~RGTMarginal();
    void setup_search(double _pmin, double _pmax, double _mmin, double _mmax);
    bool next();
    void terminate_search();
    inline const double& current_lProb() const { return lProbs[cidx]; };
    inline const double& current_mass() const { return masses[cidx]; };
    inline const Conf& current_conf() const { return confs[cidx]; };

private:
    unsigned int* alloc_and_setup_subintervals();
    unsigned int setup_subintervals(unsigned int* T, unsigned int idx, bool left);
    double* alloc_and_setup_mass_table();
    bool hard_next();






};

#endif
