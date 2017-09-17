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


class Marginal
{
private:
    bool disowned;
protected:
    const unsigned int isotopeNo;
    const unsigned int atomCnt;
    const double* const atom_masses;
    const double* const atom_lProbs;
    const Conf mode_conf;
    const double mode_lprob;


public:
    Marginal(
        const double* _masses,   // masses size = logProbs size = isotopeNo
        const double* _probs,
        int _isotopeNo,                  // No of isotope configurations.
        int _atomCnt
    );
    Marginal(Marginal& other) = delete;
    Marginal& operator= (const Marginal& other) = delete;
    Marginal(Marginal&& other);
    virtual ~Marginal();

    inline int get_isotopeNo() const { return isotopeNo; };
    double getLightestConfMass() const;
    double getHeaviestConfMass() const;
    inline double getModeLProb() const { return mode_lprob; };
};

class MarginalTrek : public Marginal
{
    int current_count;
private:
    const KeyHasher keyHasher;
    const ConfEqual equalizer;
    const ConfOrderMarginal orderMarginal;
    std::unordered_map<Conf,int,KeyHasher,ConfEqual> visited;
    std::priority_queue<Conf,std::vector<Conf>,ConfOrderMarginal> pq;
    Summator totalProb;
    Conf candidate;
    Allocator<int> allocator;
    std::vector<double> _conf_probs;
    std::vector<double> _conf_masses;
    std::vector<int*> _confs;


    bool add_next_conf();

public:
    MarginalTrek(
        Marginal&& m,
        int tabSize = 1000,
        int hashSize = 1000
    );

    inline bool probeConfigurationIdx(int idx)
    {
        while(current_count <= idx)
            if(not add_next_conf())
                return false;
        return true;
    }

    int processUntilCutoff(double cutoff);

    inline const std::vector<double>& conf_probs() const { return _conf_probs; };
    inline const std::vector<double>& conf_masses() const { return _conf_masses; };
    inline const std::vector<int*>& confs() const { return _confs; };


    virtual ~MarginalTrek();
};



class PrecalculatedMarginal : public Marginal
{
protected:
    std::vector<Conf> configurations;
    Conf* confs;
    unsigned int no_confs;
    double* masses;
    double* lProbs;
    double* eProbs;
    Allocator<int> allocator;
public:
    PrecalculatedMarginal(
        Marginal&& m,
	double lCutOff,
	bool sort = true,
	int tabSize = 1000,
	int hashSize = 1000
    );
    virtual ~PrecalculatedMarginal();
    inline bool inRange(unsigned int idx) const { return idx < no_confs; };
    inline const double& get_lProb(unsigned int idx) const { return lProbs[idx]; };
    inline const double& get_eProb(unsigned int idx) const { return eProbs[idx]; };
    inline const double& get_mass(unsigned int idx) const { return masses[idx]; };
    inline const double* get_lProbs_ptr() const { return lProbs; };
    inline const double* get_masses_ptr() const { return masses; };
    inline const Conf& get_conf(unsigned int idx) const { return confs[idx]; };
    inline unsigned int get_no_confs() const { return no_confs; };
};

class SyncMarginal : public PrecalculatedMarginal
{
protected:
    char padding[64]; /*  padding[64]; */ // against fake-sharing cache lines...
    std::atomic<unsigned int> counter;
    char padding2[64]; /// likewise...
public:
    inline SyncMarginal(
        Marginal&& m,
        double lCutOff,
        int tabSize = 1000,
        int hashSize = 1000
    ) : PrecalculatedMarginal(
        std::move(m),
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
    unsigned int lower, upper, arridx, arrend, mask, current_level, cidx, gap;
    bool goingleft, going_up;
    ReverseOrder<double> rev_ord;
public:
    RGTMarginal(
        Marginal&& m,
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
    inline const double& current_eProb() const { return eProbs[cidx]; };
    inline const Conf& current_conf() const { return confs[cidx]; };
    double min_mass_above_lProb(double prob);
    double max_mass_above_lProb(double prob);

private:
    unsigned int* alloc_and_setup_subintervals();
    unsigned int setup_subintervals(unsigned int* T, unsigned int idx, bool left);
    double* alloc_and_setup_mass_table();
    bool hard_next();
};

class LayeredMarginal : public Marginal
{
    double current_threshold;
    std::vector<Conf> configurations;
    std::vector<Conf> fringe;
    Allocator<int> allocator;
    unsigned int sorted_up_to_idx;
    const ConfEqual equalizer;
    const KeyHasher keyHasher;
    const ConfOrderMarginalDescending orderMarginal;
    std::vector<double> lProbs;
    std::vector<double> eProbs;
    std::vector<double> masses;
    double* guarded_lProbs;
    const int hashSize;

public:
    LayeredMarginal(Marginal&& m, int tabSize = 1000, int hashSize = 1000);
    bool extend(double new_threshold);
    inline const double& get_lProb(unsigned int idx) const { return guarded_lProbs[idx]; }; // access to idx == -1 is valid and gives a guardian of +inf
    inline const double& get_eProb(unsigned int idx) const { return eProbs[idx]; };
    inline const double& get_mass(unsigned int idx) const { return masses[idx]; };
    inline const Conf& get_conf(unsigned int idx) const { return configurations[idx]; };
    inline unsigned int get_no_confs() const { return configurations.size(); };

};

#endif
