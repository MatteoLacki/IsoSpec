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
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <utility>
#include <iostream>
#include <string.h>
#include <limits>
#include "lang.h"
#include "marginalTrek++.h"
#include "conf.h"
#include "allocator.h"
#include "operators.h"
#include "summator.h"
#include "element_tables.h"
#include "misc.h"




Conf initialConfigure(const int atomCnt, const int isotopeNo, const double* probs, const double* lprobs)
{
    Conf res = new int[isotopeNo];

    for(int i = 0; i < isotopeNo; ++i )
    {
        res[i] = int( atomCnt * probs[i] ) + 1;
    }

    int s = 0;

    for(int i = 0; i < isotopeNo; ++i) s += res[i];

    int diff = atomCnt - s;

    // Too little: enlarging fist index.
    if( diff > 0 ){
        res[0] += diff;
    }
    // Too much: trying to redistribute the difference: hopefully the first element is the largest.
    if( diff < 0 ){
        diff = abs(diff);
        int i = 0, coordDiff = 0;

        while( diff > 0){
            coordDiff = res[i] - diff;

            if( coordDiff >= 0 ){
                res[i] -= diff;
                diff = 0;
            } else {
                res[i] = 0;
                i++;
                diff = abs(coordDiff);
            }
        }
    }

    // What we computed so far will be very close to the mode: hillclimb the rest of the way
    
    bool modified = true;
    double LP = logProb(res, lprobs, isotopeNo);
    double NLP;

    while(modified)
    {
    	modified = false;
    	for(int ii = 0; ii<isotopeNo; ii++)
	    for(int jj = 0; jj<isotopeNo; jj++)
	        if(ii != jj and res[ii] > 0)
		{
		    res[ii]--;
		    res[jj]++;
		    NLP = logProb(res, lprobs, isotopeNo);
		    if(NLP>LP)
		    {
		    	modified = true;
			LP = NLP;
		    }
		    else
		    {
		    	res[ii]++;
			res[jj]--;
		    }
		}

	    
    }
    return res;
}


#ifndef BUILDING_R
void printMarginal( const std::tuple<double*,double*,int*,int>& results, int dim)
{
    for(int i=0; i<std::get<3>(results); i++){

        std::cout << "Mass = "  << std::get<0>(results)[i] <<
        " log-prob =\t"                 << std::get<1>(results)[i] <<
        " prob =\t"                     << exp(std::get<1>(results)[i]) <<
        "\tand configuration =\t";

        for(int j=0; j<dim; j++) std::cout << std::get<2>(results)[i*dim + j] << " ";

        std::cout << std::endl;
    }
}
#endif


double* getMLogProbs(const double* probs, int isoNo)
{
    double* ret = new double[isoNo];
    for(int i = 0; i < isoNo; i++)
    {
        ret[i] = log(probs[i]);
        for(int j=0; j<NUMBER_OF_ISOTOPIC_ENTRIES; j++)
            if(elem_table_probability[j] == probs[i])
            {
                ret[i] = elem_table_log_probability[j];
                break;
            }
    }
    return ret;
}


MarginalTrek::MarginalTrek(
    const double* masses,   // masses size = logProbs size = isotopeNo
    const double* probs,
    int isotopeNo,                  // No of isotope configurations.
    int atomCnt,
    int tabSize,
    int hashSize
) : current_count(0),
_tabSize(tabSize),
_hashSize(hashSize),
_isotopeNo(isotopeNo),
_atomCnt(atomCnt),
iso_masses(array_copy<double>(masses, isotopeNo)),
logProbs(getMLogProbs(probs, isotopeNo)),
allocator(isotopeNo, _tabSize),
keyHasher(isotopeNo),
equalizer(isotopeNo),
orderMarginal(logProbs, isotopeNo),
visited(_hashSize,keyHasher,equalizer),
pq(orderMarginal),
totalProb(),
candidate(new int[isotopeNo])
{


    int*    initialConf_tmp = initialConfigure(atomCnt, isotopeNo, probs, logProbs);
    int*    initialConf = allocator.makeCopy(initialConf_tmp);
    delete [] initialConf_tmp;


    pq.push(initialConf);
    visited[initialConf] = 0;

    totalProb = Summator();

    current_count = 0;

    add_next_conf();
}


bool MarginalTrek::add_next_conf()
{
    if(pq.size() < 1) return false;

    Conf topConf = pq.top();
    pq.pop();
    ++current_count;
    visited[topConf] = current_count;

    _confs.push_back(topConf);
    _conf_masses.push_back(mass(topConf, iso_masses, _isotopeNo));
    double logprob = logProb(topConf, logProbs, _isotopeNo);
    _conf_probs.push_back(logprob);


    totalProb.add( exp( logprob ) );

    for( unsigned int i = 0; i < _isotopeNo; ++i )
    {
        for( unsigned int j = 0; j < _isotopeNo; ++j )
        {
            // Growing index different than decreasing one AND
            // Remain on simplex condition.
            if( i != j && topConf[j] > 0 ){
                copyConf(topConf, candidate, _isotopeNo);

                ++candidate[i];
                --candidate[j];

                // candidate should not have been already visited.
                if( visited.count( candidate ) == 0 )
                {
                    Conf acceptedCandidate = allocator.makeCopy(candidate);
                    pq.push(acceptedCandidate);

                    visited[acceptedCandidate] = 0;
                }
            }
        }
    }

    return true;
}

int MarginalTrek::processUntilCutoff(double cutoff)
{
    Summator s;
    int last_idx = -1;
    for(unsigned int i=0; i<_conf_probs.size(); i++)
    {
        s.add(_conf_probs[i]);
        if(s.get() >= cutoff)
        {
            last_idx = i;
            break;
        }
    }
    if(last_idx > -1)
        return last_idx;

    while(totalProb.get() < cutoff && add_next_conf()) {}
    return _conf_probs.size();
}

double MarginalTrek::getLightestConfMass()
{
    double ret_mass = std::numeric_limits<double>::infinity();
    for(unsigned int ii=0; ii < _isotopeNo; ii++)
        if( ret_mass > iso_masses[ii] )
	    ret_mass = iso_masses[ii];
    return ret_mass*_atomCnt;
}

double MarginalTrek::getHeaviestConfMass()
{
    double ret_mass = 0.0;
    for(unsigned int ii=0; ii < _isotopeNo; ii++)
        if( ret_mass < iso_masses[ii] )
            ret_mass = iso_masses[ii];
    return ret_mass*_atomCnt;
}

MarginalTrek::~MarginalTrek()
{
    delete[] logProbs;
    delete[] candidate;
    delete[] iso_masses;
}




PrecalculatedMarginal::PrecalculatedMarginal(
        const double* _masses,
        const double* _lprobs,
        int _isotopeNo,
        int _atomCnt,
	double lCutOff,
	bool sort,
        int tabSize,
        int hashSize
) : 
isotopeNo(_isotopeNo),
atomCnt(_atomCnt),
isoMasses(array_copy<double>(_masses, isotopeNo)),
isoLProbs(getMLogProbs(_lprobs, isotopeNo)),
allocator(isotopeNo,tabSize)
{
    const ConfEqual equalizer(isotopeNo);
    const KeyHasher keyHasher(isotopeNo);
    const ConfOrderMarginal orderMarginal(isoLProbs, isotopeNo);

    std::unordered_set<Conf,KeyHasher,ConfEqual> visited(hashSize,keyHasher,equalizer);


    Conf currentConf = initialConfigure(atomCnt, isotopeNo, _lprobs, isoLProbs);
    Conf tmpConf = currentConf;
    if(logProb(currentConf, isoLProbs, _isotopeNo) >= lCutOff)
    {
        configurations.push_back(allocator.makeCopy(currentConf));
        visited.insert(currentConf);
    }

    unsigned int idx = 0;

    while(idx < configurations.size())
    {
        memcpy(currentConf, configurations[idx], sizeof(int)*isotopeNo);
	idx++;
        for( int ii = 0; ii < _isotopeNo; ii++ )
            for( int jj = 0; jj < _isotopeNo; jj++ )
                if( ii != jj and currentConf[jj] > 0)
		{
		    currentConf[ii]++;
		    currentConf[jj]--;

		    if (visited.count(currentConf) == 0 and logProb(currentConf, isoLProbs, _isotopeNo) >= lCutOff)
		    {
		    	 visited.insert(currentConf);
                         configurations.push_back(allocator.makeCopy(currentConf));
	            }

		    currentConf[ii]--;
		    currentConf[jj]++;

                }
    }

    delete[] tmpConf;

    if(sort)
        std::sort(configurations.begin(), configurations.end(), orderMarginal);


    confs  = &configurations[0];
    no_confs = configurations.size();
    lProbs = new double[no_confs];
    masses = new double[no_confs];


    for(unsigned int ii=0; ii < no_confs; ii++)
    {
        lProbs[ii] = logProb(confs[ii], isoLProbs, _isotopeNo);
	masses[ii] = mass(confs[ii], isoMasses, _isotopeNo);
    }

}


PrecalculatedMarginal::~PrecalculatedMarginal()
{
    if(lProbs != nullptr)
    	delete[] lProbs;
    if(masses != nullptr)
    	delete[] masses;
    delete[] isoMasses;
    delete[] isoLProbs;
}



RGTMarginal::RGTMarginal(
	const double* masses,
        const double* probs,
        int isotopeNo,
        int atomCnt,
        double lCutOff,
        int tabSize,
        int hashSize) : 
	PrecalculatedMarginal(masses, probs, isotopeNo, atomCnt, lCutOff, true, tabSize, hashSize),
	tree_overhead(next_pow2(no_confs)-1),
	tree_size(tree_overhead+no_confs),
//	probs_tree(RGTMarginal::alloc_and_construct_ptree()),
	subtree_sizes(alloc_and_setup_subtree_sizes()),
	mass_layer_size(compute_total_no_masses()),
	mass_order(isoMasses),
	subtree_locations(alloc_and_setup_subtree_locations()),
	subintervals(alloc_and_setup_subintervals())
{}
/*
double* RGTMarginal::alloc_and_construct_ptree()
{
	double* ret = new double[tree_size];
	construct_ptree(ret, 0, lProbs, no_confs);
	return ret;
}

void RGTMarginal::construct_ptree(double* new_tree, unsigned int tidx, double* pstart, unsigned int howmany)
{
	if(howmany == 1)
	{
		new_tree[tidx] = *pstart;
		return;
	}
	unsigned int pp = prev_pow2(howmany);
	new_tree[tidx] = pstart[pp];
	construct_ptree(new_tree, tidx*2+1, pstart, pp);
	construct_ptree(new_tree, tidx*2+2, pstart+pp, howmany-pp);
}
*/
unsigned int* RGTMarginal::alloc_and_setup_subtree_sizes()
{
        std::cout << "tree_overhead " << tree_overhead << "\t confs_no: " << no_confs << std::endl;
	unsigned int* ret = new unsigned int[tree_overhead];
        ret--; // 1-indexed array has easier (and faster) child <-> parent idx computations
	setup_subtree_sizes(ret, 1);
	return ret;
}

unsigned int RGTMarginal::setup_subtree_sizes(unsigned int* T, unsigned int idx)
{
    if(idx > tree_overhead)
	return 1;

    idx *= 2;
    T[idx] = setup_subtree_sizes(T, idx) + setup_subtree_sizes(T, idx+1);
    return T[idx];
}

unsigned int RGTMarginal::compute_total_no_masses()
{
	// FIXME: unnecessary pass. Next function could just retrun two values...
	unsigned int ret = 0;
	for(unsigned int ii=1; ii<=tree_size; ii++)
	    ret += subtree_sizes[ii];
	return ret;
}

unsigned int* RGTMarginal::alloc_and_setup_subtree_locations()
{
	// FIXME change to double*
	unsigned int* ret = new unsigned int[tree_overhead];
        ret--;
	setup_subtree_locations(ret, 1, 0);
	return ret;
}

unsigned int RGTMarginal::setup_subtree_locations(unsigned int* T, unsigned int idx, unsigned int csum)
{
	if(idx > tree_overhead)
	    return csum+1;

        
	T[idx] = setup_subtree_locations(T, 2*idx, csum);
	return setup_subtree_locations(T, 2*idx+1, T[idx]+subtree_sizes[idx]);
}

unsigned int* RGTMarginal::alloc_and_setup_subintervals()
{
	unsigned int* ret = new unsigned int[mass_layer_size+2];
	setup_subintervals(ret, 1, true);
	return ret;
}


unsigned int RGTMarginal::setup_subintervals(unsigned int* T, unsigned int idx, bool left)
{
        if(idx > tree_overhead)
        {
                T[idx] = idx - tree_overhead;
                return T[idx];
        }
	unsigned int ileft  = setup_subintervals(T, idx*2, true);
	unsigned int iright = setup_subintervals(T, idx*2+1, false);
	for(unsigned int ii = ileft; ii <= iright; ii++)
	{
		T[idx] = ii;
		idx++;
	}

	std::sort(T, T+mass_layer_size, mass_order);

	if(left)
		return ileft;
	else
		return iright;

}

double* RGTMarginal::alloc_and_setup_mass_table()
{
	// Trading off memory for speed here... One less pointer dereference per conf, and better cache coherency.
	double* ret = new double[mass_layer_size];
	for(unsigned int ii=0; ii<mass_layer_size; ii++)
		ret[ii] = masses[subintervals[ii]];
	return ret;
}

/*
 * -------- end of constructor -------------
 */


#if 0
void RGTMarginal::setup_search(double _pmin, double _pmax, double _mmin, double _mmax)
{
	splitidx = 0;
	pmin = _pmin;
	pmax = _pmax;
	mmin = _mmin;
	mmax = _mmax;
        /*
	while (splitidx < tree_size)
	if(probs_tree[splitidx]<pmin)
		    splitidx = splitidx*2+1;
		else
		    if(probs_tree[splitidx] > pmax)
		        splitidx = splitidx*2+2;
		    else
		        break;

	goingleft = true;
        currentIdx = splitidx;
        */
        if(isinf(pmin) and signbit(pmin))
            // Negative infinity
            lower = 0;
        else
        {
            lower = std::lower_bound(lProbs, lProbs+no_confs, pmin);
            if(lower == no_confs)
            {
                terminate();
                return;
            }
        }

        if (pmax >= 0.0)
            upper = no_confs - 1;
        else
        {
            upper = std::upper_bound(lProbs, lProbs+no_confs, pmax) - 1;

            if(upper == -1)
            {
                terminate();
                return;
            }
        }



        unsigned int crnt = mass_layer_size;
        if(mmin <= masses[lower] and masses[lower] <= mmax)
        {
            crnt++;
            subintervals[crnt] = lower;
        }
        if(mmin <= masses[upper] and masses[upper] <= mmax)
        {
            crnt++;
            subintervals[crnt] = upper;
        }

        arridx = mass_layer_size - 1;
        arrend = crnt;

        break_idx = next_pow2(tree_overhead) - tree_overhead;
}


bool RGTMarginal::next()
{
    // TODO: move to .h and inline this.
    arridx++;
    if(arridx < arrend and mass_table[arridx] <= mmax)
    {
        cidx = subintervals[arridx];
        return true;
    }

    return hard_next();
}

bool RGTMarginal::hard_next()
{
    if(lower == upper)
    {
        if(lower < tree_overhead)
            return false;
        cidx = lower;
        lower = -1;
        upper = -1;
        return true;
    }

    if(upper > lower)
    {
        if(upper >= tree_overhead)
        {
            // We're outside the tree, just return last element, and make sure we go into hard case next time too.
            cidx = upper;
            arridx = 10;
            arrend = 0;
            upper = (upper-1)/2;
            if masses[cidx] 
        }


        }
        if (upper % 2 == 1)
        {
            upper = upper/2;
            // left child
            return hard_next();
        }
        else
        {
            upper = (upper-2)/2;
            arridx = 
        

    }
}

#endif /* 0 */


