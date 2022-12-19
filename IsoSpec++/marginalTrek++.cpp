/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
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
#include <cstdlib>
#include <queue>
#include <utility>
#include <cstring>
#include <string>
#include <limits>
#include <memory>
#include "platform.h"
#include "marginalTrek++.h"
#include "conf.h"
#include "allocator.h"
#include "operators.h"
#include "summator.h"
#include "element_tables.h"
#include "misc.h"


namespace IsoSpec
{

//! Find one of the most probable subisotopologues.
/*!
    The algorithm uses the hill-climbing algorithm.
    It starts from a subisotopologue close to the mean of the underlying multinomial distribution.
    There might be more than one modes, in case of which this function will return only one of them, close to the mean.

    \param atomCnt

*/
void writeInitialConfiguration(const int atomCnt, const int isotopeNo, const double* lprobs, int* res)
{
    /*!
    Here we perform hill climbing to the mode of the marginal distribution (the subisotopologue distribution).
    We start from the point close to the mean of the underlying multinomial distribution.
    */

    // This approximates the mode (heuristics: the mean is close to the mode).
    for(int i = 0; i < isotopeNo; ++i)
        res[i] = static_cast<int>( atomCnt * exp(lprobs[i]) ) + 1;

    // The number of assigned atoms above.
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
        int i = 0;

        while( diff > 0){
            int coordDiff = res[i] - diff;

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
    double LP = unnormalized_logProb(res, lprobs, isotopeNo);
    double NLP;

    while(modified)
    {
        modified = false;
        for(int ii = 0; ii < isotopeNo; ii++)
            for(int jj = 0; jj < isotopeNo; jj++)
                if(ii != jj && res[ii] > 0)
                {
                    res[ii]--;
                    res[jj]++;
                    NLP = unnormalized_logProb(res, lprobs, isotopeNo);
                    if(NLP > LP || (NLP == LP && ii > jj))
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
}


double* getMLogProbs(const double* probs, int isoNo)
{
    /*!
    Here we order the processor to round the numbers up rather than down.
    Rounding down could result in the algorithm falling in an infinite loop
    because of the numerical instability of summing.
    */
    for(int ii = 0; ii < isoNo; ii++)
        if(probs[ii] <= 0.0 || probs[ii] > 1.0)
            throw std::invalid_argument("All isotope probabilities p must fulfill: 0.0 < p <= 1.0");

    double* ret = new double[isoNo];

    // here we change the table of probabilities and log it.
    for(int i = 0; i < isoNo; i++)
    {
        ret[i] = log(probs[i]);
        for(int j = 0; j < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; j++)
            if(elem_table_probability[j] == probs[i])
            {
                ret[i] = elem_table_log_probability[j];
                break;
            }
    }
    return ret;
}

double get_loggamma_nominator(int x)
{
    // calculate log gamma of the nominator calculated in the binomial exression.
    double ret = lgamma(x+1);
    return ret;
}

int verify_atom_cnt(int atomCnt)
{
    #if !ISOSPEC_BUILDING_OPENMS
    if(ISOSPEC_G_FACT_TABLE_SIZE-1 <= atomCnt)
        throw std::length_error("Subisotopologue too large, size limit (that is, the maximum number of atoms of a single element in a molecule) is: " + std::to_string(ISOSPEC_G_FACT_TABLE_SIZE-1));
    #endif
    return atomCnt;
}

Marginal::Marginal(
    const double* _masses,
    const double* _probs,
    int _isotopeNo,
    int _atomCnt
) :
disowned(false),
isotopeNo(_isotopeNo),
atomCnt(verify_atom_cnt(_atomCnt)),
atom_lProbs(getMLogProbs(_probs, isotopeNo)),
atom_masses(array_copy<double>(_masses, _isotopeNo)),
loggamma_nominator(get_loggamma_nominator(_atomCnt)),
mode_conf(nullptr)
// Deliberately not initializing mode_lprob
{}

Marginal::Marginal(const Marginal& other) :
disowned(false),
isotopeNo(other.isotopeNo),
atomCnt(other.atomCnt),
atom_lProbs(array_copy<double>(other.atom_lProbs, isotopeNo)),
atom_masses(array_copy<double>(other.atom_masses, isotopeNo)),
loggamma_nominator(other.loggamma_nominator)
{
    if(other.mode_conf == nullptr)
    {
        mode_conf = nullptr;
        // Deliberately not initializing mode_lprob. In this state other.mode_lprob is uninitialized too.
    }
    else
    {
        mode_conf = array_copy<int>(other.mode_conf, isotopeNo);
        mode_lprob = other.mode_lprob;
    }
}


// the move-constructor: used in the specialization of the marginal.
Marginal::Marginal(Marginal&& other) :
disowned(other.disowned),
isotopeNo(other.isotopeNo),
atomCnt(other.atomCnt),
atom_lProbs(other.atom_lProbs),
atom_masses(other.atom_masses),
loggamma_nominator(other.loggamma_nominator)
{
    other.disowned = true;
    if(other.mode_conf == nullptr)
    {
        mode_conf = nullptr;
        // Deliberately not initializing mode_lprob. In this state other.mode_lprob is uninitialized too.
    }
    else
    {
        mode_conf = other.mode_conf;
        mode_lprob = other.mode_lprob;
    }
}

Marginal::~Marginal()
{
    if(!disowned)
    {
        delete[] atom_masses;
        delete[] atom_lProbs;
        delete[] mode_conf;
    }
}


Conf Marginal::computeModeConf() const
{
    Conf res = new int[isotopeNo];
    writeInitialConfiguration(atomCnt, isotopeNo, atom_lProbs, res);
    return res;
}

void Marginal::setupMode()
{
    ISOSPEC_IMPOSSIBLE(mode_conf != nullptr);
    mode_conf = computeModeConf();
    mode_lprob = logProb(mode_conf);
}


double Marginal::getLightestConfMass() const
{
    double ret_mass = std::numeric_limits<double>::infinity();
    for(unsigned int ii = 0; ii < isotopeNo; ii++)
        if( ret_mass > atom_masses[ii] )
            ret_mass = atom_masses[ii];
    return ret_mass*atomCnt;
}

double Marginal::getHeaviestConfMass() const
{
    double ret_mass = 0.0;
    for(unsigned int ii = 0; ii < isotopeNo; ii++)
        if( ret_mass < atom_masses[ii] )
            ret_mass = atom_masses[ii];
    return ret_mass*atomCnt;
}

double Marginal::getMonoisotopicConfMass() const
{
    double found_prob = -std::numeric_limits<double>::infinity();
    double found_mass = 0.0;  // to avoid uninitialized var warning
    for(unsigned int ii = 0; ii < isotopeNo; ii++)
        if( found_prob < atom_lProbs[ii] )
        {
            found_prob = atom_lProbs[ii];
            found_mass = atom_masses[ii];
        }
    return found_mass*atomCnt;
}

double Marginal::getAtomAverageMass() const
{
    double ret = 0.0;
    for(unsigned int ii = 0; ii < isotopeNo; ii++)
        ret += exp(atom_lProbs[ii]) * atom_masses[ii];
    return ret;
}

double Marginal::variance() const
{
    double ret = 0.0;
    double mean = getAtomAverageMass();
    for(size_t ii = 0; ii < isotopeNo; ii++)
    {
        double msq = atom_masses[ii] - mean;
        ret += exp(atom_lProbs[ii]) * msq * msq;
    }
    return ret * atomCnt;
}

double Marginal::getLogSizeEstimate(double logEllipsoidRadius) const
{
    if(isotopeNo <= 1)
        return -std::numeric_limits<double>::infinity();

    const double i = static_cast<double>(isotopeNo);
    const double k = i - 1.0;
    const double n = static_cast<double>(atomCnt);

    double sum_lprobs = 0.0;
    for(int jj = 0; jj < i; jj++)
        sum_lprobs += atom_lProbs[jj];

    double log_V_simplex = k * log(n) - lgamma(i);
    double log_N_simplex = lgamma(n+i) - lgamma(n+1.0) - lgamma(i);
    double log_V_ellipsoid = (k * (log(n) + logpi + logEllipsoidRadius) + sum_lprobs) * 0.5 - lgamma((i+1)*0.5);

    return log_N_simplex + log_V_ellipsoid - log_V_simplex;
}


// this is roughly an equivalent of IsoSpec-Threshold-Generator
MarginalTrek::MarginalTrek(
    Marginal&& m,
    int tabSize,
    int
) :
Marginal(std::move(m)),
current_count(0),
orderMarginal(atom_lProbs, isotopeNo),
pq(),
allocator(isotopeNo, tabSize),
min_lprob(*std::min_element(atom_lProbs, atom_lProbs+isotopeNo))
{
    int* initialConf = allocator.makeCopy(mode_conf);

    pq.push({mode_lprob, initialConf});

    current_count = 0;

    fringe.resize_and_wipe(1);

    current_bucket = 0;
    initialized_until = 1;

    add_next_conf();
}


bool MarginalTrek::add_next_conf()
{
    /*!
    Add next configuration.
    If visited all, return false.
    */
    if(pq.empty())
    {
        current_bucket++;
        while(current_bucket < initialized_until && fringe[current_bucket].empty())
        {
//            std::cout << "EMPTY bucket, id: " << current_bucket << std::endl;
            current_bucket++;
        }

//        std::cout << "Entering bucket, size: " << fringe[current_bucket].size() << std::endl;

        if(current_bucket >= initialized_until)
            return false;

    //    std::cout << "Fringe size at pop: " << fringe[current_bucket].size() << std::endl;
        pq = std::priority_queue<ProbAndConfPtr, pod_vector<ProbAndConfPtr> >(std::less<ProbAndConfPtr>(), pod_vector<ProbAndConfPtr>(std::move(fringe[current_bucket])));
    };

    double logprob = pq.top().first;
    Conf topConf = pq.top().second;

    pq.pop();
    ++current_count;

    _confs.push_back(topConf);

    _conf_masses.push_back(calc_mass(topConf, atom_masses, isotopeNo));
    _conf_lprobs.push_back(logprob);

    for( unsigned int j = 0; j < isotopeNo; ++j )
    {
        if( topConf[j] > mode_conf[j])
            continue;

        if( topConf[j] > 0 )
        {
            for( unsigned int i = 0; i < isotopeNo; ++i )
            {
                if( topConf[i] < mode_conf[i] )
                    continue;
                // Growing index different than decreasing one AND
                // Remain on simplex condition.
                if( i != j ){
                    Conf acceptedCandidate = allocator.makeCopy(topConf);

                    ++acceptedCandidate[i];
                    --acceptedCandidate[j];

                    double new_prob = logProb(acceptedCandidate);
                    size_t bucket_nr = bucket_no(new_prob);

                    if(bucket_nr >= initialized_until)
                    {
                    //    std::cout << "Extending to: " << bucket_nr << std::endl;
                        initialized_until = bucket_nr+1;
                        fringe.resize_and_wipe(initialized_until);
                    }

                    ISOSPEC_IMPOSSIBLE(bucket_nr < current_bucket);
                    if(bucket_nr == current_bucket)
                        pq.push({new_prob, acceptedCandidate});
                    else
                        fringe[bucket_nr].push_back({new_prob, acceptedCandidate});

                }

                if( topConf[i] > mode_conf[i] )
                    break;
            }
        }
        if( topConf[j] < mode_conf[j] )
            break;
    }

    return true;
}


MarginalTrek::~MarginalTrek()
{
    const size_t fringe_size = fringe.size();
    for(size_t ii = 0; ii < fringe_size; ii++)
        fringe[ii].clear();
}



PrecalculatedMarginal::PrecalculatedMarginal(Marginal&& m,
    double lCutOff,
    bool sort,
    int tabSize,
    int
) : Marginal(std::move(m)),
allocator(isotopeNo, tabSize)
{
    Conf currentConf = allocator.makeCopy(mode_conf);
    if(logProb(currentConf) >= lCutOff)
    {
        configurations.push_back(currentConf);
        lProbs.push_back(mode_lprob);
    }

    unsigned int idx = 0;

    std::unique_ptr<double[]> prob_partials(new double[isotopeNo]);
    std::unique_ptr<double[]> prob_part_acc(new double[isotopeNo+1]);
    prob_part_acc[0] = loggamma_nominator;

    while(idx < configurations.size())
    {
        currentConf = configurations[idx];
        idx++;

        for(size_t ii = 0; ii < isotopeNo; ii++)
            prob_partials[ii] = minuslogFactorial(currentConf[ii]) + currentConf[ii] * atom_lProbs[ii];

        for(unsigned int ii = 0; ii < isotopeNo; ii++ )
        {
            if(currentConf[ii] > mode_conf[ii])
                continue;

            if(currentConf[ii] != 0)
            {
                double prev_partial_ii = prob_partials[ii];
                currentConf[ii]--;
                prob_partials[ii] = minuslogFactorial(currentConf[ii]) + currentConf[ii] * atom_lProbs[ii];

                for(unsigned int jj = 0; jj < isotopeNo; jj++ )
                {
                    prob_part_acc[jj+1] = prob_part_acc[jj] + prob_partials[jj];

                    if(currentConf[jj] < mode_conf[jj])
                        continue;

                    if( ii != jj )
                    {
                        double logp = prob_part_acc[jj] + minuslogFactorial(1+currentConf[jj]) + (1+currentConf[jj]) * atom_lProbs[jj];
                        for(size_t kk = jj+1; kk < isotopeNo; kk++)
                            logp += prob_partials[kk];

                        if (logp >= lCutOff)
                        {
                            auto tmp = allocator.makeCopy(currentConf);
                            tmp[jj]++;
                            configurations.push_back(tmp);
                            lProbs.push_back(logp);
                        }
                    }
                    else
                        prob_part_acc[jj+1] = prob_part_acc[jj] + prob_partials[jj];

                    if (currentConf[jj] > mode_conf[jj])
                        break;
                }
                currentConf[ii]++;
                prob_partials[ii] = prev_partial_ii;
            }

            if(currentConf[ii] < mode_conf[ii])
                break;
        }
    }

    no_confs = configurations.size();
    confs  = configurations.data();

    if(sort && no_confs > 0)
    {
            std::unique_ptr<size_t[]> order_arr(get_inverse_order(lProbs.data(), no_confs));
            impose_order(order_arr.get(), no_confs, lProbs.data(), confs);
    }

    probs = new double[no_confs];
    masses = new double[no_confs];


    for(unsigned int ii = 0; ii < no_confs; ii++)
    {
        probs[ii] = exp(lProbs[ii]);
        masses[ii] = calc_mass(confs[ii], atom_masses, isotopeNo);
    }

    lProbs.push_back(-std::numeric_limits<double>::infinity());
}


PrecalculatedMarginal::~PrecalculatedMarginal()
{
    if(masses != nullptr)
        delete[] masses;
    if(probs != nullptr)
        delete[] probs;
}







LayeredMarginal::LayeredMarginal(Marginal&& m, int tabSize, int)
: Marginal(std::move(m)), current_threshold(1.0), allocator(isotopeNo, tabSize),
equalizer(isotopeNo), keyHasher(isotopeNo)
{
    fringe.push_back(mode_conf);
    lProbs.push_back(std::numeric_limits<double>::infinity());
    fringe_unn_lprobs.push_back(unnormalized_logProb(mode_conf));
    lProbs.push_back(-std::numeric_limits<double>::infinity());
    guarded_lProbs = lProbs.data()+1;
}

bool LayeredMarginal::extend(double new_threshold, bool do_sort)
{
    new_threshold -= loggamma_nominator;
    if(fringe.empty())
        return false;

    lProbs.pop_back();  // Remove the +inf guardian

    pod_vector<Conf> new_fringe;
    pod_vector<double> new_fringe_unn_lprobs;

    while(!fringe.empty())
    {
        Conf currentConf = fringe.back();
        fringe.pop_back();

        double opc = fringe_unn_lprobs.back();

        fringe_unn_lprobs.pop_back();
        if(opc < new_threshold)
        {
            new_fringe.push_back(currentConf);
            new_fringe_unn_lprobs.push_back(opc);
        }

        else
        {
            configurations.push_back(currentConf);
            lProbs.push_back(opc+loggamma_nominator);
            for(unsigned int ii = 0; ii < isotopeNo; ii++ )
            {
                if(currentConf[ii] > mode_conf[ii])
                    continue;

                if(currentConf[ii] > 0)
                {
                    currentConf[ii]--;
                    for(unsigned int jj = 0; jj < isotopeNo; jj++ )
                    {
                        if(currentConf[jj] < mode_conf[jj])
                            continue;

                        if( ii != jj )
                        {
                            Conf nc = allocator.makeCopy(currentConf);
                            nc[jj]++;

                            double lpc = unnormalized_logProb(nc);
                            if(lpc >= new_threshold)
                            {
                                fringe.push_back(nc);
                                fringe_unn_lprobs.push_back(lpc);
                            }
                            else
                            {
                                new_fringe.push_back(nc);
                                new_fringe_unn_lprobs.push_back(lpc);
                            }
                        }

                        if(currentConf[jj] > mode_conf[jj])
                            break;
                    }
                    currentConf[ii]++;
                }

                if(currentConf[ii] < mode_conf[ii])
                    break;
            }
        }
    }

    current_threshold = new_threshold;
    fringe.swap(new_fringe);
    fringe_unn_lprobs.swap(new_fringe_unn_lprobs);

    if(do_sort)
    {
        size_t to_sort_size = configurations.size() - probs.size();
        if(to_sort_size > 0)
        {
            std::unique_ptr<size_t[]> order_arr(get_inverse_order(lProbs.data()+1+probs.size(), to_sort_size));
            double* P = lProbs.data()+1+probs.size();
            Conf* C = configurations.data()+probs.size();
            size_t* O = order_arr.get();
            impose_order(O, to_sort_size, P, C);
        }
    }

    if(probs.capacity() * 2 < configurations.size() + 2)
    {
        // Reserve space for new values
        probs.reserve(configurations.size());
        masses.reserve(configurations.size());
    }  // Otherwise we're growing slowly enough that standard reallocations on push_back work better - we waste some extra memory
       // but don't reallocate on every call

//    printVector(lProbs);
    for(unsigned int ii = probs.size(); ii < configurations.size(); ii++)
    {
        probs.push_back(exp(lProbs[ii+1]));
        masses.push_back(calc_mass(configurations[ii], atom_masses, isotopeNo));
    }

    lProbs.push_back(-std::numeric_limits<double>::infinity());  // Restore guardian

    guarded_lProbs = lProbs.data()+1;  // Vector might have reallocated its backing storage

    return true;
}


double LayeredMarginal::get_min_mass() const
{
    double ret = std::numeric_limits<double>::infinity();
    for(pod_vector<double>::const_iterator it = masses.cbegin(); it != masses.cend(); ++it)
        if(*it < ret)
            ret = *it;
    return ret;
}


double LayeredMarginal::get_max_mass() const
{
    double ret = -std::numeric_limits<double>::infinity();
    for(pod_vector<double>::const_iterator it = masses.cbegin(); it != masses.cend(); ++it)
        if(*it > ret)
            ret = *it;
    return ret;
}

}  // namespace IsoSpec
