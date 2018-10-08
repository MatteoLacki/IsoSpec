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

#pragma once

#include "isoSpec++.h"

class Spectrum
{
private:
        Iso&& iso;
	double lowest_mass;
	const double bucket_width;
	unsigned long n_buckets;
	double* storage;
        double* ofset_store;
        pthread_t* threads;
        const double cutoff;
        PrecalculatedMarginal** PMs;
        unsigned int n_threads;
        bool absolute;
        std::atomic<unsigned int> thread_idxes;
        double** thread_storages;
        double* thread_partials;
        unsigned int* thread_numbers;
        unsigned int total_confs;
        double total_prob;
        const unsigned long ptr_diff;
        const unsigned long mmap_len;

public:
	Spectrum(Iso&& I, double bucket_width, double cutoff, bool _absolute);
	~Spectrum();
	void add_other(Spectrum& other);
        void run(unsigned int threads = 0, bool sync = true);
        void worker_thread();
        void wait();
        void calc_sum();
	inline unsigned int get_total_confs() const { return total_confs; };
        inline double get_total_prob() const { return total_prob; };
	void print(std::ostream& o = std::cout);

};
