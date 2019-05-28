#include <iostream>
#include <iomanip>
#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_threshold_simple(const char* formula1, const char* formula2, double threshold, double mass_lower, double mass_upper);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_threshold C10000H1000O1000N1000 C10000H1000O1000N1000 0.01 1000.0 1001.0" << std::endl;
		std::cout << "...will the configurations with probability above 0.01 for the sum of the above molecules" << std::endl;
		return -1;
	}

        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);

	size_t no_visited = test_threshold_simple(argv[1], argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]));
	
	std::cout << "The number of visited configurations is:" << no_visited << std::endl;

}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_threshold_simple(const char* formula1, const char* formula2, double threshold, double mass_lower, double mass_upper)
{
        size_t counted = 0;
        double acc_prob = 0.0;

        Iso iso1(formula1);
        double thr1 = exp(iso1.getModeLProb());
        ThresholdFixedEnvelope tfe1(std::move(iso1), threshold/thr1, true);
        tfe1.sort_by_mass();

        Iso iso2(formula2);
        double thr2 = exp(iso2.getModeLProb());
        ThresholdFixedEnvelope tfe2(std::move(iso2), threshold/thr2, true);
        tfe2.sort_by_mass();


        size_t last_ii2 = 0;

        const size_t confs_no1 = tfe1.confs_no();
        const size_t confs_no2 = tfe2.confs_no();

        const double* masses1 = tfe1.masses();
        const double* masses2 = tfe2.masses();

        const double* probs1 = tfe1.probs();
        const double* probs2 = tfe2.probs();

        std::cout << std::setprecision(20) << std::endl;
        std::cout << "Size 1: " << confs_no1 << std::endl;
        std::cout << "Size 2: " << confs_no2 << std::endl;

        std::cout << "Min mass 1: " << masses1[0] << " \t max mass 1: " << masses1[confs_no1-1] << std::endl;
        std::cout << "Min mass 2: " << masses2[0] << " \t max mass 2: " << masses2[confs_no2-1] << std::endl;

        for(ssize_t ii=confs_no1-1; ii>=0; ii--)
        {
            while(last_ii2 < confs_no2 && masses1[ii] + masses2[last_ii2] < mass_lower)
                last_ii2++;


            if(last_ii2 == confs_no2)
                break;

            for(size_t ii2 = last_ii2; ii2 < confs_no2; ii2++)
                if(masses1[ii] + masses2[ii2] > mass_upper)
                    break;
                else
                {
                    double prob = probs1[ii]*probs2[ii2];
                    if(prob > threshold)
                    {
                        counted++;
                        acc_prob += prob;
                    }
                }
        }

        std::cout << "Total prob: " << acc_prob << std::endl;

	return counted;

}
