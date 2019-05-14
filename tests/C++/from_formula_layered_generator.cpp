#include <iostream>
#include "isoSpec++.h"

using namespace IsoSpec;

size_t test_layered_tabulator(const char* formula, double total_prob, bool print_confs = false, double* returned_min = nullptr);

#ifndef ISOSPEC_TESTS_SKIP_MAIN
int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_layered C10000H1000O1000N1000 0.9999" << std::endl;
		std::cout << "...will print the configurations necessary to cover 0.9999 probability of the above molecule" << std::endl;
                std::cout << "./from_formula_layered C10000H1000O1000N1000 0.9999 false" << std::endl;
                std::cout << "will just count them" << std::endl;

		return -1;
	}

        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);

        double min_p;

	size_t no_confs = test_layered_tabulator(argv[1], atof(argv[2]), print_confs, &min_p);

	std::cout << "The number of visited configurations is:" << no_confs << std::endl;
        std::cout << "Minimum logprob visited is (hex): " << std::hexfloat << min_p << " (dec): " << std::defaultfloat << min_p << " prob: " << exp(min_p) << std::endl;
}
#endif /* #ifndef ISOSPEC_TESTS_SKIP_MAIN */

size_t test_layered_tabulator(const char* formula, double total_prob, bool print_confs, double* returned_min)
{
    IsoLayeredGenerator i(formula, 1000, 1000);
    double minlp = 1.0;
    int cnt = 0;
    int* confspace = nullptr;
    double acc_prob = 0.0;
        
    if(print_confs)
        confspace = new int[i.getAllDim()];

    while(acc_prob < total_prob && i.advanceToNextConfiguration())
    {
        acc_prob += i.prob();
        minlp = std::min(minlp, i.lprob());
        cnt += 1;
        if(print_confs)
        {
                i.get_conf_signature(confspace);
                std::cout << "PROB: " << i.prob() << "  \tMASS: " << i.mass() << "\tCONF: ";
                for(int ii=0; ii<i.getAllDim(); ii++)
                    std::cout << confspace[ii] << " ";
                std::cout << std::endl;
        }
    }

    if (returned_min != nullptr)
        *returned_min = minlp;

    if(print_confs)
        delete[] confspace;

    return cnt;
}
