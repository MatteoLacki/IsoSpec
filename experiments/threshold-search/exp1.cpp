#include <iostream>
#include "../../IsoSpec++/unity-build.cpp"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_threshold_simple(const char* formula, int count, bool print_confs);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Proper usage (for example): ./from_formula_threshold C10000H1000O1000N1000 0.01" << std::endl;
		std::cout << "...will the configurations with probability above 0.01 for the above molecule" << std::endl;
		return -1;
	}

        bool print_confs = true;

        if(argc > 3)
            print_confs = (strcmp(argv[3], "true") == 0);

	size_t no_visited = test_threshold_simple(argv[1], atof(argv[2]), print_confs);
	
	std::cout << "The number of visited configurations is:" << no_visited << std::endl;

}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_threshold_simple(const char* formula, int count, bool print_confs)
{

    std::cout << "selection" << std::endl;
    Iso molecule(formula);
    double threshold = exp(molecule.getModeLProb());

    double thr_ratio = 1.25;//64.0;
    int last_confs = 0;
    int confs = 1;
    while(confs < count)
    {
        confs = 1;
        threshold = threshold / thr_ratio;
        last_confs = confs;
        confs = IsoThresholdGenerator(Iso(molecule, false), threshold, true).count_confs();
        std::cout << confs << " " << threshold << std::endl;
    }

    std::cout << confs << " " << count << " ratio: " << ((double)confs)/count << std::endl;
    std::cout << "binsearch" << std::endl;

    IsoThresholdGenerator ITG(formula, threshold);
    double left_thr = threshold * thr_ratio;
    while(confs != count)
    {
        double mid_thr = (left_thr + threshold)/2.0;
        ITG.new_threshold(mid_thr);
        ITG.reset();
        int mid_count = ITG.count_confs();
        std::cout << mid_count << " " << " " << count << std::endl;
        if(mid_count < count)
            left_thr = mid_thr;
        else
            threshold = mid_thr;
        confs = mid_count;
    }
    std::cout << "accumulation..." << std::endl;
    IsoThresholdGenerator i(formula, threshold, true, 10, 10, false);
    double total_prob = 0.0;
    double* probs = new double[confs];
    double* masses = new double[confs];

    int idx = 0;
    while(i.advanceToNextConfiguration())
    {
        probs[idx] = i.prob();
        masses[idx] = i.mass();
    }
    std::cout << "selection..." << std::endl;
    std::nth_element(probs, probs+count, probs+confs);
//        total_prob += i.prob();
//    std::cout << "Done! total prob: " << total_prob << std::endl;
    return confs;
}
