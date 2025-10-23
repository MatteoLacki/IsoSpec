#include <iostream>
//#include "isoSpec++.h"
//#include "fixedEnvelopes.h"
//#include "marginalTrek++.h"
#include "unity-build.cpp"
#include <cassert>

using namespace IsoSpec;

#ifndef ISOSPEC_TESTS_SKIP_MAIN

size_t test_stochastic(const char* formula, size_t molecules, double precision, bool print_confs);

int main(int argc, char** argv)
{
    test_stochastic("ddd", 10, 0.3, true);
    return 0;
}
#endif /* ISOSPEC_TESTS_SKIP_MAIN */


size_t test_stochastic(const char* formula, size_t molecules, double precision, bool print_confs)
{
    std::vector<double> masses = {1.0, 10.0, 100.0, 1000.0, 10000.0};
    std::vector<double> probs = {0.09, 0.19, 0.11, 0.21, 0.4};
    std::vector<int> isotope_numbers = {5};
    std::vector<int> atom_counts = {1};
    Iso iso(1, isotope_numbers.data(), atom_counts.data(), masses.data(), probs.data());

    IsoStochasticGeneratorTemplate<IsoLayeredGeneratorTemplate<SingleAtomMarginal<true>>> ISG(std::move(iso), 10000000);
    //IsoLayeredGeneratorTemplate<LoggingMarginal<SingleAtomMarginal<true>>> ISG(std::move(iso));
    //IsoLayeredGeneratorTemplate<LoggingMarginal<LayeredMarginal>> ISG(std::move(iso));
    //IsoOrderedGeneratorTemplate<LoggingMarginal<SingleAtomMarginal<false>>> ISG(std::move(iso));
    //IsoOrderedGeneratorTemplate<LoggingMarginal<MarginalTrek>> ISG(std::move(iso));

    while(ISG.advanceToNextConfiguration())
    {
        auto mass = ISG.mass();
        auto prob = ISG.prob();
        auto lprob = ISG.lprob();
        std::cout << "Configuration, mass: " << mass << ", prob: " << prob << ", lprob: " << lprob << std::endl;
    }
    return 0;
}
