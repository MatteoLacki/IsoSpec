#include <iostream>
#include "../../IsoSpec++/isoSpec++.hpp"


int main()
{
    IsoSpec* iso = IsoSpec::IsoFromFormula("H2O1", 0.9);

    iso->processConfigurationsUntilCutoff();

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << iso->getNoVisitedConfs() << " element(s)" << std::endl;

    std::tuple<double*,double*,int*,int> product = iso->getCurrentProduct();

    double* masses = std::get<0>(product);
    double* logprobs = std::get<1>(product);
    int* configs = std::get<2>(product);

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << masses[0] << std::endl;
    std::cout << "log-prob: " << logprobs[0] << std::endl;
    std::cout << "probability: " << exp(logprobs[0]) << std::endl;

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "Protium atoms: " << configs[0] << std::endl;
    std::cout << "Deuterium atoms " << configs[1] << std::endl;
    std::cout << "O16 atoms: " << configs[2] << std::endl;
    std::cout << "O17 atoms: " << configs[3] << std::endl;
    std::cout << "O18 atoms: " << configs[4] << std::endl;

    delete iso;
    delete[] masses;
    delete[] logprobs;
    delete[] configs;


    std::cout << "Now what if both isotopes of hydrogen were equally probable, while prob. of O16 was 50%, O17 at 30% and O18 at 20%?" << std::endl;

    const int elementNumber = 2;
    const int isotopeNumbers[2] = {2,3};

    const int atomCounts[2] = {2,1};


    const double hydrogen_masses[2] = {1.00782503207, 2.0141017778};
    const double oxygen_masses[3] = {15.99491461956, 16.99913170, 17.9991610};

    const double* isotope_masses[2] = {hydrogen_masses, oxygen_masses};

    const double hydrogen_probs[2] = {0.5, 0.5};
    const double oxygen_probs[3] = {0.5, 0.3, 0.2};

    const double* probs[2] = {hydrogen_probs, oxygen_probs};



    iso = new IsoSpecLayered(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs, 0.9);

    iso->processConfigurationsUntilCutoff();

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << iso->getNoVisitedConfs() << " element(s)" << std::endl;

    product = iso->getCurrentProduct();

    masses = std::get<0>(product);
    logprobs = std::get<1>(product);
    configs = std::get<2>(product);

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << masses[0] << std::endl;
    std::cout << "log-prob: " << logprobs[0] << std::endl;
    std::cout << "probability: " << exp(logprobs[0]) << std::endl;

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "Protium atoms: " << configs[0] << std::endl;
    std::cout << "Deuterium atoms " << configs[1] << std::endl;
    std::cout << "O16 atoms: " << configs[2] << std::endl;
    std::cout << "O17 atoms: " << configs[3] << std::endl;
    std::cout << "O18 atoms: " << configs[4] << std::endl;

    delete iso;
    delete[] masses;
    delete[] logprobs;
    delete[] configs;
}
