#include <iostream>
#include "../../IsoSpec++/isoSpec++.h"

using namespace IsoSpec;

int main()
{
  int configs[5];
  {
    IsoOrderedGenerator iso("H2O1");

    iso.advanceToNextConfiguration();

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << iso.mass() << std::endl;
    std::cout << "log-prob: " << iso.lprob() << std::endl;
    std::cout << "probability: " << iso.prob() << std::endl;

    iso.get_conf_signature(configs);

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "Protium atoms: " << configs[0] << std::endl;
    std::cout << "Deuterium atoms " << configs[1] << std::endl;
    std::cout << "O16 atoms: " << configs[2] << std::endl;
    std::cout << "O17 atoms: " << configs[3] << std::endl;
    std::cout << "O18 atoms: " << configs[4] << std::endl;
  }



  {
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
 
    IsoOrderedGenerator iso(Iso(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs));

    iso.advanceToNextConfiguration();

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << iso.mass() << std::endl;
    std::cout << "log-prob: " << iso.lprob() << std::endl;
    std::cout << "probability: " << iso.prob() << std::endl;

    iso.get_conf_signature(configs);

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "Protium atoms: " << configs[0] << std::endl;
    std::cout << "Deuterium atoms " << configs[1] << std::endl;
    std::cout << "O16 atoms: " << configs[2] << std::endl;
    std::cout << "O17 atoms: " << configs[3] << std::endl;
    std::cout << "O18 atoms: " << configs[4] << std::endl;

    std::cout << "Probabilities of the remaining configurations, in a (guaranteed) nonincreasing order, are: " << std::endl;

    while(iso.advanceToNextConfiguration())
        std::cout << iso.prob() << std::endl;
  }

  // TODO: demonstrate other algorithms
}
