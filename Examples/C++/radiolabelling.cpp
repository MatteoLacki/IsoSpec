#include <iostream>
#include <cassert>
#include "../../IsoSpec++/isoSpec++.h"

using namespace IsoSpec;

int main()
{
  // Calculates the isotopic distribution of isotopically labelled glucose, where one of the carbons is replaced with 14C
  // We assume that 95% of the molecules are radiolabelled.


  {
    // As one of the atoms has a nonstandard isotopic distribution we can't use the Iso(const char*) constructor (from chemical formula)
    // We have to manually grab isotopic masses and probabilities instead, then construct an artificial "element" which is carbon
    // with one added isotope (14C) and shifted isotopic distribution


    const int elementNumber = 4;
    const int isotopeNumbers[4] = {2,2,3,3};

    // Formula of radiolabeled glucose is C5H12O6(14C)1
    const int atomCounts[4] = {5,12,6,1};

    // First, deal with standard elements. Here we define them manually, they can also be retrieved from stuff in element_tables.h
    const double normal_carbon_masses[2] = {12.0, 13.0033548352};
    const double hydrogen_masses[2] = {1.00782503207, 2.0141017778};
    const double oxygen_masses[3] = {15.99491461956, 16.99913170, 17.9991610};

    // The set of masses for the radiolabelled carbon
    const double radiocarbon_masses[3] = {12.0, 13.0033548352, 14.003241989};

    const double* isotope_masses[4] = {normal_carbon_masses, hydrogen_masses, oxygen_masses, radiocarbon_masses};

    // Now, the isotopic distributions. Standard for normal elements...
    const double normal_carbon_probs[2] = {0.9892119418504669, 0.010788058149533084};
    const double hydrogen_probs[2] = {0.9998842901643079, 0.00011570983569203331};
    const double oxygen_probs[3] = {0.997567609729561, 0.00038099847600609594, 0.002051391794432822};

    // We're assuming that the labelling was only 95% efficient, that is only 95%
    // of the molecules have standard C replaced with 14C. Non-replaced molecules have standard
    // isotopic abundance (realtive to each other)
    const double radiocarbon_probs[3] = {0.05*normal_carbon_probs[0], 0.05*normal_carbon_probs[1], 0.95};

    const double* probs[4] = {normal_carbon_probs, hydrogen_probs, oxygen_probs, radiocarbon_probs};
 
    IsoLayeredGenerator iso(Iso(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs), 0.99);

    iso.advanceToNextConfiguration();

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << iso.mass() << std::endl;
    std::cout << "log-prob: " << iso.lprob() << std::endl;
    std::cout << "probability: " << iso.prob() << std::endl;

    int configs[10];
    iso.get_conf_signature(configs);

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "12C atoms: " << configs[0] + configs[7] << std::endl; // Counting the normal and unsuccesfully radiolabelled atoms
    std::cout << "13C atoms: " << configs[1] + configs[8] << std::endl;
    std::cout << "14C atoms: " << configs[9] << std::endl;
    std::cout << "Protium atoms: " << configs[2] << std::endl;
    std::cout << "Deuterium atoms " << configs[3] << std::endl;
    std::cout << "16O atoms: " << configs[4] << std::endl;
    std::cout << "17O atoms: " << configs[5] << std::endl;
    std::cout << "18O atoms: " << configs[6] << std::endl;


    std::cout << "Probabilities of the remaining configurations are: " << std::endl;

    while(iso.advanceToNextConfiguration())
        std::cout << iso.prob() << std::endl;
  }

  // TODO: demonstrate other algorithms
}
