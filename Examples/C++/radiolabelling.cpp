#include <iostream>
#include <cassert>
#include "../../IsoSpec++/fixedEnvelopes.h"

using namespace IsoSpec;

int main()
{
  // Calculates the isotopic distribution of isotopically labelled glucose, where two of the carbons are replaced with 14C
  // We assume that radiolabelling is 95% efficient (that is, there is a 95% chance for each of the radiolabel atoms to be 14C, and
  // a 5% chance of them being either 12C or 13C, with the probability of each of those proportional to their standard isotopic abundances)


  {
    // First, we construct the glucose molecule WITHOUT the radiolabel atoms - that is C4 instead of C6
    Iso i("C4H12O6");

    // Then, we add the special "element", representing the radiolabel, with custom isotopic distribution (of 3 possible isotopes) and 2 atoms
    const double radiocarbon_masses[3] = {12.0, 13.0033548352, 14.003241989};
    const double radiocarbon_probs[3] = {0.05*0.989211941850466, 0.05*0.010788058149533084, 0.95}; // The standard isotopic distribution of nonradio-carbon multiplied by 0.05, and 0.95 of 14C

    i.addElement(2, 3, radiocarbon_masses, radiocarbon_probs);

    // Additional radiolabel elements can be added by more calls to addElement

    FixedEnvelope iso = FixedEnvelope::FromTotalProb(std::move(i), 0.99, true, true);

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << iso.masses()[0] << std::endl;
    std::cout << "probability: " << iso.probs()[0] << std::endl;

    const int* configs = iso.confs();

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "12C atoms: " << configs[0] + configs[7] << std::endl; // Counting the normal and unsuccesfully radiolabelled atoms
    std::cout << "13C atoms: " << configs[1] + configs[8] << std::endl;
    std::cout << "14C atoms: " << configs[9] << std::endl;
    std::cout << "Protium atoms: " << configs[2] << std::endl;
    std::cout << "Deuterium atoms " << configs[3] << std::endl;
    std::cout << "16O atoms: " << configs[4] << std::endl;
    std::cout << "17O atoms: " << configs[5] << std::endl;
    std::cout << "18O atoms: " << configs[6] << std::endl;


    std::cout << "Probabilities of the remaining computed configurations of the distribution are: " << std::endl;

    for(int ii = 1; ii<iso.confs_no(); ii++)
        std::cout << iso.probs()[ii] << std::endl;
  }
}
