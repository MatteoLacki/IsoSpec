#include <iostream>
#include <cassert>
#include "../../IsoSpec++/isoSpec++.h"
#include "../../IsoSpec++/fixedEnvelopes.h"

using namespace IsoSpec;

int main()
{
  {
    // We shall walk through a set of configurations which covers at least 99.9% of the total
    // probability. For water we could obviously go through the entire spectrum (100%), but for 
    // larger molecules the entire spectrum has far too many configurations. Luckily, most of the 
    // likelihood is concentrated in a relatively small set of most probable isotopologues
    // - and this method allows one to quickly calculate such a set, parametrising on the
    // percentage of coverage of total probability space required.
    //
    // This is usually better than just calculating all isotopes with probabilities above a 
    // given threshold, as it allows one to directly parametrise on the accuracy of the 
    // simplified spectrum - that is, the L1 distance between the full and simplified spectrum
    // will be less than (in this example) 0.001.
    //
    // If for some reason one would wish to just calculate a set of peaks with probabilities 
    // above a given threshold - it is possible using the IsoThresholdGenerator class.
    //
    // Note: the returned set will usually contain a bit more configurations than necessary
    // to achieve the desired coverage. These configurations need to be computed anyway, however
    // it is possible to discard them using the optional trim argument.
    //
    // Note that TotalProbFixedEnvelope is just a convenience class - it computes all the necessary
    // configurations and stores them into several arrays, allowing easy access. It's reasonably
    // efficient, but for raw pedal-to-the-metal speed you should use the IsoGenerator classes
    // directly - they are slightly more efficient, and do not store the entire set in memory.
    // This allows you to get the data directly into your own structures without needless copying,
    // and allows one to implement a generator-consumer pattern, in which the spectrum is used
    // (for example: binned with some resolution) on the fly, without ever being stored in memory - 
    // which might matter for large spectra.

    FixedEnvelope iso = FixedEnvelope::FromTotalProb("H2O1", 0.999, true, true);

    for(int ii = 0; ii < iso.confs_no(); ii++)
    {
        std::cout << "Visiting configuration: " << std::endl;
        std::cout << "Mass: " << iso.mass(ii) << std::endl;
        std::cout << "probability: " << iso.prob(ii) << std::endl;

        const int* conf = iso.conf(ii);
        // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
        std::cout << "Protium atoms: " << conf[0] << std::endl;
        std::cout << "Deuterium atoms " << conf[1] << std::endl;
        std::cout << "O16 atoms: " << conf[2] << std::endl;
        std::cout << "O17 atoms: " << conf[3] << std::endl;
        std::cout << "O18 atoms: " << conf[4] << std::endl;

        ii++;
    }
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
 
    FixedEnvelope iso = FixedEnvelope::FromTotalProb(Iso(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs), 0.99, true, true);

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << iso.mass(0) << std::endl;
    std::cout << "probability: " << iso.prob(0) << std::endl;

    const int* conf = iso.conf(0);

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "Protium atoms: " << conf[0] << std::endl;
    std::cout << "Deuterium atoms " << conf[1] << std::endl;
    std::cout << "O16 atoms: " << conf[2] << std::endl;
    std::cout << "O17 atoms: " << conf[3] << std::endl;
    std::cout << "O18 atoms: " << conf[4] << std::endl;

    std::cout << "Probabilities of the remaining configurations are: " << std::endl;

    for(int ii=0; ii<iso.confs_no(); ii++)
        std::cout << iso.prob(ii) << std::endl;
  }

  // TODO: demonstrate other algorithms
}
