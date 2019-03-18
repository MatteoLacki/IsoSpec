#include <cassert>
#include <iostream>

using namespace IsoSpec;

size_t test_empty()
{
    size_t ret = 0;

    const int elementNumber = 2;
    const int isotopeNumbers[2] = {2,3};

    const int atomCounts[2] = {200,100};


    const double hydrogen_masses[2] = {1.00782503207, 2.0141017778};
    const double oxygen_masses[3] = {15.99491461956, 16.99913170, 17.9991610};

    const double* isotope_masses[2] = {hydrogen_masses, oxygen_masses};

    const double hydrogen_probs[2] = {0.5, 0.5};
    const double oxygen_probs[3] = {0.5, 0.3, 0.2};

    const double* probs[2] = {hydrogen_probs, oxygen_probs};

    Iso iso1(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs);

    IsoThresholdGenerator t1(std::move(iso1), 0.01, false);

    Iso iso2;

    iso2.addElement(atomCounts[0], isotopeNumbers[0], hydrogen_masses, hydrogen_probs);
    iso2.addElement(atomCounts[1], isotopeNumbers[1], oxygen_masses, oxygen_probs);

    IsoThresholdGenerator t2(std::move(iso1), 0.01, false);

    while(t1.advanceToNextConfiguration())
    {
        assert(t2.advanceToNextConfiguration());
        assert(t1.lprob() == t2.lprob());
        assert(t1.mass() == t2.mass());
        ret++;
    }
    assert(!t2.advanceToNextConfiguration());

    assert(ret > 0);

    return ret;

}

void test_empty_and_print()
{
    #if defined(ISOSPEC_TESTS_MEMSAN)
    test_empty();
    #else
    size_t no = test_empty();
    std::cout << "Test empty OK, " << no << " confs tested." << std::endl;
    #endif

}

#if !defined(ISOSPEC_TESTS_SKIP_MAIN)
#include "../../IsoSpec++/unity-build.cpp"
int main() { test_empty_and_print(); };
#endif
