// Edge cases and input validation.  Ports the old element_zero.cpp and
// empty_iso.cpp checks into doctest cases.

#include <stdexcept>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;

TEST_CASE("zero-probability isotope is rejected") {
    int isotopeNumbers[1] = {3};
    int atomCounts[1] = {100};
    double masses[3] = {1.0, 2.0, 3.0};
    double probs[3] = {0.0, 0.4, 0.6};  // first isotope has probability 0
    const double* IM[1] = {masses};
    const double* IP[1] = {probs};

    CHECK_THROWS_AS(
        Iso(1, isotopeNumbers, atomCounts, IM, IP),
        std::invalid_argument);
}

TEST_CASE("raw-array Iso enumerates consistently across two identical generators") {
    // Ported from empty_iso.cpp: build the same molecule twice and confirm the
    // two threshold generators walk it in lockstep with identical mass/lprob.
    const int elementNumber = 2;
    const int isotopeNumbers[2] = {2, 3};
    const int atomCounts[2] = {200, 100};

    const double hydrogen_masses[2] = {1.00782503207, 2.0141017778};
    const double oxygen_masses[3] = {15.99491461956, 16.99913170, 17.9991610};
    const double* isotope_masses[2] = {hydrogen_masses, oxygen_masses};

    const double hydrogen_probs[2] = {0.5, 0.5};
    const double oxygen_probs[3] = {0.5, 0.3, 0.2};
    const double* probs[2] = {hydrogen_probs, oxygen_probs};

    Iso iso1(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs);
    Iso iso2(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs);

    IsoThresholdGenerator t1(std::move(iso1), 0.01, false);
    IsoThresholdGenerator t2(std::move(iso2), 0.01, false);

    std::size_t n = 0;
    while (t1.advanceToNextConfiguration()) {
        REQUIRE(t2.advanceToNextConfiguration());
        CHECK(t1.lprob() == t2.lprob());
        CHECK(t1.mass() == t2.mass());
        ++n;
    }
    CHECK_FALSE(t2.advanceToNextConfiguration());
    CHECK(n > 0);
}

TEST_CASE("count_confs matches the number actually generated") {
    for (const char* f : test_helpers::small_formulas()) {
        INFO("formula=" << f);
        IsoThresholdGenerator g(f, 1e-6, true);
        std::size_t predicted = g.count_confs();
        g.reset();
        std::size_t actual = 0;
        while (g.advanceToNextConfiguration()) ++actual;
        CHECK(predicted == actual);
    }
}
