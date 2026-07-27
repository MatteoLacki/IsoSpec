// The SingleAtomMarginal generator specializations.
//
// isoSpec++.cpp explicitly instantiates the layered and ordered generators (and
// the stochastic wrappers around them) over SingleAtomMarginal<>, a marginal
// for elements present as a single atom, where every "subisotopologue" is just
// one isotope.  Nothing in the public API reaches these instantiations, so they
// would otherwise go completely untested; this file covers them.

#include <algorithm>
#include <cmath>
#include <random>
#include <set>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// One fictitious element with five well-separated isotopes, one atom of it.
const double kMasses[5] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
const double kProbs[5] = {0.09, 0.19, 0.11, 0.21, 0.4};

Iso single_atom_iso(int elements = 1) {
    std::vector<int> isotope_numbers(elements, 5);
    std::vector<int> atom_counts(elements, 1);
    std::vector<double> masses, probs;
    for (int i = 0; i < elements; ++i) {
        masses.insert(masses.end(), kMasses, kMasses + 5);
        probs.insert(probs.end(), kProbs, kProbs + 5);
    }
    return Iso(elements, isotope_numbers.data(), atom_counts.data(),
               masses.data(), probs.data());
}

}  // namespace

TEST_CASE("layered generator over SingleAtomMarginal enumerates the isotopes") {
    IsoLayeredGeneratorTemplate<SingleAtomMarginal<true>> g(single_atom_iso());

    std::vector<Peak> got;
    double last = 2.0;
    // Heap-allocated (and oversized): get_conf_by_indexes writes one int per
    // element, and a stack array of exactly that size makes GCC's -Warray-bounds
    // complain about the loop bound it cannot prove.
    std::vector<int> index_space(8);
    std::set<int> indexes;
    while (g.advanceToNextConfiguration()) {
        // Layers descend, so peaks come out in non-increasing probability.
        REQUIRE(g.prob() <= last + 1e-15);
        last = g.prob();
        g.get_conf_by_indexes(index_space.data());
        // The index is the isotope's position in the input arrays.
        REQUIRE(index_space[0] >= 0);
        REQUIRE(index_space[0] < 5);
        REQUIRE(g.mass() == doctest::Approx(kMasses[index_space[0]]));
        REQUIRE(g.prob() == doctest::Approx(kProbs[index_space[0]]));
        indexes.insert(index_space[0]);
        got.push_back({g.mass(), g.prob()});
    }

    CHECK(got.size() == 5);
    CHECK(indexes.size() == 5);  // a bijection onto the isotopes
    CHECK(total_prob(got) == doctest::Approx(1.0));
}

TEST_CASE("ordered generator over SingleAtomMarginal enumerates the isotopes") {
    IsoOrderedGeneratorTemplate<SingleAtomMarginal<false>> g(single_atom_iso());

    std::vector<Peak> got;
    double last = 2.0;
    while (g.advanceToNextConfiguration()) {
        REQUIRE(g.prob() <= last + 1e-15);
        last = g.prob();
        got.push_back({g.mass(), g.prob()});
    }

    REQUIRE(got.size() == 5);
    CHECK(total_prob(got) == doctest::Approx(1.0));

    std::vector<Peak> expected;
    for (int i = 0; i < 5; ++i) expected.push_back({kMasses[i], kProbs[i]});
    CHECK(peaks_close(got, expected, 1e-12));
}

TEST_CASE("DOC: ordered generator's get_conf_by_indexes is not a bijection") {
    // IsoOrderedGeneratorTemplate::get_conf_by_indexes() computes
    // `space[0] = max(c[0]-1, 0)`, which collapses the first two positions:
    // for a five-isotope single-atom molecule it reports 0,1,2,3,3 — the last
    // isotope is unreachable and index 3 is reported twice.  The layered
    // specialization, which goes through get_original_position(), is correct
    // (checked above).  Nothing in the public API reaches this path, so it is
    // documented rather than asserted.
    IsoOrderedGeneratorTemplate<SingleAtomMarginal<false>> g(single_atom_iso());

    std::multiset<int> reported;
    std::vector<int> index_space(8);
    while (g.advanceToNextConfiguration()) {
        g.get_conf_by_indexes(index_space.data());
        reported.insert(index_space[0]);
    }
    CHECK(reported.size() == 5);
    CHECK(std::set<int>(reported.begin(), reported.end()).size() < reported.size());

    MESSAGE("IsoOrderedGeneratorTemplate::get_conf_by_indexes reports duplicate "
            "indices for SingleAtomMarginal molecules (max(c[0]-1, 0)).");
}

TEST_CASE("stochastic sampling over the single-atom specializations") {
    const std::size_t molecules = 100000;

    SUBCASE("over the layered generator") {
        std::mt19937 rng(20260727);
        IsoStochasticGeneratorTemplate<IsoLayeredGeneratorTemplate<SingleAtomMarginal<true>>>
            g(single_atom_iso(), molecules, 0.9999, 5.0, rng);

        double total = 0.0;
        std::size_t peaks = 0;
        while (g.advanceToNextConfiguration()) {
            REQUIRE(g.prob() >= 1.0);
            REQUIRE(g.prob() == doctest::Approx(std::round(g.prob())));
            total += g.prob();
            ++peaks;
        }
        CHECK(peaks > 0);
        CHECK(peaks <= 5);
        CHECK(total == doctest::Approx(static_cast<double>(molecules)));
    }

    SUBCASE("over the ordered generator") {
        std::mt19937 rng(20260727);
        IsoStochasticGeneratorTemplate<IsoOrderedGeneratorTemplate<SingleAtomMarginal<false>>>
            g(single_atom_iso(), molecules, 0.9999, 5.0, rng);

        double total = 0.0;
        std::size_t peaks = 0;
        while (g.advanceToNextConfiguration()) {
            total += g.prob();
            ++peaks;
        }
        CHECK(peaks > 0);
        CHECK(total == doctest::Approx(static_cast<double>(molecules)));
    }
}

TEST_CASE("single-atom specialization with several elements") {
    // Two independent single-atom elements: 25 configurations, and the
    // distribution is the product of the two marginals.
    IsoLayeredGeneratorTemplate<SingleAtomMarginal<true>> g(single_atom_iso(2));

    std::vector<Peak> got;
    while (g.advanceToNextConfiguration()) got.push_back({g.mass(), g.prob()});

    CHECK(got.size() == 25);
    CHECK(total_prob(got) == doctest::Approx(1.0));

    std::vector<Peak> expected;
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            expected.push_back({kMasses[i] + kMasses[j], kProbs[i] * kProbs[j]});
    CHECK(peaks_close(merge_equal_masses(got), merge_equal_masses(expected), 1e-12));
}

TEST_CASE("single-atom layered generator against the general one") {
    // The specialization must produce the same distribution as the general
    // LayeredMarginal path for the same molecule.
    std::vector<Peak> special;
    {
        IsoLayeredGeneratorTemplate<SingleAtomMarginal<true>> g(single_atom_iso(2));
        while (g.advanceToNextConfiguration()) special.push_back({g.mass(), g.prob()});
    }

    std::vector<Peak> general;
    {
        IsoLayeredGenerator g(single_atom_iso(2));
        while (g.advanceToNextConfiguration()) general.push_back({g.mass(), g.prob()});
    }

    CHECK(special.size() == general.size());
    CHECK(peaks_close(merge_equal_masses(special), merge_equal_masses(general), 1e-12));
}
