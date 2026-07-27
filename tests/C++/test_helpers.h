// Shared helpers for the doctest-based C++ test suite.
//
// Header-only: every test translation unit includes this; the definitions are
// inline so they can be pulled into multiple TUs without ODR violations.  The
// IsoSpec symbols themselves come from unity-build.cpp, linked into the test
// binary (see the Makefile).

#ifndef ISOSPEC_TEST_HELPERS_H
#define ISOSPEC_TEST_HELPERS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "isoSpec++.h"
#include "fixedEnvelopes.h"

namespace test_helpers {

using namespace IsoSpec;

// A single materialized isotopologue peak.
struct Peak {
    double mass;
    double prob;
};

// Order peaks canonically (by mass, then prob) so two enumerations of the same
// distribution — produced in different orders by different generators — become
// element-wise comparable.
inline void sort_peaks(std::vector<Peak>& v) {
    std::sort(v.begin(), v.end(), [](const Peak& a, const Peak& b) {
        if (a.mass != b.mass) return a.mass < b.mass;
        return a.prob < b.prob;
    });
}

// Full enumeration (all configurations) via the threshold generator with a zero
// cutoff.  For molecules built from real elements every configuration has
// nonzero probability, so this is the entire support.
inline std::vector<Peak> enumerate_threshold_full(const char* formula) {
    std::vector<Peak> out;
    IsoThresholdGenerator g(formula, 0.0, true);
    while (g.advanceToNextConfiguration())
        out.push_back({g.mass(), g.prob()});
    return out;
}

// Full enumeration via the ordered generator (descending probability) run to
// exhaustion.
inline std::vector<Peak> enumerate_ordered_full(const char* formula) {
    std::vector<Peak> out;
    IsoOrderedGenerator g(formula);
    while (g.advanceToNextConfiguration())
        out.push_back({g.mass(), g.prob()});
    return out;
}

// Full enumeration via the layered generator run to exhaustion.
inline std::vector<Peak> enumerate_layered_full(const char* formula) {
    std::vector<Peak> out;
    IsoLayeredGenerator g(formula, 1000, 1000);
    while (g.advanceToNextConfiguration())
        out.push_back({g.mass(), g.prob()});
    return out;
}

inline double total_prob(const std::vector<Peak>& v) {
    double s = 0.0;
    for (const Peak& p : v) s += p.prob;
    return s;
}

// True if two peak sets describe the same distribution up to floating-point
// slop: identical count, and each canonically-ordered peak agrees within tol.
// Different generators sum/multiply marginals in different orders, so exact
// bit-identity is NOT expected across generators — hence the tolerance.
inline bool peaks_close(std::vector<Peak> a, std::vector<Peak> b, double tol = 1e-9) {
    if (a.size() != b.size()) return false;
    sort_peaks(a);
    sort_peaks(b);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i].mass - b[i].mass) > tol) return false;
        if (std::fabs(a[i].prob - b[i].prob) > tol) return false;
    }
    return true;
}

// Small, fully-enumerable molecules: used for cross-generator set-equality where
// we materialize the entire support three ways and compare.
inline const std::vector<const char*>& small_formulas() {
    static const std::vector<const char*> f = {
        "H2O1", "C2", "C10", "C6H12O6", "Fe2O3", "H2O2", "N2", "S1", "Cl2", "Br2",
        "C2H6O1", "C1H4", "N1H3", "C4H10", "Sn1", "Se1"
    };
    return f;
}

// Larger molecules: for threshold-bounded and total-prob comparisons where full
// enumeration would be intractable.
inline const std::vector<const char*>& big_formulas() {
    static const std::vector<const char*> f = {
        "C50H100N20O20S5", "C100H200N40O40", "C520H817N139O147S8",
        "C378H629N105O118S1", "H1000O500", "C1000"
    };
    return f;
}

}  // namespace test_helpers

#endif  // ISOSPEC_TEST_HELPERS_H
