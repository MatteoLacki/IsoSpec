// Shared helpers for the doctest-based C++ test suite.
//
// Header-only: every test translation unit includes this; the definitions are
// inline so they can be pulled into multiple TUs without ODR violations.  The
// IsoSpec symbols themselves come from unity-build.cpp, linked into the test
// binary (see the Makefile).

#ifndef ISOSPEC_TEST_HELPERS_H
#define ISOSPEC_TEST_HELPERS_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include "element_tables.h"

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

// ---------------------------------------------------------------------------
// Independent access to the built-in isotope tables.
//
// The library reaches the same data through parse_formula(); the helpers below
// go straight to element_tables.h instead, so a test that compares one against
// the other is checking the parser, not restating it.
// ---------------------------------------------------------------------------

struct Isotope {
    double mass;
    double prob;
};

//! All stable isotopes of the element with the given symbol, in table order.
//! Returns an empty vector for an unknown symbol.
inline std::vector<Isotope> element_isotopes(const char* symbol) {
    std::vector<Isotope> out;
    for (int i = 0; i < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; ++i) {
        if (std::strcmp(elem_table_symbol[i], symbol) != 0) continue;
        const int id = elem_table_ID[i];
        for (int j = i; j < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES && elem_table_ID[j] == id; ++j)
            out.push_back({elem_table_mass[j], elem_table_probability[j]});
        break;
    }
    return out;
}

//! One element of a molecule: how many atoms, and of which isotopes.
struct ElementSpec {
    int count;
    std::vector<Isotope> isotopes;
};

//! Parse a simple "Symbol<count>..." formula into per-element specs using the
//! raw tables (no parse_formula involved).
inline std::vector<ElementSpec> elements_of(const std::string& formula) {
    std::vector<ElementSpec> out;
    std::size_t i = 0;
    while (i < formula.size()) {
        std::size_t s = i;
        while (i < formula.size() && std::isalpha(static_cast<unsigned char>(formula[i]))) ++i;
        const std::string sym = formula.substr(s, i - s);
        std::size_t d = i;
        while (i < formula.size() && std::isdigit(static_cast<unsigned char>(formula[i]))) ++i;
        const int cnt = std::stoi(formula.substr(d, i - d));
        out.push_back({cnt, element_isotopes(sym.c_str())});
    }
    return out;
}

// ---------------------------------------------------------------------------
// Brute-force distribution oracle.
//
// Builds the isotopic distribution by explicit convolution: start from a single
// zero-mass, unit-probability peak and, for every atom in turn, convolve with
// that atom's isotope distribution.  This shares no code with IsoSpec beyond
// the isotope tables — it is an independent computation of the same quantity,
// exponential in the number of atoms, so keep the molecules tiny.
// ---------------------------------------------------------------------------

inline std::vector<Peak> brute_force_distribution(const std::vector<ElementSpec>& elems) {
    std::vector<Peak> dist = {{0.0, 1.0}};
    for (const ElementSpec& e : elems)
        for (int a = 0; a < e.count; ++a) {
            std::vector<Peak> next;
            next.reserve(dist.size() * e.isotopes.size());
            for (const Peak& p : dist)
                for (const Isotope& iso : e.isotopes)
                    next.push_back({p.mass + iso.mass, p.prob * iso.prob});
            dist.swap(next);
        }
    return dist;
}

//! Same, but merging peaks of exactly equal mass (the enumeration of a
//! molecule yields one peak per *configuration*, whereas the brute force above
//! yields one per atom-labelling; equal-mass merging is not what we want there,
//! so this is provided separately for mass-domain comparisons).
inline std::vector<Peak> merge_equal_masses(std::vector<Peak> v, double tol = 1e-9) {
    sort_peaks(v);
    std::vector<Peak> out;
    for (const Peak& p : v) {
        if (!out.empty() && std::fabs(out.back().mass - p.mass) <= tol)
            out.back().prob += p.prob;
        else
            out.push_back(p);
    }
    return out;
}

//! Envelope contents as a peak vector.
inline std::vector<Peak> envelope_peaks(const FixedEnvelope& env) {
    std::vector<Peak> out;
    out.reserve(env.confs_no());
    for (std::size_t i = 0; i < env.confs_no(); ++i)
        out.push_back({env.masses()[i], env.probs()[i]});
    return out;
}

//! Probability-weighted mean mass of a peak set.
inline double mean_mass(const std::vector<Peak>& v) {
    double m = 0.0, p = 0.0;
    for (const Peak& q : v) { m += q.mass * q.prob; p += q.prob; }
    return m / p;
}

}  // namespace test_helpers

#endif  // ISOSPEC_TEST_HELPERS_H
