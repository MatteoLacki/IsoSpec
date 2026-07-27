// Analytic oracles: independent ground truth, not another generator.
//
// For a molecule made of a single element with controlled isotope
// probabilities, the isotopologue distribution is exactly multinomial.  We
// build such an Iso from raw arrays (so the probabilities are ours, not the
// built-in element tables) and check IsoSpec reproduces the closed form.

#include <cmath>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;

namespace {

// log C(n, k)
double log_choose(int n, int k) {
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

// Enumerate all probabilities of a single-element molecule via FromThreshold(0),
// sorted ascending.
std::vector<double> isospec_probs(int isotopeNo, int atomCount,
                                  const double* masses, const double* probs) {
    int isotopeNumbers[1] = {isotopeNo};
    int atomCounts[1] = {atomCount};
    const double* IM[1] = {masses};
    const double* IP[1] = {probs};
    Iso iso(1, isotopeNumbers, atomCounts, IM, IP);
    FixedEnvelope env = FixedEnvelope::FromThreshold(std::move(iso), 0.0, true, false);
    std::vector<double> out(env.probs(), env.probs() + env.confs_no());
    std::sort(out.begin(), out.end());
    return out;
}

void check_close_multiset(std::vector<double> got, std::vector<double> want, double tol) {
    std::sort(want.begin(), want.end());
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < got.size(); ++i)
        CHECK(std::fabs(got[i] - want[i]) < tol);
}

}  // namespace

TEST_CASE("single element, 2 isotopes: binomial distribution") {
    const double masses[2] = {10.0, 11.0};
    const double probs[2] = {0.7, 0.3};

    for (int n : {1, 2, 5, 10, 20}) {
        CAPTURE(n);
        std::vector<double> want;
        for (int k = 0; k <= n; ++k)
            want.push_back(std::exp(log_choose(n, k) +
                                    (n - k) * std::log(probs[0]) +
                                    k * std::log(probs[1])));
        check_close_multiset(isospec_probs(2, n, masses, probs), want, 1e-12);

        // Distribution normalizes to 1.
        double s = 0.0;
        for (double p : want) s += p;
        CHECK(std::fabs(s - 1.0) < 1e-12);
    }
}

TEST_CASE("single element, 3 isotopes: trinomial distribution") {
    const double masses[3] = {20.0, 21.0, 22.0};
    const double probs[3] = {0.5, 0.3, 0.2};

    for (int n : {1, 2, 4, 8}) {
        CAPTURE(n);
        std::vector<double> want;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n - i; ++j) {
                int k = n - i - j;  // counts of the three isotopes
                double lp = std::lgamma(n + 1.0) - std::lgamma(i + 1.0) -
                            std::lgamma(j + 1.0) - std::lgamma(k + 1.0) +
                            i * std::log(probs[0]) + j * std::log(probs[1]) +
                            k * std::log(probs[2]);
                want.push_back(std::exp(lp));
            }
        check_close_multiset(isospec_probs(3, n, masses, probs), want, 1e-12);
    }
}
