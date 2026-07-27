// Generator semantics: the iteration contract each generator offers, and the
// knobs (thresholds, marginal reordering, table sizes, layer stepping) that are
// supposed to change performance without changing results.

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// A doctest assertion costs far more than the generator step it is checking, so
// loops that run over millions of configurations accumulate into one of these
// and assert once at the end.
struct Tally {
    std::size_t checked = 0;
    std::size_t failed = 0;
    void expect(bool ok) { ++checked; if (!ok) ++failed; }
    bool clean() const { return failed == 0 && checked > 0; }
};

// Is this signature a valid configuration of the molecule, and does it weigh
// what the generator says it does?
bool signature_ok(const std::vector<int>& sig, const std::vector<ElementSpec>& es,
                  double expected_mass) {
    std::size_t off = 0;
    double mass = 0.0;
    for (const ElementSpec& e : es) {
        int total = 0;
        for (std::size_t j = 0; j < e.isotopes.size(); ++j) {
            if (sig[off + j] < 0) return false;
            total += sig[off + j];
            mass += sig[off + j] * e.isotopes[j].mass;
        }
        if (total != e.count) return false;
        off += e.isotopes.size();
    }
    if (off != sig.size()) return false;
    return std::fabs(mass - expected_mass) <= 1e-9 * std::fabs(expected_mass);
}

}  // namespace

TEST_CASE("threshold generator: mass, prob and lprob stay consistent") {
    for (const char* f : {"H2O1", "C10H20O2", "C50H100N20O20S5"}) {
        INFO("formula=" << f);
        IsoThresholdGenerator g(Iso(f), 1e-6, true);
        Tally t;
        while (g.advanceToNextConfiguration())
            t.expect(g.prob() >= 1e-6 && g.mass() > 0.0 &&
                     std::fabs(g.prob() - std::exp(g.lprob())) <= 1e-9 * g.prob());
        CHECK(t.clean());
    }
}

TEST_CASE("threshold generator: signatures describe the reported peak") {
    const char* formula = "C10H20O2";
    const std::vector<ElementSpec> es = elements_of(formula);
    Iso iso(formula);
    const int allDim = iso.getAllDim();

    IsoThresholdGenerator g(std::move(iso), 1e-6, true);
    std::vector<int> sig(allDim);
    std::set<std::vector<int>> seen;
    Tally t;
    while (g.advanceToNextConfiguration()) {
        g.get_conf_signature(sig.data());
        // Every configuration is valid, and emitted exactly once.
        t.expect(signature_ok(sig, es, g.mass()) && seen.insert(sig).second);
    }
    CHECK(t.clean());
    CHECK(seen.size() > 10);
}

TEST_CASE("threshold generator: count_confs predicts, reset replays") {
    for (const char* f : small_formulas()) {
        INFO("formula=" << f);
        IsoThresholdGenerator g(Iso(f), 1e-10, true);

        const std::size_t predicted = g.count_confs();
        g.reset();

        std::vector<Peak> first;
        while (g.advanceToNextConfiguration()) first.push_back({g.mass(), g.prob()});
        CHECK(first.size() == predicted);

        // A reset generator replays the identical sequence, in the same order.
        g.reset();
        std::vector<Peak> second;
        while (g.advanceToNextConfiguration()) second.push_back({g.mass(), g.prob()});
        REQUIRE(second.size() == first.size());
        Tally t;
        for (std::size_t i = 0; i < first.size(); ++i)
            t.expect(first[i].mass == second[i].mass && first[i].prob == second[i].prob);
        CHECK(t.clean());
    }
}

TEST_CASE("threshold generator: terminate_search stops the iteration") {
    IsoThresholdGenerator g(Iso("C100H200N40O40"), 1e-10, true);
    for (int i = 0; i < 5; ++i) REQUIRE(g.advanceToNextConfiguration());

    g.terminate_search();
    CHECK_FALSE(g.advanceToNextConfiguration());
    // Still terminated on a second ask.
    CHECK_FALSE(g.advanceToNextConfiguration());

    // ... and reset() brings it back to life.
    g.reset();
    CHECK(g.advanceToNextConfiguration());
}

TEST_CASE("threshold generator: marginal reordering does not change the peak set") {
    for (const char* f : {"C10H20O2", "C50H100N20O20S5", "H1000O500"}) {
        INFO("formula=" << f);
        std::vector<Peak> reordered, plain;

        IsoThresholdGenerator a(Iso(f), 1e-6, true, 1000, 1000, true);
        while (a.advanceToNextConfiguration()) reordered.push_back({a.mass(), a.prob()});

        IsoThresholdGenerator b(Iso(f), 1e-6, true, 1000, 1000, false);
        while (b.advanceToNextConfiguration()) plain.push_back({b.mass(), b.prob()});

        CHECK(reordered.size() == plain.size());
        CHECK(peaks_close(reordered, plain, 1e-12));
    }
}

TEST_CASE("threshold generator: table and hash sizes do not change the result") {
    std::vector<Peak> reference = enumerate_threshold_full("C4H8O2");
    for (int tab : {1, 7, 1000, 100000}) {
        for (int hash : {1, 7, 1000}) {
            CAPTURE(tab);
            CAPTURE(hash);
            IsoThresholdGenerator g(Iso("C4H8O2"), 0.0, true, tab, hash);
            std::vector<Peak> got;
            while (g.advanceToNextConfiguration()) got.push_back({g.mass(), g.prob()});
            CHECK(peaks_close(got, reference, 1e-12));
        }
    }
}

TEST_CASE("threshold generator: relative threshold is scaled by the mode") {
    Iso iso("C50H100N20O20S5");
    const double mode = std::exp(iso.getModeLProb());

    std::vector<Peak> relative, absolute;
    IsoThresholdGenerator a(Iso("C50H100N20O20S5"), 1e-3, false);
    while (a.advanceToNextConfiguration()) relative.push_back({a.mass(), a.prob()});
    IsoThresholdGenerator b(Iso("C50H100N20O20S5"), 1e-3 * mode, true);
    while (b.advanceToNextConfiguration()) absolute.push_back({b.mass(), b.prob()});

    CHECK(relative.size() == absolute.size());
    CHECK(peaks_close(relative, absolute, 1e-15));
}

TEST_CASE("threshold generator: no_carry advance walks one marginal run") {
    // advanceToNextConfiguration_no_carry() is the inner loop the SIMD fill
    // drives: it walks the current run and stops instead of carrying, and a
    // manual carry() then continues exactly where the plain generator would.
    IsoThresholdGenerator manual(Iso("C10H20O2"), 1e-10, true);
    IsoThresholdGenerator plain(Iso("C10H20O2"), 1e-10, true);

    Tally t;
    bool more = true;
    while (more) {
        while (manual.advanceToNextConfiguration_no_carry())
            t.expect(plain.advanceToNextConfiguration() &&
                     manual.mass() == plain.mass() && manual.prob() == plain.prob());
        more = manual.carry();
        if (more)
            t.expect(plain.advanceToNextConfiguration() &&
                     manual.mass() == plain.mass() && manual.prob() == plain.prob());
    }
    CHECK_FALSE(plain.advanceToNextConfiguration());
    CHECK(t.clean());
}

TEST_CASE("layered generator: layers descend and cover everything") {
    // Run to exhaustion, so the molecule has to be one whose *entire* support
    // is enumerable — 1.2M configurations here.
    IsoLayeredGenerator g(Iso("C20H40N4O6S2"));

    double previous_threshold = std::numeric_limits<double>::infinity();
    double total = 0.0;
    Tally t;
    while (g.advanceToNextConfiguration()) {
        const double threshold = g.get_currentLThreshold();
        // The prob/exp(lprob) agreement is only meaningful above the denormal
        // range, where exp() still has full relative precision.
        t.expect(threshold <= previous_threshold + 1e-12 &&
                 (g.prob() < 1e-300 ||
                  std::fabs(g.prob() - std::exp(g.lprob())) <= 1e-9 * g.prob()));
        previous_threshold = threshold;
        total += g.prob();
    }
    CHECK(t.clean());
    CHECK(total == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("layered generator: within-layer iteration plus nextLayer covers the same set") {
    std::vector<Peak> stepwise;
    {
        IsoLayeredGenerator g(Iso("C10H20O2"));
        bool more = true;
        while (more) {
            while (g.advanceToNextConfigurationWithinLayer())
                stepwise.push_back({g.mass(), g.prob()});
            more = g.nextLayer(-2.0);
        }
    }

    std::vector<Peak> direct;
    {
        IsoLayeredGenerator g(Iso("C10H20O2"));
        while (g.advanceToNextConfiguration()) direct.push_back({g.mass(), g.prob()});
    }

    CHECK(stepwise.size() == direct.size());
    CHECK(peaks_close(stepwise, direct, 1e-15));
}

TEST_CASE("layered generator: terminate_search and reordering") {
    IsoLayeredGenerator g(Iso("C20H40N4O6S2"));
    for (int i = 0; i < 5; ++i) REQUIRE(g.advanceToNextConfiguration());
    g.terminate_search();
    CHECK_FALSE(g.advanceToNextConfiguration());

    // reorder_marginals and the probability hint are performance knobs only.
    std::vector<Peak> reference = enumerate_layered_full("C10H20O2");
    for (bool reorder : {true, false}) {
        for (double hint : {0.5, 0.9, 0.99, 0.9999}) {
            CAPTURE(reorder);
            CAPTURE(hint);
            IsoLayeredGenerator h(Iso("C10H20O2"), 1000, 1000, reorder, hint);
            std::vector<Peak> got;
            while (h.advanceToNextConfiguration()) got.push_back({h.mass(), h.prob()});
            CHECK(got.size() == reference.size());
            CHECK(peaks_close(got, reference, 1e-12));
        }
    }
}

TEST_CASE("layered generator: signatures describe the reported peak") {
    const char* formula = "C6H12O6";
    const std::vector<ElementSpec> es = elements_of(formula);
    Iso iso(formula);
    const int allDim = iso.getAllDim();

    IsoLayeredGenerator g(std::move(iso));
    std::vector<int> sig(allDim);
    std::set<std::vector<int>> seen;
    Tally t;
    while (g.advanceToNextConfiguration()) {
        g.get_conf_signature(sig.data());
        t.expect(signature_ok(sig, es, g.mass()) && seen.insert(sig).second);
    }
    CHECK(t.clean());
    CHECK(seen.size() > 10);
}

TEST_CASE("ordered generator: strictly descending, complete, and signed") {
    const char* formula = "C6H12O6";
    const std::vector<ElementSpec> es = elements_of(formula);
    Iso iso(formula);
    const int allDim = iso.getAllDim();

    IsoOrderedGenerator g(std::move(iso));
    std::vector<int> sig(allDim);
    double last = std::numeric_limits<double>::infinity();
    double total = 0.0;
    Tally t;
    while (g.advanceToNextConfiguration()) {
        g.get_conf_signature(sig.data());
        t.expect(g.prob() <= last + 1e-15 &&
                 std::fabs(g.prob() - std::exp(g.lprob())) <= 1e-9 * g.prob() &&
                 signature_ok(sig, es, g.mass()));
        last = g.prob();
        total += g.prob();
    }
    CHECK(t.clean());
    CHECK(total == doctest::Approx(1.0).epsilon(1e-9));

    // The first configuration is the mode.
    IsoOrderedGenerator h{Iso(formula)};
    REQUIRE(h.advanceToNextConfiguration());
    CHECK(h.lprob() == doctest::Approx(Iso(formula).getModeLProb()));
}

TEST_CASE("ordered generator: table and hash sizes do not change the result") {
    std::vector<Peak> reference = enumerate_ordered_full("C4H8O2");
    for (int tab : {1, 13, 10000}) {
        for (int hash : {1, 13, 10000}) {
            CAPTURE(tab);
            CAPTURE(hash);
            IsoOrderedGenerator g(Iso("C4H8O2"), tab, hash);
            std::vector<Peak> got;
            while (g.advanceToNextConfiguration()) got.push_back({g.mass(), g.prob()});
            CHECK(peaks_close(got, reference, 1e-12));
        }
    }
}

TEST_CASE("stochastic generator: counts sum to the sample size") {
    for (std::size_t molecules : {std::size_t(1), std::size_t(10), std::size_t(100000)}) {
        CAPTURE(molecules);
        std::mt19937 rng(20260727);
        IsoStochasticGenerator g(Iso("C100H200N40O40"), molecules, 0.9999, 5.0, rng);

        double total = 0.0;
        Tally t;
        while (g.advanceToNextConfiguration()) {
            const double count = g.prob();
            t.expect(count >= 1.0 && count == std::round(count) &&
                     g.count() == static_cast<std::size_t>(count) &&
                     std::fabs(g.lprob() - std::log(count)) < 1e-9 &&
                     g.mass() > 0.0);
            total += count;
        }
        CHECK(t.clean());
        CHECK(total == doctest::Approx(static_cast<double>(molecules)));
    }
}

TEST_CASE("stochastic generator: the same seed gives the same sample") {
    auto sample = [](unsigned seed) {
        std::mt19937 rng(seed);
        IsoStochasticGenerator g(Iso("C100H200N40O40"), 10000, 0.9999, 5.0, rng);
        std::vector<Peak> out;
        while (g.advanceToNextConfiguration()) out.push_back({g.mass(), g.prob()});
        return out;
    };

    std::vector<Peak> a = sample(42);
    std::vector<Peak> b = sample(42);
    REQUIRE(a.size() == b.size());
    Tally t;
    for (std::size_t i = 0; i < a.size(); ++i)
        t.expect(a[i].mass == b[i].mass && a[i].prob == b[i].prob);
    CHECK(t.clean());

    // A different seed gives a different sample (with overwhelming probability).
    std::vector<Peak> c = sample(43);
    bool identical = c.size() == a.size();
    if (identical)
        for (std::size_t i = 0; i < a.size() && identical; ++i)
            identical = (a[i].prob == c[i].prob);
    CHECK_FALSE(identical);
}

TEST_CASE("stochastic generator: the sample tracks the true distribution") {
    // With 200k draws the empirical mean must be very close to the theoretical
    // one, and every sampled mass must be a real isotopologue mass.
    Iso iso("C100H200N40O40");
    const double theoretical = iso.getTheoreticalAverageMass();
    const double lightest = iso.getLightestPeakMass();
    const double heaviest = iso.getHeaviestPeakMass();

    std::mt19937 rng(7);
    IsoStochasticGenerator g(Iso("C100H200N40O40"), 200000, 0.9999, 5.0, rng);
    double sum = 0.0, count = 0.0;
    Tally t;
    while (g.advanceToNextConfiguration()) {
        t.expect(g.mass() >= lightest - 1e-9 && g.mass() <= heaviest + 1e-9);
        sum += g.mass() * g.prob();
        count += g.prob();
    }
    CHECK(t.clean());
    CHECK(count == doctest::Approx(200000.0));
    CHECK(sum / count == doctest::Approx(theoretical).epsilon(1e-4));
}

TEST_CASE("stochastic sampling on top of the ordered generator") {
    // IsoStochasticGenerator wraps the layered generator, but the template is
    // also instantiated over the ordered one; that combination samples in
    // descending-probability order.
    std::mt19937 rng(20260727);
    IsoStochasticGeneratorTemplate<IsoOrderedGeneratorTemplate<MarginalTrek>>
        g(Iso("C10H20O2"), 10000, 0.9999, 5.0, rng);

    double total = 0.0;
    std::size_t peaks = 0;
    Tally t;
    while (g.advanceToNextConfiguration()) {
        const double count = g.prob();
        t.expect(count >= 1.0 && count == std::round(count) && g.mass() > 0.0);
        total += count;
        ++peaks;
    }
    CHECK(t.clean());
    CHECK(peaks > 0);
    CHECK(total == doctest::Approx(10000.0));
}

TEST_CASE("generators agree on a single-atom molecule") {
    // The degenerate case: one element, one atom — the support is the isotope
    // list itself.
    for (const char* f : {"C1", "S1", "Sn1", "F1"}) {
        INFO("formula=" << f);
        const std::vector<Isotope> isotopes = elements_of(f)[0].isotopes;

        std::vector<Peak> expected;
        for (const Isotope& i : isotopes) expected.push_back({i.mass, i.prob});

        CHECK(peaks_close(enumerate_threshold_full(f), expected, 1e-12));
        CHECK(peaks_close(enumerate_ordered_full(f), expected, 1e-12));
        CHECK(peaks_close(enumerate_layered_full(f), expected, 1e-12));
    }
}

TEST_CASE("generators handle a single-isotope molecule") {
    // Fluorine is monoisotopic: exactly one configuration, probability 1.
    for (const char* f : {"F1", "F100", "Na10P3"}) {
        INFO("formula=" << f);
        std::vector<Peak> thr = enumerate_threshold_full(f);
        REQUIRE(thr.size() == 1);
        CHECK(thr[0].prob == doctest::Approx(1.0));
        CHECK(thr[0].mass == doctest::Approx(Iso(f).getMonoisotopicPeakMass()));

        CHECK(enumerate_ordered_full(f).size() == 1);
        CHECK(enumerate_layered_full(f).size() == 1);
    }
}

TEST_CASE("threshold generator with a threshold above every peak yields nothing") {
    IsoThresholdGenerator g(Iso("C100H200N40O40"), 0.5, true);
    CHECK(g.count_confs() == 0);
    g.reset();
    CHECK_FALSE(g.advanceToNextConfiguration());
}
