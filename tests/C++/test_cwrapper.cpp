// The extern "C" ABI (cwrapper.h) — the surface the Python (cffi) and other
// non-C++ bindings actually use.
//
// Two things are checked throughout: that each entry point computes the same
// thing as the C++ API it wraps, and that bad input comes back as the
// documented out-of-band value (NULL / NaN / 0) instead of an exception
// unwinding across the ABI.

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "doctest.h"
#include "cwrapper.h"
#include "fasta.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// A C-side Iso handle that always gets deleted.
struct IsoHandle {
    void* p;
    explicit IsoHandle(void* q) : p(q) {}
    ~IsoHandle() { deleteIso(p); }
    IsoHandle(const IsoHandle&) = delete;
    IsoHandle& operator=(const IsoHandle&) = delete;
};

// Build a C-side Iso for a formula, going through the raw tables (the C ABI has
// no formula parser of its own).
void* c_iso_for(const char* formula, std::vector<int>& isotopeNumbers,
                std::vector<int>& atomCounts, std::vector<double>& masses,
                std::vector<double>& probs) {
    const std::vector<ElementSpec> es = elements_of(formula);
    for (const ElementSpec& e : es) {
        isotopeNumbers.push_back(static_cast<int>(e.isotopes.size()));
        atomCounts.push_back(e.count);
        for (const Isotope& i : e.isotopes) { masses.push_back(i.mass); probs.push_back(i.prob); }
    }
    return setupIso(static_cast<int>(es.size()), isotopeNumbers.data(), atomCounts.data(),
                    masses.data(), probs.data());
}

}  // namespace

TEST_CASE("setupIso and the Iso accessors mirror the C++ class") {
    std::vector<int> isoNo, atomCnt;
    std::vector<double> masses, probs;
    IsoHandle h(c_iso_for("C6H12O6", isoNo, atomCnt, masses, probs));
    REQUIRE(h.p != nullptr);

    Iso reference("C6H12O6");

    CHECK(getLightestPeakMassIso(h.p) == doctest::Approx(reference.getLightestPeakMass()));
    CHECK(getHeaviestPeakMassIso(h.p) == doctest::Approx(reference.getHeaviestPeakMass()));
    CHECK(getMonoisotopicPeakMassIso(h.p) == doctest::Approx(reference.getMonoisotopicPeakMass()));
    CHECK(getModeMassIso(h.p) == doctest::Approx(reference.getModeMass()));
    CHECK(getTheoreticalAverageMassIso(h.p) == doctest::Approx(reference.getTheoreticalAverageMass()));

    CHECK(getLightestPeakLProbIso(h.p) == doctest::Approx(reference.getLightestPeakLProb()));
    CHECK(getHeaviestPeakLProbIso(h.p) == doctest::Approx(reference.getHeaviestPeakLProb()));
    CHECK(getMonoisotopicPeakLProbIso(h.p) == doctest::Approx(reference.getMonoisotopicPeakLProb()));
    CHECK(getModeLProbIso(h.p) == doctest::Approx(reference.getModeLProb()));

    CHECK(getIsoVariance(h.p) == doctest::Approx(reference.variance()));
    CHECK(getIsoStddev(h.p) == doctest::Approx(reference.stddev()));
    CHECK(getIsoStddev(h.p) == doctest::Approx(std::sqrt(getIsoVariance(h.p))));

    // Signatures.
    const int allDim = reference.getAllDim();
    std::vector<int> c_sig(allDim, -1), cpp_sig(allDim, -1);

    getLightestPeakSignature(h.p, c_sig.data());
    reference.getLightestPeakSignature(cpp_sig.data());
    CHECK(c_sig == cpp_sig);

    getHeaviestPeakSignature(h.p, c_sig.data());
    reference.getHeaviestPeakSignature(cpp_sig.data());
    CHECK(c_sig == cpp_sig);

    getMonoisotopicPeakSignature(h.p, c_sig.data());
    reference.getMonoisotopicPeakSignature(cpp_sig.data());
    CHECK(c_sig == cpp_sig);
}

TEST_CASE("getMarginalLogSizeEstimates returns a freeable array") {
    std::vector<int> isoNo, atomCnt;
    std::vector<double> masses, probs;
    IsoHandle h(c_iso_for("C100H2O1", isoNo, atomCnt, masses, probs));
    REQUIRE(h.p != nullptr);

    double* est = getMarginalLogSizeEstimates(h.p, 0.99);
    REQUIRE(est != nullptr);

    Iso reference("C100H2O1");
    std::vector<double> expected(reference.getDimNumber());
    reference.saveMarginalLogSizeEstimates(expected.data(), 0.99);
    for (int i = 0; i < reference.getDimNumber(); ++i)
        CHECK(est[i] == doctest::Approx(expected[i]));

    freeReleasedArray(est);
}

TEST_CASE("isoFromFasta matches Iso::FromFASTA") {
    void* h = isoFromFasta("PEPTIDE", false, true);
    REQUIRE(h != nullptr);
    Iso reference = Iso::FromFASTA("PEPTIDE", false, true);
    CHECK(getMonoisotopicPeakMassIso(h) == doctest::Approx(reference.getMonoisotopicPeakMass()));
    CHECK(getTheoreticalAverageMassIso(h) == doctest::Approx(reference.getTheoreticalAverageMass()));
    deleteIso(h);

    void* nominal = isoFromFasta("PEPTIDE", true, false);
    REQUIRE(nominal != nullptr);
    Iso ref_nominal = Iso::FromFASTA("PEPTIDE", true, false);
    CHECK(getMonoisotopicPeakMassIso(nominal) ==
          doctest::Approx(ref_nominal.getMonoisotopicPeakMass()));
    deleteIso(nominal);
}

TEST_CASE("bad input comes back as NULL, not as an exception") {
    const int isoNumbers[1] = {2};
    const int atomCounts[1] = {10};
    const double masses[2] = {1.0, 2.0};
    const double zero_prob[2] = {0.0, 1.0};

    void* bad = nullptr;
    CHECK_NOTHROW(bad = setupIso(1, isoNumbers, atomCounts, masses, zero_prob));
    CHECK(bad == nullptr);
    CHECK_NOTHROW(deleteIso(nullptr));

    // A zero-dimensional Iso is constructible, but no generator accepts it.
    void* empty = setupIso(0, nullptr, nullptr, nullptr, nullptr);
    REQUIRE(empty != nullptr);
    void* gen = nullptr;
    CHECK_NOTHROW(gen = setupIsoThresholdGenerator(empty, 0.01, true, 1000, 1000, true));
    CHECK(gen == nullptr);
    CHECK_NOTHROW(deleteIsoThresholdGenerator(nullptr));
    deleteIso(empty);
}

TEST_CASE("threshold generator through the C ABI") {
    std::vector<int> isoNo, atomCnt;
    std::vector<double> masses, probs;
    void* iso = c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs);
    REQUIRE(iso != nullptr);

    void* gen = setupIsoThresholdGenerator(iso, 1e-6, true, 1000, 1000, true);
    REQUIRE(gen != nullptr);

    std::vector<Peak> got;
    const int allDim = Iso("C10H20O2").getAllDim();
    std::vector<int> sig(allDim);
    while (advanceToNextConfigurationIsoThresholdGenerator(gen)) {
        const double m = massIsoThresholdGenerator(gen);
        const double p = probIsoThresholdGenerator(gen);
        const double lp = lprobIsoThresholdGenerator(gen);
        CHECK(p == doctest::Approx(std::exp(lp)));
        CHECK(p >= 1e-6);
        get_conf_signatureIsoThresholdGenerator(gen, sig.data());
        int total = 0;
        for (int v : sig) { CHECK(v >= 0); total += v; }
        CHECK(total == 10 + 20 + 2);
        got.push_back({m, p});
    }
    deleteIsoThresholdGenerator(gen);
    deleteIso(iso);

    // Must match the C++ generator exactly.
    std::vector<Peak> expected;
    IsoThresholdGenerator cpp_gen(Iso("C10H20O2"), 1e-6, true);
    while (cpp_gen.advanceToNextConfiguration())
        expected.push_back({cpp_gen.mass(), cpp_gen.prob()});
    CHECK(got.size() == expected.size());
    CHECK(peaks_close(got, expected, 1e-12));
}

TEST_CASE("layered, ordered and stochastic generators through the C ABI") {
    std::vector<int> isoNo, atomCnt;
    std::vector<double> masses, probs;
    const int allDim = Iso("C10H20O2").getAllDim();
    std::vector<int> sig(allDim);

    SUBCASE("layered") {
        void* iso = c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs);
        void* gen = setupIsoLayeredGenerator(iso, 1000, 1000, true, 0.99);
        REQUIRE(gen != nullptr);
        std::size_t n = 0;
        double total = 0.0;
        while (advanceToNextConfigurationIsoLayeredGenerator(gen)) {
            total += probIsoLayeredGenerator(gen);
            CHECK(massIsoLayeredGenerator(gen) > 0.0);
            CHECK(lprobIsoLayeredGenerator(gen) <= 0.0);
            get_conf_signatureIsoLayeredGenerator(gen, sig.data());
            ++n;
        }
        CHECK(n > 0);
        CHECK(total == doctest::Approx(1.0).epsilon(1e-9));
        deleteIsoLayeredGenerator(gen);
        deleteIso(iso);
    }

    SUBCASE("ordered") {
        void* iso = c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs);
        void* gen = setupIsoOrderedGenerator(iso, 1000, 1000);
        REQUIRE(gen != nullptr);
        std::size_t n = 0;
        double last = 2.0, total = 0.0;
        while (advanceToNextConfigurationIsoOrderedGenerator(gen)) {
            const double p = probIsoOrderedGenerator(gen);
            CHECK(p <= last + 1e-15);
            last = p;
            total += p;
            CHECK(massIsoOrderedGenerator(gen) > 0.0);
            CHECK(lprobIsoOrderedGenerator(gen) <= 0.0);
            get_conf_signatureIsoOrderedGenerator(gen, sig.data());
            ++n;
        }
        CHECK(n > 0);
        CHECK(total == doctest::Approx(1.0).epsilon(1e-9));
        deleteIsoOrderedGenerator(gen);
        deleteIso(iso);
    }

    SUBCASE("stochastic") {
        void* iso = c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs);
        void* gen = setupIsoStochasticGenerator(iso, 10000, 0.9999, 5.0);
        REQUIRE(gen != nullptr);
        double total = 0.0;
        std::size_t n = 0;
        while (advanceToNextConfigurationIsoStochasticGenerator(gen)) {
            const double count = probIsoStochasticGenerator(gen);
            CHECK(count >= 1.0);
            CHECK(count == doctest::Approx(std::round(count)));
            CHECK(lprobIsoStochasticGenerator(gen) == doctest::Approx(std::log(count)));
            CHECK(massIsoStochasticGenerator(gen) > 0.0);
            get_conf_signatureIsoStochasticGenerator(gen, sig.data());
            total += count;
            ++n;
        }
        CHECK(n > 0);
        CHECK(total == doctest::Approx(10000.0));
        deleteIsoStochasticGenerator(gen);
        deleteIso(iso);
    }
}

TEST_CASE("fixed envelopes through the C ABI") {
    std::vector<int> isoNo, atomCnt;
    std::vector<double> masses, probs;

    SUBCASE("threshold envelope") {
        IsoHandle h(c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs));
        void* env = setupThresholdFixedEnvelope(h.p, 1e-6, true, false);
        REQUIRE(env != nullptr);

        FixedEnvelope reference = FixedEnvelope::FromThreshold(Iso("C10H20O2"), 1e-6, true, false);
        const std::size_t n = confs_noFixedEnvelope(env);
        CHECK(n == reference.confs_no());
        CHECK(getTotalProbOfEnvelope(env) == doctest::Approx(reference.get_total_prob()));
        CHECK(empiricAverageMass(env) == doctest::Approx(reference.empiric_average_mass()));
        CHECK(empiricVariance(env) == doctest::Approx(reference.empiric_variance()));
        CHECK(empiricStddev(env) == doctest::Approx(reference.empiric_stddev()));

        // massesFixedEnvelope/probsFixedEnvelope *release* the arrays to the
        // caller: a second call returns NULL and the caller must free them.
        const double* m = massesFixedEnvelope(env);
        const double* p = probsFixedEnvelope(env);
        REQUIRE(m != nullptr);
        REQUIRE(p != nullptr);
        CHECK(massesFixedEnvelope(env) == nullptr);
        std::vector<Peak> got;
        for (std::size_t i = 0; i < n; ++i) got.push_back({m[i], p[i]});
        CHECK(peaks_close(got, envelope_peaks(reference), 1e-12));
        freeReleasedArray(const_cast<double*>(m));
        freeReleasedArray(const_cast<double*>(p));
        deleteFixedEnvelope(env, false);
    }

    SUBCASE("threshold envelope with configurations") {
        IsoHandle h(c_iso_for("H2O1", isoNo, atomCnt, masses, probs));
        void* env = setupThresholdFixedEnvelope(h.p, 0.0, true, true);
        REQUIRE(env != nullptr);
        const std::size_t n = confs_noFixedEnvelope(env);
        CHECK(n > 0);
        const int* confs = confsFixedEnvelope(env);
        REQUIRE(confs != nullptr);
        const int allDim = 5;  // H:2 + O:3
        for (std::size_t i = 0; i < n; ++i) {
            int h_total = confs[i * allDim] + confs[i * allDim + 1];
            int o_total = confs[i * allDim + 2] + confs[i * allDim + 3] + confs[i * allDim + 4];
            CHECK(h_total == 2);
            CHECK(o_total == 1);
        }
        freeReleasedArray(const_cast<int*>(confs));
        deleteFixedEnvelope(env, false);
    }

    SUBCASE("total-prob envelope") {
        IsoHandle h(c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs));
        void* env = setupTotalProbFixedEnvelope(h.p, 0.99, true, false);
        REQUIRE(env != nullptr);
        CHECK(getTotalProbOfEnvelope(env) >= 0.99);
        deleteFixedEnvelope(env, false);
    }

    SUBCASE("stochastic envelope") {
        IsoHandle h(c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs));
        void* env = setupStochasticFixedEnvelope(h.p, 5000, 0.9999, 5.0, false);
        REQUIRE(env != nullptr);
        CHECK(getTotalProbOfEnvelope(env) == doctest::Approx(5000.0));
        deleteFixedEnvelope(env, false);
    }

    SUBCASE("binned envelope") {
        IsoHandle h(c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs));
        void* env = setupBinnedFixedEnvelope(h.p, 0.99, 1.0, 0.0);
        REQUIRE(env != nullptr);
        CHECK(getTotalProbOfEnvelope(env) >= 0.99);
        const std::size_t n = confs_noFixedEnvelope(env);
        const double* m = massesFixedEnvelope(env);
        REQUIRE(m != nullptr);
        for (std::size_t i = 0; i < n; ++i) CHECK(m[i] == doctest::Approx(std::round(m[i])));
        freeReleasedArray(const_cast<double*>(m));
        deleteFixedEnvelope(env, false);
    }
}

TEST_CASE("envelopes built from caller-owned arrays") {
    double masses[3] = {1.0, 2.0, 3.0};
    double probs[3] = {0.2, 0.3, 0.5};

    void* env = setupFixedEnvelope(masses, probs, 3, false, false, NAN);
    REQUIRE(env != nullptr);
    CHECK(confs_noFixedEnvelope(env) == 3);
    CHECK(getTotalProbOfEnvelope(env) == doctest::Approx(1.0));
    CHECK(empiricAverageMass(env) == doctest::Approx(1.0 * 0.2 + 2.0 * 0.3 + 3.0 * 0.5));

    // Mutating operations write through to the caller's arrays.
    scaleEnvelope(env, 2.0);
    CHECK(probs[0] == doctest::Approx(0.4));
    normalizeEnvelope(env);
    CHECK(probs[0] == doctest::Approx(0.2));
    shiftMassEnvelope(env, 10.0);
    CHECK(masses[0] == doctest::Approx(11.0));
    sortEnvelopeByProb(env);
    CHECK(probs[0] <= probs[1]);
    sortEnvelopeByMass(env);
    CHECK(masses[0] <= masses[1]);

    // release_everything: the envelope must not free what it does not own.
    deleteFixedEnvelope(env, true);

    int confs[3] = {1, 2, 3};
    void* env2 = setupFixedEnvelopeWithConfs(masses, probs, confs, 3, 1, true, false, 1.0);
    REQUIRE(env2 != nullptr);
    CHECK(confs_noFixedEnvelope(env2) == 3);
    CHECK(getTotalProbOfEnvelope(env2) == doctest::Approx(1.0));  // taken from the argument
    deleteFixedEnvelope(env2, true);
}

TEST_CASE("envelope arithmetic through the C ABI") {
    double m1[2] = {1.0, 2.0}, p1[2] = {0.5, 0.5};
    double m2[2] = {10.0, 20.0}, p2[2] = {0.25, 0.75};
    void* a = setupFixedEnvelope(m1, p1, 2, false, false, NAN);
    void* b = setupFixedEnvelope(m2, p2, 2, false, false, NAN);
    REQUIRE(a != nullptr);
    REQUIRE(b != nullptr);

    void* sum = addEnvelopes(a, b);
    REQUIRE(sum != nullptr);
    CHECK(confs_noFixedEnvelope(sum) == 4);
    CHECK(getTotalProbOfEnvelope(sum) == doctest::Approx(2.0));

    void* conv = convolveEnvelopes(a, b);
    REQUIRE(conv != nullptr);
    CHECK(confs_noFixedEnvelope(conv) == 4);
    CHECK(getTotalProbOfEnvelope(conv) == doctest::Approx(1.0));
    CHECK(empiricAverageMass(conv) ==
          doctest::Approx(empiricAverageMass(a) + empiricAverageMass(b)));

    void* copy = copyFixedEnvelope(conv);
    REQUIRE(copy != nullptr);
    CHECK(confs_noFixedEnvelope(copy) == confs_noFixedEnvelope(conv));
    CHECK(getTotalProbOfEnvelope(copy) == doctest::Approx(getTotalProbOfEnvelope(conv)));

    void* const envelopes[2] = {a, b};
    const double intensities[2] = {2.0, 3.0};
    void* combo = linearCombination(envelopes, intensities, 2);
    REQUIRE(combo != nullptr);
    CHECK(confs_noFixedEnvelope(combo) == 4);
    CHECK(getTotalProbOfEnvelope(combo) == doctest::Approx(2.0 + 3.0));

    void* binned = binnedEnvelope(conv, 1.0, 0.0);
    REQUIRE(binned != nullptr);
    CHECK(getTotalProbOfEnvelope(binned) == doctest::Approx(getTotalProbOfEnvelope(conv)));

    deleteFixedEnvelope(binned, false);
    deleteFixedEnvelope(combo, false);
    deleteFixedEnvelope(copy, false);
    deleteFixedEnvelope(conv, false);
    deleteFixedEnvelope(sum, false);
    deleteFixedEnvelope(b, true);
    deleteFixedEnvelope(a, true);
}

TEST_CASE("distances through the C ABI") {
    double m1[2] = {0.0, 1.0}, p1[2] = {0.5, 0.5};
    double m2[2] = {2.0, 3.0}, p2[2] = {0.5, 0.5};
    void* a = setupFixedEnvelope(m1, p1, 2, false, false, NAN);
    void* b = setupFixedEnvelope(m2, p2, 2, false, false, NAN);

    CHECK(wassersteinDistance(a, b) == doctest::Approx(2.0));
    CHECK(orientedWassersteinDistance(a, b) ==
          doctest::Approx(-orientedWassersteinDistance(b, a)));
    CHECK(abyssalWassersteinDistance(a, b, 10.0, 1.0) == doctest::Approx(2.0));

    struct ws_match_res res = wassersteinMatch(a, b, 0.5, 1.0);
    CHECK(res.flow == doctest::Approx(0.0));
    CHECK(res.res1 == doctest::Approx(1.0));
    CHECK(res.res2 == doctest::Approx(1.0));

    // An unnormalized pair makes WassersteinDistance throw; the guard must turn
    // that into a NaN rather than let it cross the ABI.
    double m3[1] = {0.0};
    double p3[1] = {17.0};
    void* c = setupFixedEnvelope(m3, p3, 1, false, false, NAN);
    double d = 0.0;
    CHECK_NOTHROW(d = wassersteinDistance(a, c));
    CHECK(std::isnan(d));

    deleteFixedEnvelope(c, true);
    deleteFixedEnvelope(b, true);
    deleteFixedEnvelope(a, true);
}

TEST_CASE("resample through the C ABI") {
    std::vector<int> isoNo, atomCnt;
    std::vector<double> masses, probs;
    IsoHandle h(c_iso_for("C10H20O2", isoNo, atomCnt, masses, probs));
    void* env = setupTotalProbFixedEnvelope(h.p, 0.999, true, false);
    REQUIRE(env != nullptr);
    normalizeEnvelope(env);

    resampleEnvelope(env, 100000, 1.0);
    CHECK(getTotalProbOfEnvelope(env) == doctest::Approx(100000.0));

    const std::size_t n = confs_noFixedEnvelope(env);
    const double* p = probsFixedEnvelope(env);
    REQUIRE(p != nullptr);
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(p[i] >= 0.0);
        CHECK(p[i] == doctest::Approx(std::round(p[i])));
        total += p[i];
    }
    CHECK(total == doctest::Approx(100000.0));
    freeReleasedArray(const_cast<double*>(p));
    deleteFixedEnvelope(env, false);
}

TEST_CASE("array helpers") {
    double a[4] = {1.0, 2.0, 3.0, 4.0};
    array_add(a, 4, 10.0);
    CHECK(a[0] == doctest::Approx(11.0));
    CHECK(a[3] == doctest::Approx(14.0));

    array_mul(a, 4, 2.0);
    CHECK(a[0] == doctest::Approx(22.0));
    CHECK(a[3] == doctest::Approx(28.0));

    array_fma(a, 4, 0.5, 1.0);
    CHECK(a[0] == doctest::Approx(12.0));
    CHECK(a[3] == doctest::Approx(15.0));

    // Zero-length operations must be no-ops, not crashes.
    CHECK_NOTHROW(array_add(a, 0, 1.0));
    CHECK_NOTHROW(array_mul(a, 0, 1.0));
    CHECK_NOTHROW(array_fma(a, 0, 1.0, 1.0));
    CHECK(a[0] == doctest::Approx(12.0));
}

TEST_CASE("parse_fasta_c matches the C++ parser") {
    for (const char* seq : {"", "A", "PEPTIDE", "MKWVTFISLLLLFSSAYSRGV", "ae-da"}) {
        INFO("sequence='" << seq << "'");
        int c_counts[6], cpp_counts[6];
        parse_fasta_c(seq, c_counts);
        parse_fasta(seq, cpp_counts);
        for (int i = 0; i < 6; ++i) CHECK(c_counts[i] == cpp_counts[i]);
    }
}

TEST_CASE("the algorithm-selection constants are stable") {
    // These values are baked into the bindings; renumbering them silently would
    // make every binding pick the wrong algorithm.
    CHECK(ISOSPEC_ALGO_LAYERED == 0);
    CHECK(ISOSPEC_ALGO_ORDERED == 1);
    CHECK(ISOSPEC_ALGO_THRESHOLD_ABSOLUTE == 2);
    CHECK(ISOSPEC_ALGO_THRESHOLD_RELATIVE == 3);
    CHECK(ISOSPEC_ALGO_LAYERED_ESTIMATE == 4);
}
