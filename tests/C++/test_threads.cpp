// Thread safety.
//
// The library does no parallel work of its own, but it must stand being used
// from several threads at once: independent objects on separate threads, and
// const access to one shared Iso or FixedEnvelope from many. Each test runs the
// same work on several threads and compares every thread's answer with a
// single-threaded reference, so that a broken build fails here in any
// configuration; under `make tsan` the same tests also prove there is no data
// race on the way to those answers.
//
// Assertions are made on the main thread only, after joining, rather than
// trusting the test framework's own thread safety.

#include <atomic>
#include <cstdio>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"
#include "cwrapper.h"
#include "fasta_mods.h"
#include "isa_kernels.h"

using namespace IsoSpec;

namespace {

const unsigned kThreads = 8;

const char* const kFormulas[] = {"C100H202", "C2H6O1", "C520H817N139O147S8", "Se2Fe3C40H60", "H2O1"};

// Everything one thread computes for one formula, through every generator and
// every envelope constructor.
struct Summary {
    size_t threshold_count = 0;
    double threshold_mass_sum = 0.0;
    size_t layered_count = 0;
    double ordered_first_lprobs[16] = {};
    size_t total_prob_count = 0;
    double total_prob_prob = 0.0;
    size_t binned_count = 0;
    double binned_prob = 0.0;
    double mode_lprob = 0.0;
    double mono_mass = 0.0;
    size_t stochastic_molecules = 0;
};

Summary summarise(const char* formula) {
    Summary s;
    {
        Iso iso(formula);
        s.mode_lprob = iso.getModeLProb();
        s.mono_mass = iso.getMonoisotopicPeakMass();
    }
    {
        IsoThresholdGenerator gen(Iso(formula), 1e-4, false);
        while (gen.advanceToNextConfiguration()) {
            s.threshold_count++;
            s.threshold_mass_sum += gen.mass();
        }
    }
    {
        IsoLayeredGenerator gen{Iso(formula)};
        double covered = 0.0;
        while (covered < 0.99 && gen.advanceToNextConfiguration()) {
            s.layered_count++;
            covered += gen.prob();
        }
    }
    {
        IsoOrderedGenerator gen{Iso(formula)};
        for (double& lp : s.ordered_first_lprobs)
            lp = gen.advanceToNextConfiguration() ? gen.lprob() : 0.0;
    }
    {
        FixedEnvelope env = FixedEnvelope::FromTotalProb(Iso(formula), 0.999, true, true);
        s.total_prob_count = env.confs_no();
        s.total_prob_prob = env.get_total_prob();
    }
    {
        FixedEnvelope env = FixedEnvelope::Binned(Iso(formula), 0.999, 1.0);
        s.binned_count = env.confs_no();
        s.binned_prob = env.get_total_prob();
    }
    {
        // Random, so only its bookkeeping is comparable: every molecule asked
        // for is accounted for exactly once.
        IsoStochasticGenerator gen(Iso(formula), 10000);
        while (gen.advanceToNextConfiguration())
            s.stochastic_molecules += static_cast<size_t>(gen.count());
    }
    return s;
}

void check_same(const Summary& got, const Summary& want) {
    CHECK(got.threshold_count == want.threshold_count);
    CHECK(got.threshold_mass_sum == want.threshold_mass_sum);
    CHECK(got.layered_count == want.layered_count);
    for (size_t ii = 0; ii < 16; ii++)
        CHECK(got.ordered_first_lprobs[ii] == want.ordered_first_lprobs[ii]);
    // FromTotalProb's trim picks quickselect pivots at random, so with tied
    // probabilities at the cut the kept set may differ -- never its size.
    CHECK(got.total_prob_count == want.total_prob_count);
    CHECK(got.total_prob_prob == doctest::Approx(want.total_prob_prob).epsilon(1e-12));
    CHECK(got.binned_count == want.binned_count);
    CHECK(got.binned_prob == doctest::Approx(want.binned_prob).epsilon(1e-12));
    CHECK(got.mode_lprob == want.mode_lprob);
    CHECK(got.mono_mass == want.mono_mass);
    CHECK(got.stochastic_molecules == want.stochastic_molecules);
}

template<typename F>
void run_on_threads(F&& body) {
    std::vector<std::thread> threads;
    for (unsigned t = 0; t < kThreads; t++)
        threads.emplace_back(body, t);
    for (auto& th : threads)
        th.join();
}

}  // namespace

TEST_CASE("threads: independent Isos through every generator and envelope") {
    std::vector<Summary> want;
    for (const char* f : kFormulas)
        want.push_back(summarise(f));

    std::vector<std::vector<Summary>> got(kThreads);
    run_on_threads([&](unsigned t) {
        // Each thread walks the formulas from a different starting point, so
        // that different molecules really are being built at the same moment.
        const size_t n = sizeof(kFormulas) / sizeof(kFormulas[0]);
        got[t].resize(n);
        for (size_t ii = 0; ii < n; ii++) {
            const size_t idx = (ii + t) % n;
            got[t][idx] = summarise(kFormulas[idx]);
        }
    });

    for (unsigned t = 0; t < kThreads; t++)
        for (size_t ii = 0; ii < want.size(); ii++) {
            INFO("thread " << t << ", formula " << kFormulas[ii]);
            check_same(got[t][ii], want[ii]);
        }
}

TEST_CASE("threads: one shared Iso, read and copied concurrently") {
    // Freshly built and never queried before the threads start: anything the
    // const accessors set up lazily would be set up by the threads, racing.
    const Iso shared("C520H817N139O147S8");

    std::vector<double> mode(kThreads), mode_mass(kThreads), mono(kThreads), avg(kThreads), var(kThreads);
    std::vector<size_t> env_count(kThreads), binned_count(kThreads);
    run_on_threads([&](unsigned t) {
        mode[t] = shared.getModeLProb();
        mode_mass[t] = shared.getModeMass();
        mono[t] = shared.getMonoisotopicPeakMass();
        avg[t] = shared.getTheoreticalAverageMass();
        var[t] = shared.variance();
        // The const-Iso overloads copy it, so this reads the shared one.
        env_count[t] = FixedEnvelope::FromTotalProb(shared, 0.99, false).confs_no();
        binned_count[t] = FixedEnvelope::Binned(shared, 0.99, 1.0).confs_no();
    });

    for (unsigned t = 1; t < kThreads; t++) {
        INFO("thread " << t);
        CHECK(mode[t] == mode[0]);
        CHECK(mode_mass[t] == mode_mass[0]);
        CHECK(mono[t] == mono[0]);
        CHECK(avg[t] == avg[0]);
        CHECK(var[t] == var[0]);
        CHECK(env_count[t] == env_count[0]);
        CHECK(binned_count[t] == binned_count[0]);
    }
    CHECK(mode[0] == Iso("C520H817N139O147S8").getModeLProb());
}

TEST_CASE("threads: one shared FixedEnvelope, combined concurrently") {
    const FixedEnvelope a = FixedEnvelope::FromThreshold(Iso("C20H30O5"), 1e-6, false);
    const FixedEnvelope b = FixedEnvelope::FromThreshold(Iso("H2O1"), 1e-6, false);
    const size_t want_sum = (a + b).confs_no();
    const size_t want_conv = (a * b).confs_no();

    std::vector<size_t> sum(kThreads), conv(kThreads);
    run_on_threads([&](unsigned t) {
        sum[t] = (a + b).confs_no();
        conv[t] = (a * b).confs_no();
    });

    for (unsigned t = 0; t < kThreads; t++) {
        CHECK(sum[t] == want_sum);
        CHECK(conv[t] == want_conv);
    }
}

TEST_CASE("threads: Unimod-modified peptides, embedded and override tables at once") {
    // Two distinct override files, so the path cache is filled and read from
    // several threads with more than one key -- the case a single-entry cache
    // would get wrong.
    const std::string path_a = "build/test_threads_unimod_a.csv";
    const std::string path_b = "build/test_threads_unimod_b.csv";
    for (const std::string& path : {path_a, path_b}) {
        std::ofstream f(path);
        f << "id,name,mono_mass,composition\n";
        // Both list Carbamidomethyl, but with different compositions, so a
        // thread handed the wrong table gets a different answer.
        f << (path == path_a ? "4,Carbamidomethyl,57.021464,H3C2N1O1\n" : "4,Heavier,58.0,H4C2N1O1\n");
    }

    const char* const sequence = "MPEPTC[UNIMOD:4]DEKC[UNIMOD:4]";
    const char* const paths[] = {nullptr, path_a.c_str(), path_b.c_str()};
    // Hydrogen count as the C entry point reports it, per table.
    auto hydrogens = [&](void* comp) {
        const char* const* symbols = compositionSymbolsC(comp);
        const int* counts = compositionCountsC(comp);
        for (size_t ii = 0; ii < compositionSizeC(comp); ii++)
            if (std::string(symbols[ii]) == "H")
                return counts[ii];
        return -1;
    };

    double want[3];
    int want_h[3];
    for (int ii = 0; ii < 3; ii++) {
        want[ii] = Iso::FromFASTAWithMods(sequence, false, true, paths[ii]).getMonoisotopicPeakMass();
        void* comp = parseFastaWithModsC(sequence, paths[ii]);
        REQUIRE(comp != nullptr);
        want_h[ii] = hydrogens(comp);
        deleteCompositionC(comp);
    }
    CHECK(want[1] != want[2]);
    CHECK(want_h[2] == want_h[1] + 2);

    const int kRounds = 50;
    std::vector<std::vector<double>> got(kThreads);
    std::vector<int> c_api_ok(kThreads, 1);
    run_on_threads([&](unsigned t) {
        for (int round = 0; round < kRounds; round++) {
            const int which = static_cast<int>((round + t) % 3);
            got[t].push_back(Iso::FromFASTAWithMods(sequence, false, true, paths[which]).getMonoisotopicPeakMass());

            // The C entry point, which used to share one static scratch buffer.
            void* comp = parseFastaWithModsC(sequence, paths[which]);
            if (comp == nullptr) {
                c_api_ok[t] = 0;
                continue;
            }
            if (hydrogens(comp) != want_h[which])
                c_api_ok[t] = 0;
            deleteCompositionC(comp);
        }
    });

    for (unsigned t = 0; t < kThreads; t++) {
        INFO("thread " << t);
        CHECK(c_api_ok[t] == 1);
        for (int round = 0; round < kRounds; round++)
            CHECK(got[t][round] == want[(round + t) % 3]);
    }

    std::remove(path_a.c_str());
    std::remove(path_b.c_str());
}

TEST_CASE("threads: ISA level switched while envelopes are being computed") {
    const IsaLevel saved = current_isa_level();
    const IsaLevel best = max_supported_isa_level();

    const char* const formula = "C520H817N139O147S8";
    const size_t want = FixedEnvelope::Binned(Iso(formula), 0.999, 1.0).confs_no();

    std::atomic<bool> stop{false};
    std::thread switcher([&] {
        bool flip = false;
        while (!stop.load()) {
            set_isa_level(flip ? best : ISA_LEVEL_BASELINE);
            flip = !flip;
        }
    });

    std::vector<std::vector<size_t>> got(kThreads);
    run_on_threads([&](unsigned t) {
        for (int round = 0; round < 5; round++)
            got[t].push_back(FixedEnvelope::Binned(Iso(formula), 0.999, 1.0).confs_no());
    });
    stop.store(true);
    switcher.join();
    set_isa_level(saved);

    for (unsigned t = 0; t < kThreads; t++)
        for (size_t n : got[t])
            CHECK(n == want);
}
