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
#include <cmath>
#include <cstdio>
#include <fstream>
#include <future>
#include <memory>
#include <stdexcept>
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

TEST_CASE("threads: the lazily filled log-factorial tail, filled by many threads at once") {
    // Atom counts past the eagerly filled prefix (ISOSPEC_LFACT_EAGER_FILL)
    // reach slots computed on first use. These counts appear in no other test,
    // and the reference is computed only *after* the threads, so the threads
    // are the first to touch those slots and really do race to fill them.
    static_assert(ISOSPEC_LFACT_EAGER_FILL < 70000, "these counts must land in the lazy tail");
    const char* const formulas[] = {"C71233H5", "C70411N3", "O72917H2"};

    auto work = [&](std::vector<double>& out) {
        for (const char* f : formulas) {
            out.push_back(Iso(f).getModeLProb());
            out.push_back(static_cast<double>(FixedEnvelope::FromThreshold(Iso(f), 1e-3, false).confs_no()));
        }
    };

    std::vector<std::vector<double>> got(kThreads);
    run_on_threads([&](unsigned t) { work(got[t]); });

    std::vector<double> want;
    work(want);
    for (unsigned t = 0; t < kThreads; t++) {
        INFO("thread " << t);
        CHECK(got[t] == want);
    }
}

TEST_CASE("threads: the C ABI end to end, one set of handles per thread") {
    // Everything a binding does: composition -> isotope arrays -> Iso ->
    // generator and envelopes -> zero-copy array hand-over -> distances.
    // Threshold envelopes throughout, since they are deterministic.
    const char* const symbols[] = {"C", "H", "N", "O", "S"};
    const int counts[] = {120, 190, 34, 37, 2};

    struct Result {
        double mono = 0.0;
        size_t gen_count = 0;
        size_t env_size = 0;
        double mass_sum = 0.0;
        double distance = 0.0;
        double fasta_mono = 0.0;
        bool ok = true;
    };

    auto work = [&]() {
        Result r;
        void* expanded = expandCompositionC(symbols, counts, 5, false);
        if (expanded == nullptr) { r.ok = false; return r; }
        auto make_iso = [&]() {
            return setupIso(5, expandedIsotopeNumbersC(expanded), expandedAtomCountsC(expanded),
                            expandedIsotopeMassesC(expanded), expandedIsotopeProbabilitiesC(expanded));
        };

        void* iso = make_iso();
        r.mono = getMonoisotopicPeakMassIso(iso);
        void* env = setupThresholdFixedEnvelope(iso, 1e-4, false, true);
        void* env2 = setupThresholdFixedEnvelope(iso, 1e-2, false, false);
        void* gen = setupIsoThresholdGenerator(iso, 1e-4, false, 1000, 1000, true);  // consumes iso
        while (advanceToNextConfigurationIsoThresholdGenerator(gen))
            r.gen_count++;
        deleteIsoThresholdGenerator(gen);
        deleteIso(iso);

        normalizeEnvelope(env);
        normalizeEnvelope(env2);
        r.distance = wassersteinDistance(env, env2);
        // `size` is the deleter's argument, not the element count.
        r.env_size = confs_noFixedEnvelope(env);
        size_t size = 0;
        IsoSpecArrayDeleter deleter = nullptr;
        double* masses = massesFixedEnvelopeWithDeleter(env, &size, &deleter);
        for (size_t ii = 0; ii < r.env_size; ii++)
            r.mass_sum += masses[ii];
        freeReleasedArrayWithDeleter(masses, size, deleter);
        deleteFixedEnvelope(env, false);
        deleteFixedEnvelope(env2, false);
        deleteExpandedCompositionC(expanded);

        void* fasta = isoFromFasta("MPEPTIDEKR", false, true);
        r.fasta_mono = getMonoisotopicPeakMassIso(fasta);
        deleteIso(fasta);
        return r;
    };

    const Result want = work();
    REQUIRE(want.ok);
    std::vector<std::vector<Result>> got(kThreads);
    run_on_threads([&](unsigned t) {
        for (int round = 0; round < 10; round++)
            got[t].push_back(work());
    });

    for (unsigned t = 0; t < kThreads; t++)
        for (const Result& r : got[t]) {
            INFO("thread " << t);
            CHECK(r.ok);
            CHECK(r.mono == want.mono);
            CHECK(r.gen_count == want.gen_count);
            CHECK(r.env_size == want.env_size);
            CHECK(r.mass_sum == want.mass_sum);
            CHECK(r.distance == want.distance);
            CHECK(r.fasta_mono == want.fasta_mono);
        }
}

TEST_CASE("threads: objects handed from thread to thread") {
    // No object may be tied to the thread that made it: a binding (Rust's
    // Send, a Python thread pool) builds on one thread and finishes on
    // another. Each stage below runs on its own thread, started only after the
    // previous one ended.
    const char* const formula = "C520H817N139O147S8";

    auto on_new_thread = [](auto&& fn) {
        std::thread th(std::forward<decltype(fn)>(fn));
        th.join();
    };

    std::unique_ptr<Iso> iso;
    on_new_thread([&] { iso = std::make_unique<Iso>(formula); });

    std::unique_ptr<IsoThresholdGenerator> gen;
    on_new_thread([&] { gen = std::make_unique<IsoThresholdGenerator>(std::move(*iso), 1e-4, false); });

    size_t count = 0;
    on_new_thread([&] { while (gen->advanceToNextConfiguration()) count++; });
    on_new_thread([&] { gen.reset(); });

    std::unique_ptr<FixedEnvelope> env;
    on_new_thread([&] { env = std::make_unique<FixedEnvelope>(FixedEnvelope::FromTotalProb(Iso(formula), 0.99, true)); });
    double total = 0.0;
    on_new_thread([&] { env->normalize(); env->sort_by_mass(); total = env->get_total_prob(); });
    on_new_thread([&] { env.reset(); });

    size_t want_count = 0;
    {
        IsoThresholdGenerator ref(Iso(formula), 1e-4, false);
        while (ref.advanceToNextConfiguration()) want_count++;
    }
    CHECK(count == want_count);
    CHECK(total == doctest::Approx(1.0));
}

TEST_CASE("threads: a stochastic generator keeps working after the thread that made it") {
    // The stochastic generator draws random numbers. Made on one thread, used
    // on another while the first keeps drawing for itself, and finished after
    // the first is gone: it must never touch the maker's random-number state.
    const char* const formula = "C520H817N139O147S8";
    const size_t molecules = 100000;

    for (int round = 0; round < 4; round++) {
        std::promise<std::unique_ptr<IsoStochasticGenerator>> handoff;
        auto received = handoff.get_future();
        std::atomic<bool> maker_may_stop{false};

        std::thread maker([&] {
            handoff.set_value(std::make_unique<IsoStochasticGenerator>(Iso(formula), molecules));
            // Keep this thread's own generator busy while the other thread
            // samples, so that a shared random-number state would be raced on.
            while (!maker_may_stop.load())
                FixedEnvelope::FromStochastic(Iso("C100H202"), 1000);
        });

        std::unique_ptr<IsoStochasticGenerator> gen = received.get();
        size_t seen = 0;
        std::thread user([&] {
            for (int ii = 0; ii < 50 && gen->advanceToNextConfiguration(); ii++)
                seen += gen->count();
        });
        user.join();
        maker_may_stop.store(true);
        maker.join();

        // The maker is gone now; finish the sample here.
        while (gen->advanceToNextConfiguration())
            seen += gen->count();
        CHECK(seen == molecules);
    }
}

TEST_CASE("threads: one shared FixedEnvelope, copied and linearly combined concurrently") {
    const FixedEnvelope shared = FixedEnvelope::FromThreshold(Iso("C200H300N50O60S2"), 1e-5, false);
    const FixedEnvelope other = FixedEnvelope::FromThreshold(Iso("C190H310N52O58S3"), 1e-5, false);

    // Per-thread copies are free to be mutated: sorted, normalised,
    // resampled, measured, binned.
    auto work = [&](std::vector<double>& out) {
        FixedEnvelope mine(shared);
        FixedEnvelope theirs(other);
        out.push_back(mine.WassersteinDistance(theirs));
        mine.normalize();
        out.push_back(mine.empiric_average_mass());
        out.push_back(static_cast<double>(mine.bin(1.0).confs_no()));
        const FixedEnvelope* parts[] = {&shared, &other};
        const double weights[] = {0.25, 0.75};
        out.push_back(FixedEnvelope::LinearCombination(parts, weights, 2).get_total_prob());
        FixedEnvelope sampled(shared);
        sampled.resample(1000);
        out.push_back(sampled.get_total_prob());
    };

    std::vector<double> want;
    work(want);
    std::vector<std::vector<double>> got(kThreads);
    run_on_threads([&](unsigned t) { work(got[t]); });

    for (unsigned t = 0; t < kThreads; t++) {
        INFO("thread " << t);
        REQUIRE(got[t].size() == want.size());
        CHECK(got[t][0] == want[0]);
        CHECK(got[t][1] == want[1]);
        CHECK(got[t][2] == want[2]);
        CHECK(got[t][3] == doctest::Approx(want[3]));
        CHECK(got[t][4] == doctest::Approx(1000.0));  // resample conserves the ion count
    }
}

TEST_CASE("threads: errors raised on many threads at once") {
    // Bad input from every thread at the same moment: exceptions in C++, NULL
    // or NaN through the C ABI, and the Unimod cache refusing a missing file
    // (which must not leave its lock held).
    std::vector<int> failures_seen(kThreads, 0);
    run_on_threads([&](unsigned t) {
        int seen = 0;
        for (int round = 0; round < 20; round++) {
            try { Iso iso("C1C1"); } catch (const std::invalid_argument&) { seen++; }
            try { Iso iso("Xx5"); } catch (const std::invalid_argument&) { seen++; }
            try {
                Iso::FromFASTAWithMods("PEPC[UNIMOD:4]", false, true,
                                       ("build/no_such_unimod_" + std::to_string(t % 2) + ".csv").c_str());
            } catch (const std::invalid_argument&) { seen++; }
            if (parseFastaWithModsC("PEP[UNIMOD:999999]TIDE", nullptr) == nullptr) seen++;
            if (parseFastaWithModsC("PEP[UNIMOD:4", nullptr) == nullptr) seen++;
            const char* bad_symbol[] = {"Qq"};
            const int one[] = {1};
            if (expandCompositionC(bad_symbol, one, 1, false) == nullptr) seen++;
        }
        failures_seen[t] = seen;
    });

    for (unsigned t = 0; t < kThreads; t++)
        CHECK(failures_seen[t] == 20 * 6);
    // And the cache still works afterwards, on this thread.
    CHECK(Iso::FromFASTAWithMods("PEPC[UNIMOD:4]").getMonoisotopicPeakMass() > 0.0);
}

TEST_CASE("threads: many short-lived threads") {
    // Thread-local state (the random-number generator) is created and torn
    // down with every thread; churn through many of them.
    const char* const formula = "C100H202";
    const size_t want = FixedEnvelope::FromThreshold(Iso(formula), 1e-3, false).confs_no();
    std::atomic<int> wrong{0};
    for (int batch = 0; batch < 8; batch++) {
        std::vector<std::thread> threads;
        for (int ii = 0; ii < 16; ii++)
            threads.emplace_back([&] {
                FixedEnvelope sampled = FixedEnvelope::FromStochastic(Iso(formula), 500);
                if (std::lround(sampled.get_total_prob()) != 500) wrong++;
                if (FixedEnvelope::FromThreshold(Iso(formula), 1e-3, false).confs_no() != want) wrong++;
            });
        for (auto& th : threads)
            th.join();
    }
    CHECK(wrong.load() == 0);
}
