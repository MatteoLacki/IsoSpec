// Reproducible randomness: seeded_rng() and the rng parameters of the sampling
// API (FromStochastic, resample, the stochastic generator), and the
// deterministic pivots that make FromTotalProb's output the same on every run.

#include <cstdint>
#include <random>
#include <thread>
#include <utility>
#include <vector>

#include "doctest.h"
#include "test_helpers.h"

using namespace IsoSpec;

namespace {

std::pair<std::vector<double>, std::vector<double>> contents(const FixedEnvelope& env)
{
    return {std::vector<double>(env.masses(), env.masses() + env.confs_no()),
            std::vector<double>(env.probs(), env.probs() + env.confs_no())};
}

template<typename F> auto on_another_thread(F f)
{
    decltype(f()) result;
    std::thread t([&] { result = f(); });
    t.join();
    return result;
}

}  // namespace

TEST_CASE("seeded_rng: the stream is fixed by the whole 64-bit seed") {
    // std::seed_seq and std::mt19937 are fully specified by the standard, so
    // these values hold on every platform and standard library.
    CHECK(seeded_rng(0)() == 2579818388u);
    CHECK(seeded_rng(0x123456789abcdefULL)() == 3070560217u);

    std::mt19937 a = seeded_rng(12345);
    std::mt19937 b = seeded_rng(12345);
    for (int i = 0; i < 100; i++)
        CHECK(a() == b());
    // The high half of the seed counts too.
    CHECK(seeded_rng(12345)() != seeded_rng(12345 + (std::uint64_t{1} << 32))());
}

TEST_CASE("FromStochastic: same seed, same sample, whichever thread draws it") {
    auto sample = [](std::uint64_t seed) {
        std::mt19937 rng = seeded_rng(seed);
        return contents(FixedEnvelope::FromStochastic(Iso("C520H817N139O147S8"), 100000, 0.9999, 5.0, false, rng));
    };
    const auto first = sample(7);
    CHECK(sample(7) == first);
    CHECK(on_another_thread([&] { return sample(7); }) == first);
    CHECK(sample(8) != first);
}

TEST_CASE("the stochastic generator draws from the rng it is given") {
    auto draw = [](std::uint64_t seed) {
        std::mt19937 rng = seeded_rng(seed);
        IsoStochasticGenerator gen(Iso("C100H202"), 50000, 0.9999, 5.0, rng);
        std::vector<std::pair<double, double>> out;
        while (gen.advanceToNextConfiguration())
            out.emplace_back(gen.mass(), gen.prob());
        return out;
    };
    const auto first = draw(11);
    CHECK(draw(11) == first);
    CHECK(draw(12) != first);
}

TEST_CASE("resample: same seed, same counts") {
    auto counts = [](std::uint64_t seed) {
        FixedEnvelope env = FixedEnvelope::FromTotalProb(Iso("C100H200N40O40"), 0.999, true, false);
        env.normalize();
        std::mt19937 rng = seeded_rng(seed);
        env.resample(100000, 1.0, rng);
        return contents(env);
    };
    const auto first = counts(3);
    CHECK(counts(3) == first);
    CHECK(counts(4) != first);
}

TEST_CASE("FromTotalProb gives the same envelope, in the same order, every time") {
    // Its trimming quickselect used to pick pivots from the randomly seeded
    // random_gen, so the order of the result differed between runs (and
    // between calls in one run). Both an optimized and a plain cut, on a
    // molecule with plenty of equally probable configurations.
    for (bool optimize : {true, false}) {
        CAPTURE(optimize);
        auto run = [&] { return contents(FixedEnvelope::FromTotalProb(Iso("C200H400O20"), 0.9999, optimize, false)); };
        const auto first = run();
        CHECK(run() == first);
        CHECK(on_another_thread(run) == first);
    }
}
