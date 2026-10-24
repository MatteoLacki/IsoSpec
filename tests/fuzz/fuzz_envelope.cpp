// libFuzzer target: FixedEnvelope::FromArrays on arbitrary IEEE doubles -- NaN,
// infinities, denormals, values near DBL_MAX -- then a fuzzer-chosen sequence of
// envelope operations with arbitrary scalar arguments. What FromArrays accepts
// must stay true through every operation, including one that throws (which must
// leave the envelope as it was): every mass finite, every probability finite
// and nonnegative. A broken invariant aborts, like a sanitizer report.

#include <fuzzer/FuzzedDataProvider.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <random>
#include <vector>

#include "fixedEnvelopes.h"
#include "isoSpec++.h"

namespace
{

// Raw bit patterns half the time, so that NaN, infinities, denormals and huge
// values turn up; plausible values the other half, so that inputs get past the
// checks and into the operations.
double any_double(FuzzedDataProvider& in, double lo, double hi)
{
    return in.ConsumeBool() ? in.ConsumeFloatingPoint<double>() : in.ConsumeFloatingPointInRange<double>(lo, hi);
}

IsoSpec::FixedEnvelope make_envelope(FuzzedDataProvider& in)
{
    const size_t n = in.ConsumeIntegralInRange<size_t>(0, 48);
    std::vector<double> masses, probs;
    for(size_t ii = 0; ii < n; ii++)
    {
        masses.push_back(any_double(in, -10.0, 2000.0));
        probs.push_back(any_double(in, 0.0, 1.0));
    }
    return IsoSpec::FixedEnvelope::FromArrays(masses.data(), probs.data(), n);
}

void check(const IsoSpec::FixedEnvelope& e)
{
    for(size_t ii = 0; ii < e.confs_no(); ii++)
        if(!std::isfinite(e.mass(ii)) || !(e.prob(ii) >= 0.0) || !std::isfinite(e.prob(ii)))
            std::abort();
}

}  // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    FuzzedDataProvider in(data, size);
    std::mt19937 rng = IsoSpec::seeded_rng(size);

    IsoSpec::FixedEnvelope a, b;
    try
    {
        a = make_envelope(in);
        b = make_envelope(in);
    }
    catch(const std::exception&)
    {
        return 0;
    }
    check(a);
    check(b);

    for(int step = 0; step < 16 && in.remaining_bytes() > 0; step++)
    {
        IsoSpec::FixedEnvelope& x = in.ConsumeBool() ? a : b;
        IsoSpec::FixedEnvelope& y = in.ConsumeBool() ? a : b;
        try
        {
            volatile double sink = 0.0;
            switch(in.ConsumeIntegralInRange<int>(0, 15))
            {
            case 0: x.scale(any_double(in, 0.0, 4.0)); break;
            case 1: x.shift_mass(any_double(in, -100.0, 100.0)); break;
            case 2: x.normalize(); break;
            case 3: x.resample(in.ConsumeIntegralInRange<size_t>(0, 5000), any_double(in, 0.0, 10.0), rng); break;
            case 4: x = x.bin(any_double(in, 0.0, 2.0), any_double(in, -1.0, 1.0)); break;
            case 5: x = x + y; break;
            case 6:
                if(x.confs_no() * y.confs_no() <= 4096)
                    x = x * y;
                break;
            case 7:
            {
                const IsoSpec::FixedEnvelope* spectra[2] = {&x, &y};
                const double intensities[2] = {any_double(in, 0.0, 3.0), any_double(in, 0.0, 3.0)};
                x = IsoSpec::FixedEnvelope::LinearCombination(spectra, intensities, 2);
                break;
            }
            case 8: sink = x.WassersteinDistance(y); break;
            case 9: sink = x.OrientedWassersteinDistance(y); break;
            case 10: sink = x.AbyssalWassersteinDistance(y, any_double(in, 0.0, 10.0), any_double(in, 0.0, 3.0)); break;
            case 11:
            {
                auto res = x.WassersteinMatch(y, any_double(in, 0.0, 10.0), any_double(in, 0.0, 3.0));
                sink = std::get<0>(res) + std::get<1>(res) + std::get<2>(res);
                break;
            }
            case 12: x.sort_by_mass(); break;
            case 13: x.sort_by_prob(); break;
            case 14: sink = x.empiric_average_mass() + x.empiric_variance() + x.get_total_prob(); break;
            case 15: x = IsoSpec::FixedEnvelope(x); break;
            }
            (void)sink;
        }
        catch(const std::exception&)
        {
        }
        check(a);
        check(b);
    }
    return 0;
}
