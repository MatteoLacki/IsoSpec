// libFuzzer target: FixedEnvelope::Binned on raw-array molecules, with
// fuzzer-chosen bin width, middle and target probability -- IEEE doubles of
// every kind, and widths derived from the masses so that bin indices land on
// either side of the 2^48 limit, where the shifted grid origin rounds the most.
// The accumulator is mmap'ed, beyond ASan's sight: what catches a slot outside
// it is the bin_slot()/bin_run() check, live in this ISOSPEC_DEBUG build. The
// output must be finite, nonnegative, sum to at most 1, and sit on bins
// between the lightest and heaviest masses; a broken invariant aborts.
//
// Inputs whose work is legitimately large -- millions of bins, or (nearly) full
// enumerations of millions of configurations -- are skipped, so that each runs
// in milliseconds; the size checks themselves are unit-tested.

#include <fuzzer/FuzzedDataProvider.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <vector>

#include "fixedEnvelopes.h"
#include "isoSpec++.h"

namespace
{

// Raw bit patterns a third of the time, so that NaN, infinities, denormals and
// huge values turn up; plausible values otherwise.
double any_double(FuzzedDataProvider& in, double lo, double hi)
{
    return in.ConsumeIntegralInRange<int>(0, 2) == 0 ? in.ConsumeFloatingPoint<double>()
                                                     : in.ConsumeFloatingPointInRange<double>(lo, hi);
}

// Configurations of one element: C(atoms + isotopes - 1, isotopes - 1), capped.
double configurations(int atoms, int isotopes)
{
    double n = 1.0;
    for(int ii = 1; ii < isotopes; ii++)
        n = std::min(n * (atoms + ii) / ii, 1e300);
    return n;
}

}  // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    FuzzedDataProvider in(data, size);

    const int dim = in.ConsumeIntegralInRange<int>(0, 3);
    std::vector<int> isotope_numbers, atom_counts;
    std::vector<double> masses, probs;
    double enumeration = 1.0;
    for(int ii = 0; ii < dim; ii++)
    {
        const int n = in.ConsumeIntegralInRange<int>(1, 4);
        const int atoms = in.ConsumeIntegralInRange<int>(0, 40);
        isotope_numbers.push_back(n);
        atom_counts.push_back(atoms);
        enumeration *= configurations(atoms, n);
        // Masses either independent or a cluster of near-equal values (a few
        // ulps apart), the latter so that bins can be finer than their spacing.
        const double base = any_double(in, -5.0, 300.0);
        const bool cluster = in.ConsumeBool();
        double m = base;
        for(int jj = 0; jj < n; jj++)
        {
            if(cluster)
                for(int kk = in.ConsumeIntegralInRange<int>(1, 8); kk > 0; kk--)
                    m = std::nextafter(m, HUGE_VAL);
            else
                m = any_double(in, -5.0, 300.0);
            masses.push_back(m);
            probs.push_back(in.ConsumeFloatingPointInRange<double>(0.0, 1.0));
        }
    }
    masses.resize(masses.size() + 4, 1.0);  // dim may claim more than was filled
    probs.resize(probs.size() + 4, 0.5);

    try
    {
        IsoSpec::Iso iso(dim, isotope_numbers.data(), atom_counts.data(), masses.data(), probs.data());
        const double lightest = iso.getLightestPeakMass();
        const double heaviest = iso.getHeaviestPeakMass();
        const double scale = std::max({std::fabs(lightest), std::fabs(heaviest), 1e-300});

        double width;
        switch(in.ConsumeIntegralInRange<int>(0, 2))
        {
        case 0: width = any_double(in, 1e-6, 10.0); break;
        // Around the 2^48 index limit, from both sides.
        case 1: width = std::ldexp(scale, -in.ConsumeIntegralInRange<int>(40, 56)); break;
        default: width = (heaviest - lightest) / in.ConsumeIntegralInRange<int>(1, 100000); break;
        }
        double middle;
        switch(in.ConsumeIntegralInRange<int>(0, 2))
        {
        case 0: middle = any_double(in, -1.0, 1.0); break;
        case 1: middle = std::ldexp(in.ConsumeBool() ? 1.0 : -1.0, in.ConsumeIntegralInRange<int>(0, 1023)) +
                         in.ConsumeFloatingPointInRange<double>(-1.0, 1.0); break;
        default: middle = lightest + in.ConsumeFloatingPointInRange<double>(-2.0, 2.0) * width; break;
        }
        const double target = in.ConsumeBool() ? any_double(in, 0.0, 1.0) : 1.0;

        // Legitimate but slow: skip. (Past the 2^48 index limit Binned throws at
        // once, however many bins that would have been.)
        const double bins = (heaviest - lightest) / width;
        if(bins > 2e6 && scale / width < 0x1p49)
            return 0;
        // A target near 1 makes the layered generator walk nearly everything too.
        if(target > 0.99 && enumeration > 2e5)
            return 0;

        IsoSpec::FixedEnvelope env = IsoSpec::FixedEnvelope::Binned(std::move(iso), target, width, middle);

        const double slack = 3.0 * width + 1e-9 * scale;
        double total = 0.0;
        for(size_t ii = 0; ii < env.confs_no(); ii++)
        {
            const double m = env.mass(ii);
            const double p = env.prob(ii);
            if(!std::isfinite(m) || !std::isfinite(p) || !(p > 0.0) ||
               !(m >= lightest - slack && m <= heaviest + slack))
                std::abort();
            total += p;
        }
        if(!(total <= 1.0 + 1e-9))
            std::abort();
    }
    catch(const std::exception&)
    {
    }
    return 0;
}
