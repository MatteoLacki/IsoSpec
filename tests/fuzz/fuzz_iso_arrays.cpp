// libFuzzer target: the raw-array Iso constructor (what setupIso and Python's
// Iso(atomCounts=..., isotopeMasses=..., isotopeProbabilities=...) reach),
// with arbitrary shapes, counts and IEEE doubles -- NaN, infinities, denormals
// -- then the envelope operations on whatever it accepts. Atom counts are kept
// small so that each input runs in milliseconds; the size limits themselves
// are plain integer checks, covered by the unit tests.

#include <fuzzer/FuzzedDataProvider.h>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <vector>

#include "fixedEnvelopes.h"
#include "isoSpec++.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    FuzzedDataProvider in(data, size);
    const int dim = in.ConsumeIntegralInRange<int>(-1, 4);
    std::vector<int> isotope_numbers, atom_counts;
    std::vector<double> masses, probs;
    for(int ii = 0; ii < (dim > 0 ? dim : 0); ii++)
    {
        const int n = in.ConsumeIntegralInRange<int>(-1, 6);
        isotope_numbers.push_back(n);
        atom_counts.push_back(in.ConsumeIntegralInRange<int>(-2, 40));
        for(int jj = 0; jj < n; jj++)
        {
            // Raw bit patterns half the time, so NaN, infinities and denormals
            // turn up; plausible values the other half, so inputs get past the
            // checks and into the algorithms.
            if(in.ConsumeBool())
            {
                masses.push_back(in.ConsumeFloatingPoint<double>());
                probs.push_back(in.ConsumeFloatingPoint<double>());
            }
            else
            {
                masses.push_back(in.ConsumeFloatingPointInRange<double>(-5.0, 300.0));
                probs.push_back(in.ConsumeFloatingPointInRange<double>(0.0, 1.0));
            }
        }
    }
    isotope_numbers.resize(4, 0);  // dim may exceed what was filled
    atom_counts.resize(4, 0);
    masses.resize(masses.size() + 8, 1.0);
    probs.resize(probs.size() + 8, 0.5);

    try
    {
        IsoSpec::Iso iso(dim, isotope_numbers.data(), atom_counts.data(), masses.data(), probs.data());
        volatile double sink = iso.getMonoisotopicPeakMass() + iso.getModeLProb() + iso.getLightestPeakMass() +
                               iso.getHeaviestPeakMass() + iso.getTheoreticalAverageMass() + iso.variance();
        (void)sink;

        IsoSpec::FixedEnvelope total = IsoSpec::FixedEnvelope::FromTotalProb(IsoSpec::Iso(iso, true), 0.99, true, true);
        IsoSpec::FixedEnvelope thresh = IsoSpec::FixedEnvelope::FromThreshold(IsoSpec::Iso(iso, true), 1e-3, false, true);
        std::mt19937 rng = IsoSpec::seeded_rng(size);
        IsoSpec::FixedEnvelope stoch = IsoSpec::FixedEnvelope::FromStochastic(IsoSpec::Iso(iso, true), 1000, 0.99, 5.0, true, rng);
        IsoSpec::FixedEnvelope binned = IsoSpec::FixedEnvelope::Binned(IsoSpec::Iso(iso, true), 0.99, 0.5);
        total.sort_by_mass();
        thresh.sort_by_prob();
        sink = total.empiric_average_mass() + thresh.get_total_prob() + stoch.get_total_prob() + binned.get_total_prob();
        if(total.confs_no() > 0 && thresh.confs_no() > 0)
        {
            IsoSpec::FixedEnvelope sum = total + thresh;
            total.normalize();
            thresh.normalize();
            sink = total.WassersteinDistance(thresh) + sum.get_total_prob();
            total.resample(500, 1.0, rng);
        }
    }
    catch(const std::exception&)
    {
    }
    return 0;
}
