// libFuzzer target: a formula as raw bytes, through Iso(std::string_view) -- the
// entry the Rust binding uses. Any std::exception is a correct rejection; a
// crash, a sanitizer report or a hang is a bug.

#include <cstddef>
#include <cstdint>
#include <exception>
#include <string_view>
#include <vector>

#include "isoSpec++.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    // Exactly the input bytes, NULs included and no terminator: an
    // exactly-sized heap copy, so that ASan reports any read past the end.
    const std::vector<char> formula(data, data + size);
    try
    {
        IsoSpec::Iso iso{std::string_view(formula.data(), formula.size())};
        volatile double sink = iso.getMonoisotopicPeakMass() + iso.getModeLProb() + iso.getTheoreticalAverageMass();
        (void)sink;
    }
    catch(const std::exception&)
    {
    }
    return 0;
}
