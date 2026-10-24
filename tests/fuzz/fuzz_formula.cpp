// libFuzzer target: a formula string, as Iso(const char*) and the C ABI's
// callers receive it. Any std::exception is a correct rejection; a crash, a
// sanitizer report or a hang is a bug.

#include <cstddef>
#include <cstdint>
#include <exception>
#include <string>

#include "isoSpec++.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    // A C string: the parser sees everything up to the first NUL.
    const std::string formula(reinterpret_cast<const char*>(data), size);
    try
    {
        IsoSpec::Iso iso(formula.c_str());
        volatile double sink = iso.getMonoisotopicPeakMass() + iso.getModeLProb() + iso.getTheoreticalAverageMass();
        (void)sink;
    }
    catch(const std::exception&)
    {
    }
    return 0;
}
