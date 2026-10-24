// libFuzzer target: a peptide sequence, plain and with [UNIMOD:<id>]
// modifications resolved against the embedded table, through the paths
// Iso::FromFASTA, Iso::FromFASTAWithMods and the composition parser take.

#include <cstddef>
#include <cstdint>
#include <exception>
#include <string>

#include "composition.h"
#include "fasta_mods.h"
#include "isoSpec++.h"
#include "unimod.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    const std::string sequence(reinterpret_cast<const char*>(data), size);
    const bool add_water = size % 2 == 0;
    try
    {
        IsoSpec::Iso plain = IsoSpec::Iso::FromFASTA(sequence.c_str(), false, add_water);
        volatile double sink = plain.getMonoisotopicPeakMass();
        (void)sink;
    }
    catch(const std::exception&)
    {
    }
    try
    {
        IsoSpec::ElementComposition composition =
            IsoSpec::parse_fasta_with_mods(sequence.c_str(), IsoSpec::embedded_unimod_table());
        IsoSpec::Iso modified = IsoSpec::build_iso_from_composition(composition);
        volatile double sink = modified.getMonoisotopicPeakMass() + modified.getModeLProb();
        (void)sink;
    }
    catch(const std::exception&)
    {
    }
    return 0;
}
