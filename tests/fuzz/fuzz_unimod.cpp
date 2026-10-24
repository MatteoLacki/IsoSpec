// libFuzzer target: an override Unimod CSV (what unimod_db_path points at),
// then a sequence parsed against the table it produced. The input is split at
// its first NUL: CSV text before it, sequence after it.

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <string>
#include <string_view>
#include <vector>

#include "composition.h"
#include "fasta_mods.h"
#include "isoSpec++.h"
#include "unimod.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    const char* text = reinterpret_cast<const char*>(data);
    const char* nul = static_cast<const char*>(std::memchr(text, '\0', size));
    const std::string csv(text, nul == nullptr ? size : static_cast<size_t>(nul - text));
    // The sequence part as exact bytes, as in fuzz_sequence.
    const std::vector<char> sequence = nul == nullptr ? std::vector<char>{'P', 'E', 'P', '[', 'U', 'N', 'I', 'M', 'O', 'D', ':', '1', ']', 'K'}
                                                      : std::vector<char>(nul + 1, text + size);
    try
    {
        const IsoSpec::UnimodTable table = IsoSpec::parse_unimod_csv(csv);
        for(unsigned int id = 0; id < 8; id++)
            if(const IsoSpec::UnimodEntry* entry = table.lookup(id))
            {
                volatile size_t sink = entry->name.size() + entry->element_first_index.size();
                (void)sink;
            }
        IsoSpec::ElementComposition composition = IsoSpec::parse_fasta_with_mods_full(std::string_view(sequence.data(), sequence.size()), table);
        IsoSpec::Iso iso = IsoSpec::build_iso_from_composition(composition);
        volatile double sink = iso.getMonoisotopicPeakMass();
        (void)sink;
    }
    catch(const std::exception&)
    {
    }
    return 0;
}
