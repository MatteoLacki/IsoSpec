#include "isospec/rust/src/ffi.rs.h"

#include <string_view>

namespace isospec_rs {

namespace {

std::string_view view(rust::Str s) { return std::string_view(s.data(), s.size()); }

}  // namespace

std::unique_ptr<Iso> iso_from_formula(rust::Str formula, bool use_nominal_masses)
{
    return std::make_unique<Iso>(view(formula), use_nominal_masses);
}

rust::Str simd_level() { return rust::Str(IsoSpec::active_simd_level()); }

}  // namespace isospec_rs
