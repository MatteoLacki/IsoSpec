// The C++ side of the cxx bridge (rust/src/ffi.rs): free functions wrapping
// what cxx cannot bind directly -- constructors, overloads, default
// arguments, std::mt19937 -- in the library's own namespace-free shapes.
// Nothing here holds state or implements anything the library does not.

#pragma once

#include <exception>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>

#include "rust/cxx.h"

#include "isoSpec++.h"
#include "isa_kernels.h"

// How a C++ exception crosses into Rust. cxx's default keeps only what(); this
// prefixes it with a one-letter kind, which rust/src/error.rs turns back into
// ErrorKind. Most derived classes first. An exception escaping a bridge
// function not declared to return Result calls std::terminate instead, so
// every function that can throw is declared with Result.
namespace rust {
namespace behavior {

template <typename Try, typename Fail>
static void trycatch(Try&& func, Fail&& fail) noexcept try {
    func();
} catch (const std::bad_alloc& e) {
    fail((std::string("M") + e.what()).c_str());
} catch (const std::invalid_argument& e) {
    fail((std::string("I") + e.what()).c_str());
} catch (const std::domain_error& e) {
    fail((std::string("I") + e.what()).c_str());
} catch (const std::length_error& e) {
    fail((std::string("L") + e.what()).c_str());
} catch (const std::out_of_range& e) {
    fail((std::string("R") + e.what()).c_str());
} catch (const std::overflow_error& e) {
    fail((std::string("O") + e.what()).c_str());
} catch (const std::exception& e) {
    fail((std::string("X") + e.what()).c_str());
} catch (...) {
    fail("Xunknown C++ exception");
}

}  // namespace behavior
}  // namespace rust

namespace isospec_rs {

using IsoSpec::Iso;

std::unique_ptr<Iso> iso_from_formula(rust::Str formula, bool use_nominal_masses);

rust::Str simd_level();

}  // namespace isospec_rs
