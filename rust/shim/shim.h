// The C++ side of the cxx bridge (rust/src/ffi.rs): thin free functions and
// wrapper classes around what cxx cannot bind directly -- constructors,
// overloads, default arguments, class templates, std::mt19937, arrays of
// pointers. Nothing here implements anything the library does not; the
// argument checks the C++ cannot make itself (that a slice is as long as a
// count says) are made on the Rust side before any call lands here.

#pragma once

#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <new>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "rust/cxx.h"

#include "composition.h"
#include "fixedEnvelopes.h"
#include "isa_kernels.h"
#include "isoMath.h"
#include "isoSpec++.h"
#include "unimod.h"

// How a C++ exception crosses into Rust. cxx's default keeps only what(); this
// prefixes it with a one-letter kind, which rust/src/error.rs turns back into
// ErrorKind. Most derived classes first. An exception escaping a bridge
// function not declared to return Result would call std::terminate, so every
// function that can throw is declared with Result.
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

using IsoSpec::FixedEnvelope;
using IsoSpec::Iso;
using IsoSpec::UnimodTable;

struct PeakData;
struct MatchResult;

// ---- Iso ------------------------------------------------------------------

std::unique_ptr<Iso> iso_from_formula(rust::Str formula, bool use_nominal_masses);
std::unique_ptr<Iso> iso_from_fasta(rust::Str sequence, bool use_nominal_masses, bool add_water);
std::unique_ptr<Iso> iso_from_peptide(rust::Str sequence, const UnimodTable& mods, bool use_nominal_masses,
                                      bool add_water);
std::unique_ptr<Iso> iso_from_composition(rust::Slice<const int32_t> element_first_index,
                                          rust::Slice<const int32_t> counts, bool use_nominal_masses);
std::unique_ptr<Iso> iso_from_isotopes(rust::Slice<const int32_t> isotope_numbers,
                                       rust::Slice<const int32_t> atom_counts, rust::Slice<const double> masses,
                                       rust::Slice<const double> probs);
std::unique_ptr<Iso> iso_clone(const Iso& iso);
void iso_add_element(Iso& iso, int32_t atom_count, rust::Slice<const double> masses,
                     rust::Slice<const double> probs);

rust::Slice<const int32_t> iso_isotope_numbers(const Iso& iso);
rust::Slice<const int32_t> iso_atom_counts(const Iso& iso);
rust::Slice<const double> iso_element_masses(const Iso& iso, int32_t element);
rust::Slice<const double> iso_element_lprobs(const Iso& iso, int32_t element);
void iso_lightest_conf(const Iso& iso, rust::Slice<int32_t> out);
void iso_heaviest_conf(const Iso& iso, rust::Slice<int32_t> out);
void iso_monoisotopic_conf(const Iso& iso, rust::Slice<int32_t> out);
void iso_marginal_log_size_estimates(const Iso& iso, rust::Slice<double> out, double target_total_prob);

// ---- Compositions and the element table ------------------------------------

void composition_from_formula(rust::Str formula, rust::Vec<int32_t>& element_first_index,
                              rust::Vec<int32_t>& counts);
void composition_from_peptide(rust::Str sequence, const UnimodTable& mods, bool add_water,
                              rust::Vec<int32_t>& element_first_index, rust::Vec<int32_t>& counts);
void composition_from_fasta(rust::Str sequence, bool add_water, rust::Slice<int32_t> chnosse);
void composition_expand(rust::Slice<const int32_t> element_first_index, rust::Slice<const int32_t> counts,
                        bool use_nominal_masses, rust::Vec<int32_t>& isotope_numbers,
                        rust::Vec<int32_t>& atom_counts, rust::Vec<double>& masses, rust::Vec<double>& probs);

int32_t element_find(rust::Str symbol);
int32_t element_isotope_count(int32_t first_index);
size_t element_table_len();
rust::Str element_table_symbol(size_t row);
rust::Str element_table_name(size_t row);
int32_t element_table_atomic_number(size_t row);
double element_table_mass(size_t row);
double element_table_probability(size_t row);
double element_table_mass_number(size_t row);

// ---- Unimod ----------------------------------------------------------------

const UnimodTable& unimod_embedded();
const UnimodTable& unimod_for_path(rust::Str path);
std::unique_ptr<UnimodTable> unimod_from_csv(rust::Str csv);
size_t unimod_size(const UnimodTable& table);
bool unimod_lookup(const UnimodTable& table, uint32_t id, rust::String& name, double& mono_mass,
                   rust::Vec<int32_t>& element_first_index, rust::Vec<int32_t>& counts);

// ---- Generators ------------------------------------------------------------
//
// Wrapper classes rather than the library's own: the generators are class
// templates behind `using` aliases, and the stochastic one keeps a reference
// to the std::mt19937 it is given, which must therefore live next to it.

struct ThresholdGenerator {
    IsoSpec::IsoThresholdGenerator g;
    ThresholdGenerator(Iso&& iso, double threshold, bool absolute, bool reorder_marginals)
        : g(std::move(iso), threshold, absolute, IsoSpec::default_tab_size, 1000, reorder_marginals) {}
};

struct LayeredGenerator {
    IsoSpec::IsoLayeredGenerator g;
    LayeredGenerator(Iso&& iso, bool reorder_marginals, double t_prob_hint)
        : g(std::move(iso), IsoSpec::default_tab_size, 1000, reorder_marginals, t_prob_hint) {}
};

struct OrderedGenerator {
    IsoSpec::IsoOrderedGenerator g;
    explicit OrderedGenerator(Iso&& iso) : g(std::move(iso)) {}
};

struct StochasticGenerator {
    // Declared before g, so constructed before and destroyed after it.
    std::unique_ptr<std::mt19937> rng;
    IsoSpec::IsoStochasticGenerator g;
    StochasticGenerator(Iso&& iso, size_t no_molecules, double precision, double beta_bias,
                        std::unique_ptr<std::mt19937> seeded);
};

std::unique_ptr<ThresholdGenerator> threshold_generator_new(std::unique_ptr<Iso> iso, double threshold,
                                                            bool absolute, bool reorder_marginals);
bool threshold_generator_next(ThresholdGenerator& gen, PeakData& out);
size_t threshold_generator_fill(ThresholdGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs);
void threshold_generator_conf(const ThresholdGenerator& gen, rust::Slice<int32_t> out);
size_t threshold_generator_count_confs(ThresholdGenerator& gen);
void threshold_generator_reset(ThresholdGenerator& gen);

std::unique_ptr<LayeredGenerator> layered_generator_new(std::unique_ptr<Iso> iso, bool reorder_marginals,
                                                        double t_prob_hint);
bool layered_generator_next(LayeredGenerator& gen, PeakData& out);
size_t layered_generator_fill(LayeredGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs);
void layered_generator_conf(const LayeredGenerator& gen, rust::Slice<int32_t> out);

std::unique_ptr<OrderedGenerator> ordered_generator_new(std::unique_ptr<Iso> iso);
bool ordered_generator_next(OrderedGenerator& gen, PeakData& out);
size_t ordered_generator_fill(OrderedGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs);
void ordered_generator_conf(const OrderedGenerator& gen, rust::Slice<int32_t> out);

std::unique_ptr<StochasticGenerator> stochastic_generator_new(std::unique_ptr<Iso> iso, size_t no_molecules,
                                                              double precision, double beta_bias, bool seeded,
                                                              uint64_t seed);
bool stochastic_generator_next(StochasticGenerator& gen, PeakData& out);
size_t stochastic_generator_fill(StochasticGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs);
void stochastic_generator_conf(const StochasticGenerator& gen, rust::Slice<int32_t> out);

// ---- FixedEnvelope -----------------------------------------------------------

std::unique_ptr<FixedEnvelope> envelope_threshold(const Iso& iso, double threshold, bool absolute, bool confs);
std::unique_ptr<FixedEnvelope> envelope_total_prob(const Iso& iso, double target, bool optimize, bool confs);
std::unique_ptr<FixedEnvelope> envelope_stochastic(const Iso& iso, size_t no_molecules, double precision,
                                                   double beta_bias, bool confs, bool seeded, uint64_t seed);
std::unique_ptr<FixedEnvelope> envelope_binned(const Iso& iso, double target, double bin_width,
                                               double bin_middle);
std::unique_ptr<FixedEnvelope> envelope_from_arrays(rust::Slice<const double> masses,
                                                    rust::Slice<const double> probs);
std::unique_ptr<FixedEnvelope> envelope_empty();
std::unique_ptr<FixedEnvelope> envelope_clone(const FixedEnvelope& env);
std::unique_ptr<FixedEnvelope> envelope_add(const FixedEnvelope& a, const FixedEnvelope& b);
std::unique_ptr<FixedEnvelope> envelope_convolve(const FixedEnvelope& a, const FixedEnvelope& b);
std::unique_ptr<FixedEnvelope> envelope_bin(FixedEnvelope& env, double bin_width, double middle);

//! LinearCombination takes arrays of pointers, which cxx cannot pass: the
//! Rust side fills one of these instead.
struct EnvelopeList {
    std::vector<const FixedEnvelope*> spectra;
    std::vector<double> intensities;
};
std::unique_ptr<EnvelopeList> envelope_list_new();
void envelope_list_push(EnvelopeList& list, const FixedEnvelope& env, double intensity);
std::unique_ptr<FixedEnvelope> envelope_linear_combination(const EnvelopeList& list);

size_t envelope_len(const FixedEnvelope& env);
int32_t envelope_all_dim(const FixedEnvelope& env);
rust::Slice<const double> envelope_masses(const FixedEnvelope& env);
rust::Slice<const double> envelope_probs(const FixedEnvelope& env);
rust::Slice<const int32_t> envelope_confs(const FixedEnvelope& env);

void envelope_sort_by_mass(FixedEnvelope& env);
void envelope_sort_by_prob(FixedEnvelope& env);
double envelope_total_prob_of(FixedEnvelope& env);
void envelope_scale(FixedEnvelope& env, double factor);
void envelope_normalize(FixedEnvelope& env);
void envelope_shift_mass(FixedEnvelope& env, double shift);
void envelope_resample(FixedEnvelope& env, size_t ionic_current, double beta_bias, bool seeded, uint64_t seed);
double envelope_average_mass(FixedEnvelope& env);
double envelope_variance(FixedEnvelope& env);
double envelope_wasserstein(FixedEnvelope& a, FixedEnvelope& b);
double envelope_oriented_wasserstein(FixedEnvelope& a, FixedEnvelope& b);
double envelope_abyssal_wasserstein(FixedEnvelope& a, FixedEnvelope& b, double abyss_depth, double other_scale);
MatchResult envelope_wasserstein_match(FixedEnvelope& a, FixedEnvelope& b, double flow_distance,
                                       double other_scale);

// ---- Runtime ISA selection ---------------------------------------------------

rust::Str simd_level();
int32_t isa_current_level();
int32_t isa_max_supported_level();
bool isa_level_available(int32_t level);
bool isa_set_level(int32_t level);
rust::Str isa_level_name(int32_t level);

}  // namespace isospec_rs
