#include "isospec/rust/src/ffi.rs.h"

#include <string_view>
#include <utility>

#include "element_lookup.h"
#include "element_tables.h"
#include "fasta.h"
#include "fasta_mods.h"

namespace isospec_rs {

namespace {

std::string_view view(rust::Str s) { return std::string_view(s.data(), s.size()); }

// A rust::Slice may not be built from a null pointer, which an empty envelope
// or composition can hand out: an empty slice is returned instead.
template <typename T>
rust::Slice<const T> slice(const T* data, size_t len)
{
    if (data == nullptr || len == 0)
        return rust::Slice<const T>();
    return rust::Slice<const T>(data, len);
}

IsoSpec::ElementComposition composition(rust::Slice<const int32_t> element_first_index,
                                        rust::Slice<const int32_t> counts)
{
    IsoSpec::ElementComposition c;
    c.element_first_index.assign(element_first_index.begin(), element_first_index.end());
    c.count.assign(counts.begin(), counts.end());
    return c;
}

void copy_out(const IsoSpec::ElementComposition& c, rust::Vec<int32_t>& element_first_index,
              rust::Vec<int32_t>& counts)
{
    element_first_index.clear();
    counts.clear();
    element_first_index.reserve(c.element_first_index.size());
    counts.reserve(c.count.size());
    for (int idx : c.element_first_index)
        element_first_index.push_back(idx);
    for (int cnt : c.count)
        counts.push_back(cnt);
}

template <typename Gen>
bool next_peak(Gen& g, PeakData& out)
{
    if (!g.advanceToNextConfiguration())
        return false;
    out.mass = g.mass();
    out.prob = g.prob();
    out.lprob = g.lprob();
    return true;
}

// Up to masses.size() peaks (all three slices are that long); fewer means
// the generator is exhausted. The output pointers are restrict: otherwise
// every store through them may alias the generator's own state, which is then
// reloaded from memory after each one.
template <typename Gen>
size_t fill_peaks_impl(Gen& g, double* __restrict masses, double* __restrict probs, double* __restrict lprobs,
                       size_t n)
{
    size_t ii = 0;
    for (; ii < n && g.advanceToNextConfiguration(); ii++)
    {
        masses[ii] = g.mass();
        probs[ii] = g.prob();
        lprobs[ii] = g.lprob();
    }
    return ii;
}

template <typename Gen>
size_t fill_peaks(Gen& g, rust::Slice<double> masses, rust::Slice<double> probs, rust::Slice<double> lprobs)
{
    return fill_peaks_impl(g, masses.data(), probs.data(), lprobs.data(), masses.size());
}

std::unique_ptr<std::mt19937> make_rng(bool seeded, uint64_t seed)
{
    if (!seeded)
        return nullptr;
    return std::make_unique<std::mt19937>(IsoSpec::seeded_rng(seed));
}

std::unique_ptr<FixedEnvelope> boxed(FixedEnvelope&& env)
{
    return std::make_unique<FixedEnvelope>(std::move(env));
}

}  // namespace

// ---- Iso ------------------------------------------------------------------

std::unique_ptr<Iso> iso_from_formula(rust::Str formula, bool use_nominal_masses)
{
    return std::make_unique<Iso>(view(formula), use_nominal_masses);
}

std::unique_ptr<Iso> iso_from_fasta(rust::Str sequence, bool use_nominal_masses, bool add_water)
{
    return std::make_unique<Iso>(Iso::FromFASTA(view(sequence), use_nominal_masses, add_water));
}

std::unique_ptr<Iso> iso_from_peptide(rust::Str sequence, const UnimodTable& mods, bool use_nominal_masses,
                                      bool add_water)
{
    return std::make_unique<Iso>(Iso::FromFASTAWithMods(view(sequence), mods, use_nominal_masses, add_water));
}

std::unique_ptr<Iso> iso_from_composition(rust::Slice<const int32_t> element_first_index,
                                          rust::Slice<const int32_t> counts, bool use_nominal_masses)
{
    return std::make_unique<Iso>(
        IsoSpec::build_iso_from_composition(composition(element_first_index, counts), use_nominal_masses));
}

std::unique_ptr<Iso> iso_from_isotopes(rust::Slice<const int32_t> isotope_numbers,
                                       rust::Slice<const int32_t> atom_counts, rust::Slice<const double> masses,
                                       rust::Slice<const double> probs)
{
    // Lengths checked on the Rust side: isotope_numbers and atom_counts alike,
    // masses and probs both the sum of isotope_numbers.
    return std::make_unique<Iso>(static_cast<int>(isotope_numbers.size()), isotope_numbers.data(),
                                 atom_counts.data(), masses.data(), probs.data());
}

std::unique_ptr<Iso> iso_clone(const Iso& iso) { return std::make_unique<Iso>(iso, true); }

void iso_add_element(Iso& iso, int32_t atom_count, rust::Slice<const double> masses,
                     rust::Slice<const double> probs)
{
    iso.addElement(atom_count, static_cast<int>(masses.size()), masses.data(), probs.data());
}

rust::Slice<const int32_t> iso_isotope_numbers(const Iso& iso)
{
    return slice(iso.getIsotopeNumbers(), static_cast<size_t>(iso.getDimNumber()));
}

rust::Slice<const int32_t> iso_atom_counts(const Iso& iso)
{
    return slice(iso.getAtomCounts(), static_cast<size_t>(iso.getDimNumber()));
}

rust::Slice<const double> iso_element_masses(const Iso& iso, int32_t element)
{
    return slice(iso.getElementIsotopeMasses(element), static_cast<size_t>(iso.getIsotopeNumbers()[element]));
}

rust::Slice<const double> iso_element_lprobs(const Iso& iso, int32_t element)
{
    return slice(iso.getElementIsotopeLProbs(element), static_cast<size_t>(iso.getIsotopeNumbers()[element]));
}

void iso_lightest_conf(const Iso& iso, rust::Slice<int32_t> out) { iso.getLightestPeakSignature(out.data()); }
void iso_heaviest_conf(const Iso& iso, rust::Slice<int32_t> out) { iso.getHeaviestPeakSignature(out.data()); }
void iso_monoisotopic_conf(const Iso& iso, rust::Slice<int32_t> out)
{
    iso.getMonoisotopicPeakSignature(out.data());
}

void iso_marginal_log_size_estimates(const Iso& iso, rust::Slice<double> out, double target_total_prob)
{
    iso.saveMarginalLogSizeEstimates(out.data(), target_total_prob);
}

// ---- Compositions and the element table ------------------------------------

void composition_from_formula(rust::Str formula, rust::Vec<int32_t>& element_first_index,
                              rust::Vec<int32_t>& counts)
{
    copy_out(IsoSpec::parse_formula_composition(view(formula)), element_first_index, counts);
}

void composition_from_peptide(rust::Str sequence, const UnimodTable& mods, bool add_water,
                              rust::Vec<int32_t>& element_first_index, rust::Vec<int32_t>& counts)
{
    IsoSpec::ElementComposition c;
    if (add_water)
        IsoSpec::parse_fasta_with_mods_full_into(view(sequence), c, mods);
    else
        IsoSpec::parse_fasta_with_mods_into(view(sequence), c, mods);
    copy_out(c, element_first_index, counts);
}

void composition_from_fasta(rust::Str sequence, bool add_water, rust::Slice<int32_t> chnosse)
{
    // Length 6 checked on the Rust side.
    if (add_water)
        IsoSpec::parse_fasta_full(view(sequence), chnosse.data());
    else
        IsoSpec::parse_fasta(view(sequence), chnosse.data());
}

void composition_expand(rust::Slice<const int32_t> element_first_index, rust::Slice<const int32_t> counts,
                        bool use_nominal_masses, rust::Vec<int32_t>& isotope_numbers,
                        rust::Vec<int32_t>& atom_counts, rust::Vec<double>& masses, rust::Vec<double>& probs)
{
    IsoSpec::ExpandedComposition out;
    IsoSpec::expand_composition_into(composition(element_first_index, counts), out, use_nominal_masses);
    isotope_numbers.clear();
    atom_counts.clear();
    masses.clear();
    probs.clear();
    for (int v : out.isotopeNumbers)
        isotope_numbers.push_back(v);
    for (int v : out.atomCounts)
        atom_counts.push_back(v);
    for (double v : out.isotope_masses)
        masses.push_back(v);
    for (double v : out.isotope_probabilities)
        probs.push_back(v);
}

int32_t element_find(rust::Str symbol)
{
    return IsoSpec::find_element_table_first_index(symbol.data(), symbol.size());
}

int32_t element_isotope_count(int32_t first_index)
{
    if (first_index < 0 || static_cast<size_t>(first_index) >= IsoSpec::isospec_number_of_isotopic_entries)
        return 0;
    return IsoSpec::element_isotope_count(first_index);
}

size_t element_table_len() { return IsoSpec::isospec_number_of_isotopic_entries; }
rust::Str element_table_symbol(size_t row) { return rust::Str(IsoSpec::elem_table_symbol[row]); }
rust::Str element_table_name(size_t row) { return rust::Str(IsoSpec::elem_table_element[row]); }
int32_t element_table_atomic_number(size_t row) { return IsoSpec::elem_table_atomicNo[row]; }
double element_table_mass(size_t row) { return IsoSpec::elem_table_mass[row]; }
double element_table_probability(size_t row) { return IsoSpec::elem_table_probability[row]; }
double element_table_mass_number(size_t row) { return IsoSpec::elem_table_massNo[row]; }

// ---- Unimod ----------------------------------------------------------------

const UnimodTable& unimod_embedded() { return IsoSpec::embedded_unimod_table(); }

const UnimodTable& unimod_for_path(rust::Str path)
{
    return IsoSpec::unimod_table_for_path(std::string(path));
}

std::unique_ptr<UnimodTable> unimod_from_csv(rust::Str csv)
{
    return std::make_unique<UnimodTable>(IsoSpec::parse_unimod_csv(std::string(csv)));
}

size_t unimod_size(const UnimodTable& table) { return table.size(); }

bool unimod_lookup(const UnimodTable& table, uint32_t id, rust::String& name, double& mono_mass,
                   rust::Vec<int32_t>& element_first_index, rust::Vec<int32_t>& counts)
{
    const IsoSpec::UnimodEntry* entry = table.lookup(id);
    if (entry == nullptr)
        return false;
    name = rust::String::lossy(entry->name);
    mono_mass = entry->mono_mass;
    element_first_index.clear();
    counts.clear();
    for (int v : entry->element_first_index)
        element_first_index.push_back(v);
    for (int v : entry->element_delta_count)
        counts.push_back(v);
    return true;
}

// ---- Generators ------------------------------------------------------------

StochasticGenerator::StochasticGenerator(Iso&& iso, size_t no_molecules, double precision, double beta_bias,
                                         std::unique_ptr<std::mt19937> seeded)
    : rng(std::move(seeded)),
      // Unseeded: the library's own default, the thread_local generator of
      // whichever thread is sampling.
      g(std::move(iso), no_molecules, precision, beta_bias, rng ? *rng : IsoSpec::random_gen)
{
}

std::unique_ptr<ThresholdGenerator> threshold_generator_new(std::unique_ptr<Iso> iso, double threshold,
                                                            bool absolute, bool reorder_marginals)
{
    return std::make_unique<ThresholdGenerator>(std::move(*iso), threshold, absolute, reorder_marginals);
}
bool threshold_generator_next(ThresholdGenerator& gen, PeakData& out) { return next_peak(gen.g, out); }
size_t threshold_generator_fill(ThresholdGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs)
{
    return fill_peaks(gen.g, masses, probs, lprobs);
}
void threshold_generator_conf(const ThresholdGenerator& gen, rust::Slice<int32_t> out)
{
    gen.g.get_conf_signature(out.data());
}
size_t threshold_generator_count_confs(ThresholdGenerator& gen) { return gen.g.count_confs(); }
void threshold_generator_reset(ThresholdGenerator& gen) { gen.g.reset(); }

std::unique_ptr<LayeredGenerator> layered_generator_new(std::unique_ptr<Iso> iso, bool reorder_marginals,
                                                        double t_prob_hint)
{
    return std::make_unique<LayeredGenerator>(std::move(*iso), reorder_marginals, t_prob_hint);
}
bool layered_generator_next(LayeredGenerator& gen, PeakData& out) { return next_peak(gen.g, out); }
size_t layered_generator_fill(LayeredGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs)
{
    return fill_peaks(gen.g, masses, probs, lprobs);
}
void layered_generator_conf(const LayeredGenerator& gen, rust::Slice<int32_t> out)
{
    gen.g.get_conf_signature(out.data());
}

std::unique_ptr<OrderedGenerator> ordered_generator_new(std::unique_ptr<Iso> iso)
{
    return std::make_unique<OrderedGenerator>(std::move(*iso));
}
bool ordered_generator_next(OrderedGenerator& gen, PeakData& out) { return next_peak(gen.g, out); }
size_t ordered_generator_fill(OrderedGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs)
{
    return fill_peaks(gen.g, masses, probs, lprobs);
}
void ordered_generator_conf(const OrderedGenerator& gen, rust::Slice<int32_t> out)
{
    gen.g.get_conf_signature(out.data());
}

std::unique_ptr<StochasticGenerator> stochastic_generator_new(std::unique_ptr<Iso> iso, size_t no_molecules,
                                                              double precision, double beta_bias, bool seeded,
                                                              uint64_t seed)
{
    return std::make_unique<StochasticGenerator>(std::move(*iso), no_molecules, precision, beta_bias,
                                                 make_rng(seeded, seed));
}
bool stochastic_generator_next(StochasticGenerator& gen, PeakData& out) { return next_peak(gen.g, out); }
size_t stochastic_generator_fill(StochasticGenerator& gen, rust::Slice<double> masses, rust::Slice<double> probs,
                               rust::Slice<double> lprobs)
{
    return fill_peaks(gen.g, masses, probs, lprobs);
}
void stochastic_generator_conf(const StochasticGenerator& gen, rust::Slice<int32_t> out)
{
    gen.g.get_conf_signature(out.data());
}

// ---- FixedEnvelope -----------------------------------------------------------

std::unique_ptr<FixedEnvelope> envelope_threshold(const Iso& iso, double threshold, bool absolute, bool confs)
{
    return boxed(FixedEnvelope::FromThreshold(iso, threshold, absolute, confs));
}

std::unique_ptr<FixedEnvelope> envelope_total_prob(const Iso& iso, double target, bool optimize, bool confs)
{
    return boxed(FixedEnvelope::FromTotalProb(iso, target, optimize, confs));
}

std::unique_ptr<FixedEnvelope> envelope_stochastic(const Iso& iso, size_t no_molecules, double precision,
                                                   double beta_bias, bool confs, bool seeded, uint64_t seed)
{
    if (seeded)
    {
        std::mt19937 rng = IsoSpec::seeded_rng(seed);
        return boxed(FixedEnvelope::FromStochastic(iso, no_molecules, precision, beta_bias, confs, rng));
    }
    return boxed(FixedEnvelope::FromStochastic(iso, no_molecules, precision, beta_bias, confs));
}

std::unique_ptr<FixedEnvelope> envelope_binned(const Iso& iso, double target, double bin_width,
                                               double bin_middle)
{
    return boxed(FixedEnvelope::Binned(iso, target, bin_width, bin_middle));
}

std::unique_ptr<FixedEnvelope> envelope_from_arrays(rust::Slice<const double> masses,
                                                    rust::Slice<const double> probs)
{
    // Equal lengths checked on the Rust side.
    return boxed(FixedEnvelope::FromArrays(masses.data(), probs.data(), masses.size()));
}

std::unique_ptr<FixedEnvelope> envelope_empty() { return std::make_unique<FixedEnvelope>(); }

std::unique_ptr<FixedEnvelope> envelope_clone(const FixedEnvelope& env)
{
    return std::make_unique<FixedEnvelope>(env);
}

std::unique_ptr<FixedEnvelope> envelope_add(const FixedEnvelope& a, const FixedEnvelope& b) { return boxed(a + b); }

std::unique_ptr<FixedEnvelope> envelope_convolve(const FixedEnvelope& a, const FixedEnvelope& b)
{
    return boxed(a * b);
}

std::unique_ptr<FixedEnvelope> envelope_bin(FixedEnvelope& env, double bin_width, double middle)
{
    return boxed(env.bin(bin_width, middle));
}

std::unique_ptr<EnvelopeList> envelope_list_new() { return std::make_unique<EnvelopeList>(); }

void envelope_list_push(EnvelopeList& list, const FixedEnvelope& env, double intensity)
{
    list.spectra.push_back(&env);
    list.intensities.push_back(intensity);
}

std::unique_ptr<FixedEnvelope> envelope_linear_combination(const EnvelopeList& list)
{
    return boxed(FixedEnvelope::LinearCombination(list.spectra, list.intensities));
}

size_t envelope_len(const FixedEnvelope& env) { return env.confs_no(); }
int32_t envelope_all_dim(const FixedEnvelope& env) { return env.getAllDim(); }
rust::Slice<const double> envelope_masses(const FixedEnvelope& env) { return slice(env.masses(), env.confs_no()); }
rust::Slice<const double> envelope_probs(const FixedEnvelope& env) { return slice(env.probs(), env.confs_no()); }
rust::Slice<const int32_t> envelope_confs(const FixedEnvelope& env)
{
    return slice(env.confs(), env.confs_no() * static_cast<size_t>(env.getAllDim()));
}

void envelope_sort_by_mass(FixedEnvelope& env) { env.sort_by_mass(); }
void envelope_sort_by_prob(FixedEnvelope& env) { env.sort_by_prob(); }
double envelope_total_prob_of(FixedEnvelope& env) { return env.get_total_prob(); }
void envelope_scale(FixedEnvelope& env, double factor) { env.scale(factor); }
void envelope_normalize(FixedEnvelope& env) { env.normalize(); }
void envelope_shift_mass(FixedEnvelope& env, double shift) { env.shift_mass(shift); }

void envelope_resample(FixedEnvelope& env, size_t ionic_current, double beta_bias, bool seeded, uint64_t seed)
{
    if (seeded)
    {
        std::mt19937 rng = IsoSpec::seeded_rng(seed);
        env.resample(ionic_current, beta_bias, rng);
    }
    else
        env.resample(ionic_current, beta_bias);
}

double envelope_average_mass(FixedEnvelope& env) { return env.empiric_average_mass(); }
double envelope_variance(FixedEnvelope& env) { return env.empiric_variance(); }
double envelope_wasserstein(FixedEnvelope& a, FixedEnvelope& b) { return a.WassersteinDistance(b); }
double envelope_oriented_wasserstein(FixedEnvelope& a, FixedEnvelope& b)
{
    return a.OrientedWassersteinDistance(b);
}
double envelope_abyssal_wasserstein(FixedEnvelope& a, FixedEnvelope& b, double abyss_depth, double other_scale)
{
    return a.AbyssalWassersteinDistance(b, abyss_depth, other_scale);
}

MatchResult envelope_wasserstein_match(FixedEnvelope& a, FixedEnvelope& b, double flow_distance,
                                       double other_scale)
{
    const auto [unmatched_this, unmatched_other, flow] = a.WassersteinMatch(b, flow_distance, other_scale);
    return MatchResult{unmatched_this, unmatched_other, flow};
}

// ---- Runtime ISA selection ---------------------------------------------------

rust::Str simd_level() { return rust::Str(IsoSpec::active_simd_level()); }
int32_t isa_current_level() { return IsoSpec::current_isa_level(); }
int32_t isa_max_supported_level() { return IsoSpec::max_supported_isa_level(); }

bool isa_level_available(int32_t level)
{
    if (level < 0 || level >= IsoSpec::ISA_LEVEL_COUNT)
        return false;
    return IsoSpec::isa_level_available(static_cast<IsoSpec::IsaLevel>(level));
}

bool isa_set_level(int32_t level)
{
    if (level < 0 || level >= IsoSpec::ISA_LEVEL_COUNT)
        return false;
    return IsoSpec::set_isa_level(static_cast<IsoSpec::IsaLevel>(level));
}

rust::Str isa_level_name(int32_t level)
{
    return rust::Str(IsoSpec::isa_level_name(static_cast<IsoSpec::IsaLevel>(level)));
}

}  // namespace isospec_rs
