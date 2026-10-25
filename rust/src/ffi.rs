//! The cxx bridge: raw bindings to rust/shim/shim.h. Nothing here is public;
//! the safe API is built on it in the other modules, which also make every
//! check the C++ cannot make itself (slice lengths matching counts, element
//! indices in range) before calling in.

#[cxx::bridge(namespace = "isospec_rs")]
pub(crate) mod bridge {
    struct PeakData {
        mass: f64,
        prob: f64,
        lprob: f64,
    }

    struct MatchResult {
        unmatched_this: f64,
        unmatched_other: f64,
        flow: f64,
    }

    unsafe extern "C++" {
        include!("isospec/rust/shim/shim.h");

        #[namespace = "IsoSpec"]
        type Iso;
        #[namespace = "IsoSpec"]
        type FixedEnvelope;
        #[namespace = "IsoSpec"]
        type UnimodTable;
        type ThresholdGenerator;
        type LayeredGenerator;
        type OrderedGenerator;
        type StochasticGenerator;
        type EnvelopeList;

        // ---- Iso ----
        fn iso_from_formula(formula: &str, use_nominal_masses: bool) -> Result<UniquePtr<Iso>>;
        fn iso_from_fasta(
            sequence: &str,
            use_nominal_masses: bool,
            add_water: bool,
        ) -> Result<UniquePtr<Iso>>;
        fn iso_from_peptide(
            sequence: &str,
            mods: &UnimodTable,
            use_nominal_masses: bool,
            add_water: bool,
        ) -> Result<UniquePtr<Iso>>;
        fn iso_from_composition(
            element_first_index: &[i32],
            counts: &[i32],
            use_nominal_masses: bool,
        ) -> Result<UniquePtr<Iso>>;
        fn iso_from_isotopes(
            isotope_numbers: &[i32],
            atom_counts: &[i32],
            masses: &[f64],
            probs: &[f64],
        ) -> Result<UniquePtr<Iso>>;
        fn iso_clone(iso: &Iso) -> Result<UniquePtr<Iso>>;
        fn iso_add_element(
            iso: Pin<&mut Iso>,
            atom_count: i32,
            masses: &[f64],
            probs: &[f64],
        ) -> Result<()>;

        fn getDimNumber(self: &Iso) -> i32;
        fn getAllDim(self: &Iso) -> i32;
        fn getLightestPeakMass(self: &Iso) -> f64;
        fn getLightestPeakLProb(self: &Iso) -> f64;
        fn getHeaviestPeakMass(self: &Iso) -> f64;
        fn getHeaviestPeakLProb(self: &Iso) -> f64;
        fn getMonoisotopicPeakMass(self: &Iso) -> f64;
        fn getMonoisotopicPeakLProb(self: &Iso) -> f64;
        fn getModeLProb(self: &Iso) -> f64;
        fn getModeMass(self: &Iso) -> f64;
        fn getUnlikeliestPeakLProb(self: &Iso) -> f64;
        fn getTheoreticalAverageMass(self: &Iso) -> f64;
        fn variance(self: &Iso) -> f64;
        fn getMassRoundingBound(self: &Iso) -> f64;
        fn getMaxAtomMassSpan(self: &Iso) -> f64;
        fn iso_isotope_numbers(iso: &Iso) -> &[i32];
        fn iso_atom_counts(iso: &Iso) -> &[i32];
        fn iso_element_masses(iso: &Iso, element: i32) -> &[f64];
        fn iso_element_lprobs(iso: &Iso, element: i32) -> &[f64];
        fn iso_lightest_conf(iso: &Iso, out: &mut [i32]);
        fn iso_heaviest_conf(iso: &Iso, out: &mut [i32]);
        fn iso_monoisotopic_conf(iso: &Iso, out: &mut [i32]);
        fn iso_marginal_log_size_estimates(
            iso: &Iso,
            out: &mut [f64],
            target_total_prob: f64,
        ) -> Result<()>;

        // ---- Compositions and the element table ----
        fn composition_from_formula(
            formula: &str,
            element_first_index: &mut Vec<i32>,
            counts: &mut Vec<i32>,
        ) -> Result<()>;
        fn composition_from_peptide(
            sequence: &str,
            mods: &UnimodTable,
            add_water: bool,
            element_first_index: &mut Vec<i32>,
            counts: &mut Vec<i32>,
        ) -> Result<()>;
        fn composition_from_fasta(sequence: &str, add_water: bool, chnosse: &mut [i32]);
        #[allow(clippy::too_many_arguments)]
        fn composition_expand(
            element_first_index: &[i32],
            counts: &[i32],
            use_nominal_masses: bool,
            isotope_numbers: &mut Vec<i32>,
            atom_counts: &mut Vec<i32>,
            masses: &mut Vec<f64>,
            probs: &mut Vec<f64>,
        ) -> Result<()>;
        fn element_find(symbol: &str) -> i32;
        fn element_isotope_count(first_index: i32) -> i32;
        fn element_table_len() -> usize;
        fn element_table_symbol(row: usize) -> &'static str;
        fn element_table_name(row: usize) -> &'static str;
        fn element_table_atomic_number(row: usize) -> i32;
        fn element_table_mass(row: usize) -> f64;
        fn element_table_probability(row: usize) -> f64;
        fn element_table_mass_number(row: usize) -> f64;

        // ---- Unimod ----
        fn unimod_embedded() -> Result<&'static UnimodTable>;
        fn unimod_for_path(path: &str) -> Result<&'static UnimodTable>;
        fn unimod_from_csv(csv: &str) -> Result<UniquePtr<UnimodTable>>;
        fn unimod_size(table: &UnimodTable) -> usize;
        fn unimod_lookup(
            table: &UnimodTable,
            id: u32,
            name: &mut String,
            mono_mass: &mut f64,
            element_first_index: &mut Vec<i32>,
            counts: &mut Vec<i32>,
        ) -> Result<bool>;

        // ---- Generators ----
        fn threshold_generator_new(
            iso: UniquePtr<Iso>,
            threshold: f64,
            absolute: bool,
            reorder_marginals: bool,
        ) -> Result<UniquePtr<ThresholdGenerator>>;
        fn threshold_generator_next(
            g: Pin<&mut ThresholdGenerator>,
            out: &mut PeakData,
        ) -> Result<bool>;
        fn threshold_generator_fill(
            g: Pin<&mut ThresholdGenerator>,
            masses: &mut [f64],
            probs: &mut [f64],
            lprobs: &mut [f64],
        ) -> Result<usize>;
        fn threshold_generator_conf(g: &ThresholdGenerator, out: &mut [i32]);
        fn threshold_generator_count_confs(g: Pin<&mut ThresholdGenerator>) -> Result<usize>;
        fn threshold_generator_reset(g: Pin<&mut ThresholdGenerator>) -> Result<()>;

        fn layered_generator_new(
            iso: UniquePtr<Iso>,
            reorder_marginals: bool,
            t_prob_hint: f64,
        ) -> Result<UniquePtr<LayeredGenerator>>;
        fn layered_generator_next(
            g: Pin<&mut LayeredGenerator>,
            out: &mut PeakData,
        ) -> Result<bool>;
        fn layered_generator_fill(
            g: Pin<&mut LayeredGenerator>,
            masses: &mut [f64],
            probs: &mut [f64],
            lprobs: &mut [f64],
        ) -> Result<usize>;
        fn layered_generator_conf(g: &LayeredGenerator, out: &mut [i32]);

        fn ordered_generator_new(iso: UniquePtr<Iso>) -> Result<UniquePtr<OrderedGenerator>>;
        fn ordered_generator_next(
            g: Pin<&mut OrderedGenerator>,
            out: &mut PeakData,
        ) -> Result<bool>;
        fn ordered_generator_fill(
            g: Pin<&mut OrderedGenerator>,
            masses: &mut [f64],
            probs: &mut [f64],
            lprobs: &mut [f64],
        ) -> Result<usize>;
        fn ordered_generator_conf(g: &OrderedGenerator, out: &mut [i32]);

        fn stochastic_generator_new(
            iso: UniquePtr<Iso>,
            no_molecules: usize,
            precision: f64,
            beta_bias: f64,
            seeded: bool,
            seed: u64,
        ) -> Result<UniquePtr<StochasticGenerator>>;
        fn stochastic_generator_next(
            g: Pin<&mut StochasticGenerator>,
            out: &mut PeakData,
        ) -> Result<bool>;
        fn stochastic_generator_fill(
            g: Pin<&mut StochasticGenerator>,
            masses: &mut [f64],
            probs: &mut [f64],
            lprobs: &mut [f64],
        ) -> Result<usize>;
        fn stochastic_generator_conf(g: &StochasticGenerator, out: &mut [i32]);

        // ---- FixedEnvelope ----
        fn envelope_threshold(
            iso: &Iso,
            threshold: f64,
            absolute: bool,
            confs: bool,
        ) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_total_prob(
            iso: &Iso,
            target: f64,
            optimize: bool,
            confs: bool,
        ) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_stochastic(
            iso: &Iso,
            no_molecules: usize,
            precision: f64,
            beta_bias: f64,
            confs: bool,
            seeded: bool,
            seed: u64,
        ) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_binned(
            iso: &Iso,
            target: f64,
            bin_width: f64,
            bin_middle: f64,
        ) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_from_arrays(masses: &[f64], probs: &[f64]) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_empty() -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_clone(env: &FixedEnvelope) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_add(a: &FixedEnvelope, b: &FixedEnvelope) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_convolve(
            a: &FixedEnvelope,
            b: &FixedEnvelope,
        ) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_bin(
            env: Pin<&mut FixedEnvelope>,
            bin_width: f64,
            middle: f64,
        ) -> Result<UniquePtr<FixedEnvelope>>;
        fn envelope_list_new() -> Result<UniquePtr<EnvelopeList>>;
        fn envelope_list_push(
            list: Pin<&mut EnvelopeList>,
            env: &FixedEnvelope,
            intensity: f64,
        ) -> Result<()>;
        fn envelope_linear_combination(list: &EnvelopeList) -> Result<UniquePtr<FixedEnvelope>>;

        fn envelope_len(env: &FixedEnvelope) -> usize;
        fn envelope_all_dim(env: &FixedEnvelope) -> i32;
        fn envelope_masses(env: &FixedEnvelope) -> &[f64];
        fn envelope_probs(env: &FixedEnvelope) -> &[f64];
        fn envelope_confs(env: &FixedEnvelope) -> &[i32];

        fn envelope_sort_by_mass(env: Pin<&mut FixedEnvelope>) -> Result<()>;
        fn envelope_sort_by_prob(env: Pin<&mut FixedEnvelope>) -> Result<()>;
        fn envelope_total_prob_of(env: Pin<&mut FixedEnvelope>) -> f64;
        fn envelope_scale(env: Pin<&mut FixedEnvelope>, factor: f64) -> Result<()>;
        fn envelope_normalize(env: Pin<&mut FixedEnvelope>) -> Result<()>;
        fn envelope_shift_mass(env: Pin<&mut FixedEnvelope>, shift: f64) -> Result<()>;
        fn envelope_resample(
            env: Pin<&mut FixedEnvelope>,
            ionic_current: usize,
            beta_bias: f64,
            seeded: bool,
            seed: u64,
        ) -> Result<()>;
        fn envelope_average_mass(env: Pin<&mut FixedEnvelope>) -> Result<f64>;
        fn envelope_variance(env: Pin<&mut FixedEnvelope>) -> Result<f64>;
        fn envelope_wasserstein(
            a: Pin<&mut FixedEnvelope>,
            b: Pin<&mut FixedEnvelope>,
        ) -> Result<f64>;
        fn envelope_oriented_wasserstein(
            a: Pin<&mut FixedEnvelope>,
            b: Pin<&mut FixedEnvelope>,
        ) -> Result<f64>;
        fn envelope_abyssal_wasserstein(
            a: Pin<&mut FixedEnvelope>,
            b: Pin<&mut FixedEnvelope>,
            abyss_depth: f64,
            other_scale: f64,
        ) -> Result<f64>;
        fn envelope_wasserstein_match(
            a: Pin<&mut FixedEnvelope>,
            b: Pin<&mut FixedEnvelope>,
            flow_distance: f64,
            other_scale: f64,
        ) -> Result<MatchResult>;

        // ---- Runtime ISA selection ----
        fn simd_level() -> &'static str;
        fn isa_current_level() -> i32;
        fn isa_max_supported_level() -> i32;
        fn isa_level_available(level: i32) -> bool;
        fn isa_set_level(level: i32) -> bool;
        fn isa_level_name(level: i32) -> &'static str;
    }
}

// The C++ objects are safe to move between threads (nothing in them is tied
// to the thread that built them -- see AGENTS.md, Threading), and the const
// members of Iso, FixedEnvelope and UnimodTable may run concurrently. The
// safe wrappers decide what is exposed through & and what through &mut.
unsafe impl Send for bridge::Iso {}
unsafe impl Sync for bridge::Iso {}
unsafe impl Send for bridge::FixedEnvelope {}
unsafe impl Sync for bridge::FixedEnvelope {}
unsafe impl Send for bridge::UnimodTable {}
unsafe impl Sync for bridge::UnimodTable {}
unsafe impl Send for bridge::ThresholdGenerator {}
unsafe impl Send for bridge::LayeredGenerator {}
unsafe impl Send for bridge::OrderedGenerator {}
unsafe impl Send for bridge::StochasticGenerator {}
