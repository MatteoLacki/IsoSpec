//! Generators: isotopologues one at a time, without storing them.

use std::sync::Arc;

use cxx::UniquePtr;

use crate::error::Result;
use crate::ffi::bridge as sys;
use crate::{Conf, ConfRef, Iso};

/// One isotopologue as a generator reports it.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct Peak {
    /// Mass, in daltons.
    pub mass: f64,
    /// Probability. For a [`StochasticGenerator`], the number of molecules
    /// sampled at this isotopologue instead.
    pub prob: f64,
    /// Natural logarithm of `prob`.
    pub lprob: f64,
}

/// Whether a probability threshold is absolute or relative to the most
/// probable isotopologue.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Threshold {
    /// Isotopologues with probability at least the threshold.
    Absolute,
    /// Isotopologues with probability at least the threshold times that of
    /// the most probable one.
    Relative,
}

mod sealed {
    pub trait Batch {
        fn fill(
            &mut self,
            masses: &mut [f64],
            probs: &mut [f64],
            lprobs: &mut [f64],
        ) -> crate::Result<usize>;
    }
}

const BATCH: usize = 256;

/// The isotopologues of a generator, as an iterator: see the generators'
/// `IntoIterator`. Peaks come from the library 256 at a time.
pub struct Peaks<G: sealed::Batch> {
    generator: G,
    // masses, probs and lprobs, BATCH each.
    buf: Box<[f64]>,
    pos: usize,
    len: usize,
    exhausted: bool,
}

impl<G: sealed::Batch> Peaks<G> {
    fn new(generator: G) -> Self {
        Peaks {
            generator,
            buf: vec![0.0; 3 * BATCH].into_boxed_slice(),
            pos: 0,
            len: 0,
            exhausted: false,
        }
    }

    /// The generator back, positioned after the last peak fetched (which may
    /// be up to a batch past the last one this iterator returned).
    pub fn into_inner(self) -> G {
        self.generator
    }
}

impl<G: sealed::Batch> Peaks<G> {
    // Off the per-peak path, so that `next` stays small enough to inline.
    #[cold]
    #[inline(never)]
    fn refill(&mut self) -> bool {
        if self.exhausted {
            return false;
        }
        let (masses, rest) = self.buf.split_at_mut(BATCH);
        let (probs, lprobs) = rest.split_at_mut(BATCH);
        self.len = self
            .generator
            .fill(masses, probs, lprobs)
            .expect("isotopologue generator");
        self.pos = 0;
        self.exhausted = self.len < BATCH;
        self.len > 0
    }
}

impl<G: sealed::Batch> Iterator for Peaks<G> {
    type Item = Peak;

    #[inline]
    fn next(&mut self) -> Option<Peak> {
        if self.pos == self.len && !self.refill() {
            return None;
        }
        let ii = self.pos;
        self.pos += 1;
        Some(Peak {
            mass: self.buf[ii],
            prob: self.buf[BATCH + ii],
            lprob: self.buf[2 * BATCH + ii],
        })
    }
}

impl<G: sealed::Batch> std::iter::FusedIterator for Peaks<G> {}

impl<G: sealed::Batch + std::fmt::Debug> std::fmt::Debug for Peaks<G> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Peaks")
            .field("generator", &self.generator)
            .finish_non_exhaustive()
    }
}

macro_rules! generator_common {
    ($name:ident, $raw:ident, $next:ident, $fill:ident, $conf:ident) => {
        impl $name {
            fn wrap(inner: UniquePtr<sys::$raw>, isotope_numbers: Arc<[i32]>) -> Self {
                let all_dim = isotope_numbers.iter().map(|&n| n as usize).sum();
                $name {
                    inner,
                    isotope_numbers,
                    all_dim,
                    current: None,
                    done: false,
                    conf_buf: Vec::new(),
                }
            }

            /// Moves to the next isotopologue and returns it, or `None` once
            /// there are no more. An error is a C++ allocation failure.
            pub fn advance(&mut self) -> Result<Option<Peak>> {
                let mut out = sys::PeakData {
                    mass: 0.0,
                    prob: 0.0,
                    lprob: 0.0,
                };
                self.current = None;
                // Never asked again once it has said there is no more.
                if self.done {
                    return Ok(None);
                }
                if !sys::$next(self.inner.pin_mut(), &mut out)? {
                    self.done = true;
                    return Ok(None);
                }
                // Returned from the local, not read back from self.current:
                // a store and an immediate reload of a different width cannot
                // be forwarded, and that stall was most of the cost of a step.
                let peak = Peak {
                    mass: out.mass,
                    prob: out.prob,
                    lprob: out.lprob,
                };
                self.current = Some(peak);
                Ok(Some(peak))
            }

            /// The isotopologue the last [`advance`](Self::advance) (or
            /// iterator step) moved to, if any.
            pub fn current(&self) -> Option<Peak> {
                self.current
            }

            /// The configuration of the current isotopologue: `None` before
            /// the first step and after the last.
            pub fn conf(&mut self) -> Option<ConfRef<'_>> {
                self.current?;
                self.conf_buf.resize(self.all_dim, 0);
                sys::$conf(self.inner.as_ref().expect("non-null"), &mut self.conf_buf);
                Some(ConfRef::new(&self.conf_buf, &self.isotope_numbers))
            }

            /// As [`conf`](Self::conf), owned.
            pub fn conf_owned(&mut self) -> Option<Conf> {
                let counts = self.conf()?.counts().to_vec();
                Some(Conf::from_parts(counts, self.isotope_numbers.clone()))
            }

            /// Writes the current configuration's counts, flat, into `out`,
            /// which must be [`all_dim`](Self::all_dim) long. Returns `false`
            /// (and writes nothing) when there is no current isotopologue.
            /// Panics if `out` has the wrong length.
            pub fn write_conf(&self, out: &mut [i32]) -> bool {
                assert_eq!(
                    out.len(),
                    self.all_dim,
                    "configuration buffer of the wrong length"
                );
                if self.current.is_none() {
                    return false;
                }
                sys::$conf(self.inner.as_ref().expect("non-null"), out);
                true
            }

            /// Number of isotopes of each element of the molecule.
            pub fn isotope_numbers(&self) -> &[i32] {
                &self.isotope_numbers
            }

            /// Length of a configuration's flat counts.
            pub fn all_dim(&self) -> usize {
                self.all_dim
            }
        }

        impl sealed::Batch for $name {
            fn fill(
                &mut self,
                masses: &mut [f64],
                probs: &mut [f64],
                lprobs: &mut [f64],
            ) -> Result<usize> {
                self.current = None;
                if self.done {
                    return Ok(0);
                }
                let n = sys::$fill(self.inner.pin_mut(), masses, probs, lprobs)?;
                if n < masses.len() {
                    self.done = true;
                }
                Ok(n)
            }
        }

        /// Iterates over the remaining isotopologues, fetching them from the
        /// library in batches (far faster than one
        /// [`advance`]($name::advance) per peak). Panics if the C++ library
        /// fails to allocate memory.
        impl IntoIterator for $name {
            type Item = Peak;
            type IntoIter = Peaks<$name>;
            fn into_iter(self) -> Peaks<$name> {
                Peaks::new(self)
            }
        }

        impl std::fmt::Debug for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("isotope_numbers", &self.isotope_numbers)
                    .field("current", &self.current)
                    .finish_non_exhaustive()
            }
        }
    };
}

/// Every isotopologue with probability above a threshold, in no particular
/// order. The fastest generator.
///
/// ```
/// use isospec::{Iso, Threshold, ThresholdGenerator};
/// let iso = Iso::from_formula("C6H12O6")?;
/// let total: f64 = ThresholdGenerator::new(iso, 1e-6, Threshold::Absolute)?.into_iter().map(|p| p.prob).sum();
/// assert!(total > 0.9999 && total <= 1.0);
/// # Ok::<(), isospec::Error>(())
/// ```
pub struct ThresholdGenerator {
    inner: UniquePtr<sys::ThresholdGenerator>,
    isotope_numbers: Arc<[i32]>,
    all_dim: usize,
    current: Option<Peak>,
    done: bool,
    conf_buf: Vec<i32>,
}

/// Options of a [`ThresholdGenerator`]: start from
/// [`default`](Default::default) and change what you need, so that options
/// added later leave your code compiling.
///
/// ```
/// use isospec::ThresholdOptions;
/// let options = ThresholdOptions::default().with_reorder_marginals(false);
/// assert!(!options.reorder_marginals);
/// ```
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct ThresholdOptions {
    /// Reorder the elements internally for speed; the set of isotopologues is
    /// the same either way, their order is not. Default: `true`.
    pub reorder_marginals: bool,
}

impl Default for ThresholdOptions {
    fn default() -> Self {
        ThresholdOptions {
            reorder_marginals: true,
        }
    }
}

impl ThresholdOptions {
    /// Sets [`reorder_marginals`](Self::reorder_marginals).
    pub fn with_reorder_marginals(mut self, reorder: bool) -> Self {
        self.reorder_marginals = reorder;
        self
    }
}

generator_common!(
    ThresholdGenerator,
    ThresholdGenerator,
    threshold_generator_next,
    threshold_generator_fill,
    threshold_generator_conf
);

impl ThresholdGenerator {
    /// The isotopologues of `iso` at or above `threshold`.
    pub fn new(iso: Iso, threshold: f64, mode: Threshold) -> Result<Self> {
        ThresholdGenerator::with_options(iso, threshold, mode, ThresholdOptions::default())
    }

    /// As [`new`](Self::new), with options.
    pub fn with_options(
        iso: Iso,
        threshold: f64,
        mode: Threshold,
        options: ThresholdOptions,
    ) -> Result<Self> {
        let isotope_numbers = iso.isotope_numbers_arc();
        let inner = sys::threshold_generator_new(
            iso.into_raw(),
            threshold,
            mode == Threshold::Absolute,
            options.reorder_marginals,
        )?;
        Ok(ThresholdGenerator::wrap(inner, isotope_numbers))
    }

    /// Counts the isotopologues above the threshold, running through all of
    /// them, and [`reset`](Self::reset)s the generator.
    pub fn count_confs(&mut self) -> Result<usize> {
        self.current = None;
        self.done = false;
        let count = sys::threshold_generator_count_confs(self.inner.pin_mut())?;
        // The library resets after counting only when it had to walk the
        // configurations; a one-element molecule's count is a lookup.
        sys::threshold_generator_reset(self.inner.pin_mut())?;
        Ok(count)
    }

    /// Starts again from the beginning.
    pub fn reset(&mut self) -> Result<()> {
        self.current = None;
        self.done = false;
        sys::threshold_generator_reset(self.inner.pin_mut())?;
        Ok(())
    }
}

/// Isotopologues in layers of decreasing probability: each layer's are less
/// probable than the previous layer's, in no particular order within a layer.
/// Runs through the whole distribution if not stopped; stop once enough
/// probability is covered.
///
/// ```
/// use isospec::{Iso, LayeredGenerator};
/// let mut covered = 0.0;
/// let mut count = 0;
/// for peak in LayeredGenerator::new(Iso::from_formula("C100H202")?)? {
///     covered += peak.prob;
///     count += 1;
///     if covered >= 0.99 { break; }
/// }
/// assert!(count > 1);
/// # Ok::<(), isospec::Error>(())
/// ```
pub struct LayeredGenerator {
    inner: UniquePtr<sys::LayeredGenerator>,
    isotope_numbers: Arc<[i32]>,
    all_dim: usize,
    current: Option<Peak>,
    done: bool,
    conf_buf: Vec<i32>,
}

/// Options of a [`LayeredGenerator`]: start from
/// [`default`](Default::default) and change what you need.
///
/// ```
/// use isospec::LayeredOptions;
/// let options = LayeredOptions::default().with_prob_hint(0.999);
/// assert_eq!(options.prob_hint, 0.999);
/// ```
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct LayeredOptions {
    /// Reorder the elements internally for speed. Default: `true`.
    pub reorder_marginals: bool,
    /// The probability the caller expects to cover, to size the first layers.
    /// Default: `0.99`.
    pub prob_hint: f64,
}

impl Default for LayeredOptions {
    fn default() -> Self {
        LayeredOptions {
            reorder_marginals: true,
            prob_hint: 0.99,
        }
    }
}

impl LayeredOptions {
    /// Sets [`reorder_marginals`](Self::reorder_marginals).
    pub fn with_reorder_marginals(mut self, reorder: bool) -> Self {
        self.reorder_marginals = reorder;
        self
    }

    /// Sets [`prob_hint`](Self::prob_hint).
    pub fn with_prob_hint(mut self, prob_hint: f64) -> Self {
        self.prob_hint = prob_hint;
        self
    }
}

generator_common!(
    LayeredGenerator,
    LayeredGenerator,
    layered_generator_next,
    layered_generator_fill,
    layered_generator_conf
);

impl LayeredGenerator {
    /// The isotopologues of `iso`, layer by layer.
    pub fn new(iso: Iso) -> Result<Self> {
        LayeredGenerator::with_options(iso, LayeredOptions::default())
    }

    /// As [`new`](Self::new), with options.
    pub fn with_options(iso: Iso, options: LayeredOptions) -> Result<Self> {
        let isotope_numbers = iso.isotope_numbers_arc();
        let inner = sys::layered_generator_new(
            iso.into_raw(),
            options.reorder_marginals,
            options.prob_hint,
        )?;
        Ok(LayeredGenerator::wrap(inner, isotope_numbers))
    }
}

/// Isotopologues in strictly decreasing order of probability. Slower than
/// the others (a priority queue); use it when the order itself matters.
///
/// ```
/// use isospec::{Iso, OrderedGenerator};
/// let probs: Vec<f64> = OrderedGenerator::new(Iso::from_formula("H2O1")?)?.into_iter().take(3).map(|p| p.prob).collect();
/// assert!(probs[0] >= probs[1] && probs[1] >= probs[2]);
/// # Ok::<(), isospec::Error>(())
/// ```
pub struct OrderedGenerator {
    inner: UniquePtr<sys::OrderedGenerator>,
    isotope_numbers: Arc<[i32]>,
    all_dim: usize,
    current: Option<Peak>,
    done: bool,
    conf_buf: Vec<i32>,
}

generator_common!(
    OrderedGenerator,
    OrderedGenerator,
    ordered_generator_next,
    ordered_generator_fill,
    ordered_generator_conf
);

impl OrderedGenerator {
    /// The isotopologues of `iso`, most probable first.
    pub fn new(iso: Iso) -> Result<Self> {
        let isotope_numbers = iso.isotope_numbers_arc();
        let inner = sys::ordered_generator_new(iso.into_raw())?;
        Ok(OrderedGenerator::wrap(inner, isotope_numbers))
    }
}

/// A simulated measurement: `molecules` molecules sampled from the
/// distribution, reported as the isotopologues that received any, each with
/// its count in [`Peak::prob`].
///
/// ```
/// use isospec::{Iso, StochasticGenerator, StochasticOptions};
/// let options = StochasticOptions::default().with_seed(42);
/// let total: f64 = StochasticGenerator::with_options(Iso::from_formula("C100")?, 1000, options)?
///     .into_iter()
///     .map(|p| p.prob)
///     .sum();
/// assert!(total <= 1000.0);
/// # Ok::<(), isospec::Error>(())
/// ```
pub struct StochasticGenerator {
    inner: UniquePtr<sys::StochasticGenerator>,
    isotope_numbers: Arc<[i32]>,
    all_dim: usize,
    current: Option<Peak>,
    done: bool,
    conf_buf: Vec<i32>,
}

/// Options of stochastic sampling: [`StochasticGenerator`] and
/// [`Envelope::stochastic`](crate::EnvelopeBuilder::stochastic). Start from
/// [`default`](Default::default) and change what you need.
///
/// ```
/// use isospec::StochasticOptions;
/// let options = StochasticOptions::default().with_seed(42).with_precision(0.999);
/// assert_eq!(options.seed, Some(42));
/// ```
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct StochasticOptions {
    /// The share of the distribution sampled from, in (0, 1]; molecules that
    /// would fall in the remaining tail are dropped. Default: `0.9999`.
    pub precision: f64,
    /// Tuning of the sampling algorithm (nonnegative); changes speed, not the
    /// distribution of the result. Default: `5.0`.
    pub beta_bias: f64,
    /// Seed of the random number generator, for a reproducible result (the
    /// same as C++ and Python give for the same seed). `None` draws from a
    /// randomly seeded generator. Default: `None`.
    pub seed: Option<u64>,
}

impl Default for StochasticOptions {
    fn default() -> Self {
        StochasticOptions {
            precision: 0.9999,
            beta_bias: 5.0,
            seed: None,
        }
    }
}

impl StochasticOptions {
    /// Sets [`precision`](Self::precision).
    pub fn with_precision(mut self, precision: f64) -> Self {
        self.precision = precision;
        self
    }

    /// Sets [`beta_bias`](Self::beta_bias).
    pub fn with_beta_bias(mut self, beta_bias: f64) -> Self {
        self.beta_bias = beta_bias;
        self
    }

    /// Sets [`seed`](Self::seed), for a reproducible result.
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }
}

generator_common!(
    StochasticGenerator,
    StochasticGenerator,
    stochastic_generator_next,
    stochastic_generator_fill,
    stochastic_generator_conf
);

impl StochasticGenerator {
    /// `molecules` molecules of `iso`, sampled.
    pub fn new(iso: Iso, molecules: usize) -> Result<Self> {
        StochasticGenerator::with_options(iso, molecules, StochasticOptions::default())
    }

    /// As [`new`](Self::new), with options.
    pub fn with_options(iso: Iso, molecules: usize, options: StochasticOptions) -> Result<Self> {
        let isotope_numbers = iso.isotope_numbers_arc();
        let inner = sys::stochastic_generator_new(
            iso.into_raw(),
            molecules,
            options.precision,
            options.beta_bias,
            options.seed.is_some(),
            options.seed.unwrap_or(0),
        )?;
        Ok(StochasticGenerator::wrap(inner, isotope_numbers))
    }

    /// Number of molecules at the current isotopologue.
    pub fn molecule_count(&self) -> Option<u64> {
        self.current.map(|p| p.prob as u64)
    }
}
