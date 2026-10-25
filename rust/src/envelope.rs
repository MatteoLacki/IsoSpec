//! Envelopes: isotopic distributions computed whole, as arrays of peaks.

use std::fmt;
use std::ops::{Add, Mul};
use std::sync::Arc;

use cxx::UniquePtr;

use crate::error::{Error, Result};
use crate::ffi::bridge as sys;
use crate::generator::{StochasticOptions, Threshold};
use crate::{ConfRef, Iso};

/// An isotopic distribution as arrays: a mass and a probability per peak,
/// and optionally each peak's isotopologue configuration.
///
/// Built from an [`Iso`] (see [`Iso::envelope`] for the options), from
/// arrays, or from other envelopes. Peaks are in no particular order until
/// sorted.
///
/// ```
/// use isospec::{Envelope, Iso};
/// let iso = Iso::from_formula("C100H202")?;
/// let mut env = Envelope::from_total_prob(&iso, 0.99)?;
/// assert!(env.total_prob() >= 0.99);
/// env.sort_by_mass();
/// assert!(env.masses().windows(2).all(|w| w[0] <= w[1]));
/// # Ok::<(), isospec::Error>(())
/// ```
pub struct Envelope {
    inner: UniquePtr<sys::FixedEnvelope>,
    // The molecule's isotopes per element, when the peaks carry configurations.
    isotope_numbers: Option<Arc<[i32]>>,
}

/// The result of [`Envelope::wasserstein_match`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Match {
    /// Probability of this envelope left unmatched.
    pub unmatched_this: f64,
    /// Probability of the other envelope (scaled) left unmatched.
    pub unmatched_other: f64,
    /// Probability matched between the two.
    pub flow: f64,
}

/// Options for computing an [`Envelope`] from an [`Iso`]: whether to keep
/// configurations. Start one with [`Iso::envelope`].
#[derive(Clone, Copy, Debug)]
pub struct EnvelopeBuilder<'a> {
    iso: &'a Iso,
    confs: bool,
}

impl Iso {
    /// Starts computing an [`Envelope`] of this molecule.
    ///
    /// ```
    /// use isospec::{Iso, Threshold};
    /// let iso = Iso::from_formula("H2O1")?;
    /// let env = iso.envelope().with_confs().threshold(1e-3, Threshold::Absolute)?;
    /// let conf = env.conf(0).unwrap();
    /// assert_eq!(conf.len(), 2); // hydrogen and oxygen
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn envelope(&self) -> EnvelopeBuilder<'_> {
        EnvelopeBuilder {
            iso: self,
            confs: false,
        }
    }
}

impl<'a> EnvelopeBuilder<'a> {
    /// Keep each peak's configuration (more memory; see [`Envelope::conf`]).
    pub fn with_confs(mut self) -> Self {
        self.confs = true;
        self
    }

    fn finish(self, inner: UniquePtr<sys::FixedEnvelope>) -> Envelope {
        let layout = if self.confs {
            Some(self.iso.isotope_numbers_arc())
        } else {
            None
        };
        Envelope::adopt(inner, layout.as_ref())
    }

    /// Every isotopologue with probability at least `threshold` (or that
    /// times the most probable one's, with [`Threshold::Relative`]).
    pub fn threshold(self, threshold: f64, mode: Threshold) -> Result<Envelope> {
        let inner = sys::envelope_threshold(
            self.iso.raw(),
            threshold,
            mode == Threshold::Absolute,
            self.confs,
        )?;
        Ok(self.finish(inner))
    }

    /// The smallest set of isotopologues whose probabilities sum to at least
    /// `target` (in [0, 1]).
    pub fn total_prob(self, target: f64) -> Result<Envelope> {
        let inner = sys::envelope_total_prob(self.iso.raw(), target, true, self.confs)?;
        Ok(self.finish(inner))
    }

    /// Like [`total_prob`](Self::total_prob), but stopping at the end of the
    /// probability layer that reaches `target` rather than trimming to the
    /// smallest set: faster, and a few more peaks.
    pub fn total_prob_layered(self, target: f64) -> Result<Envelope> {
        let inner = sys::envelope_total_prob(self.iso.raw(), target, false, self.confs)?;
        Ok(self.finish(inner))
    }

    /// A simulated spectrum of `molecules` molecules: each peak's probability
    /// is its molecule count. See [`StochasticOptions`].
    pub fn stochastic(self, molecules: usize, options: StochasticOptions) -> Result<Envelope> {
        let inner = sys::envelope_stochastic(
            self.iso.raw(),
            molecules,
            options.precision,
            options.beta_bias,
            self.confs,
            options.seed.is_some(),
            options.seed.unwrap_or(0),
        )?;
        Ok(self.finish(inner))
    }
}

impl Envelope {
    // Keeps `layout` only if the C++ envelope really carries configurations of
    // that shape: operations that drop them, or envelopes with none, get none.
    fn adopt(inner: UniquePtr<sys::FixedEnvelope>, layout: Option<&Arc<[i32]>>) -> Envelope {
        debug_assert!(!inner.is_null());
        let env = inner.as_ref().expect("non-null");
        let len = sys::envelope_len(env);
        let all_dim = sys::envelope_all_dim(env).max(0) as usize;
        let isotope_numbers = layout.filter(|layout| {
            let layout_dim: usize = layout.iter().map(|&n| n as usize).sum();
            layout_dim == all_dim
                && sys::envelope_confs(env).len() == len * all_dim
                && (len > 0 || all_dim > 0)
        });
        Envelope {
            inner,
            isotope_numbers: isotope_numbers.cloned(),
        }
    }

    fn raw(&self) -> &sys::FixedEnvelope {
        self.inner.as_ref().expect("non-null")
    }

    /// Every isotopologue of `iso` with probability at least `threshold`,
    /// without configurations. See [`EnvelopeBuilder::threshold`].
    pub fn from_threshold(iso: &Iso, threshold: f64, mode: Threshold) -> Result<Envelope> {
        iso.envelope().threshold(threshold, mode)
    }

    /// The smallest set of isotopologues of `iso` covering `target`, without
    /// configurations. See [`EnvelopeBuilder::total_prob`].
    pub fn from_total_prob(iso: &Iso, target: f64) -> Result<Envelope> {
        iso.envelope().total_prob(target)
    }

    /// The distribution of `iso` binned on a grid of width `bin_width` (one
    /// bin centred on `bin_middle`), computed directly rather than by binning
    /// an envelope: isotopologues are enumerated until `target` of the
    /// probability is covered, and each bin reports the mass at its centre.
    pub fn binned(iso: &Iso, target: f64, bin_width: f64, bin_middle: f64) -> Result<Envelope> {
        Ok(Envelope::adopt(
            sys::envelope_binned(iso.raw(), target, bin_width, bin_middle)?,
            None,
        ))
    }

    /// An envelope of the given peaks. Masses must be finite (|mass| at most
    /// 1e100 or so), probabilities finite and nonnegative, and the two slices
    /// of equal length.
    pub fn from_arrays(masses: &[f64], probs: &[f64]) -> Result<Envelope> {
        if masses.len() != probs.len() {
            return Err(Error::invalid(format!(
                "{} masses but {} probabilities",
                masses.len(),
                probs.len()
            )));
        }
        Ok(Envelope::adopt(
            sys::envelope_from_arrays(masses, probs)?,
            None,
        ))
    }

    /// The sum of `intensity × envelope` over the given pairs.
    pub fn linear_combination(terms: &[(&Envelope, f64)]) -> Result<Envelope> {
        let mut list = sys::envelope_list_new()?;
        for (env, intensity) in terms {
            sys::envelope_list_push(list.pin_mut(), env.raw(), *intensity)?;
        }
        let raw = sys::envelope_linear_combination(list.as_ref().expect("non-null"))?;
        Ok(Envelope::adopt(raw, None))
    }

    /// Number of peaks.
    pub fn len(&self) -> usize {
        sys::envelope_len(self.raw())
    }

    /// Whether there are no peaks.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The masses of the peaks.
    pub fn masses(&self) -> &[f64] {
        sys::envelope_masses(self.raw())
    }

    /// The probabilities of the peaks (molecule counts, for a stochastic
    /// envelope).
    pub fn probs(&self) -> &[f64] {
        sys::envelope_probs(self.raw())
    }

    /// `(mass, prob)` of each peak.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (f64, f64)> + '_ {
        self.masses()
            .iter()
            .copied()
            .zip(self.probs().iter().copied())
    }

    /// Whether the peaks carry configurations.
    pub fn has_confs(&self) -> bool {
        self.isotope_numbers.is_some()
    }

    /// The configuration of peak `index`, if the envelope keeps them. Panics
    /// if `index >= self.len()`.
    pub fn conf(&self, index: usize) -> Option<ConfRef<'_>> {
        assert!(index < self.len(), "peak index out of range");
        let layout = self.isotope_numbers.as_ref()?;
        let all_dim = sys::envelope_all_dim(self.raw()) as usize;
        let confs = sys::envelope_confs(self.raw());
        Some(ConfRef::new(
            &confs[index * all_dim..(index + 1) * all_dim],
            layout,
        ))
    }

    /// The configurations of all peaks, in peak order, if kept.
    pub fn confs(&self) -> Option<impl ExactSizeIterator<Item = ConfRef<'_>>> {
        let layout = self.isotope_numbers.as_ref()?;
        let all_dim = sys::envelope_all_dim(self.raw()) as usize;
        let confs = sys::envelope_confs(self.raw());
        Some(
            (0..self.len())
                .map(move |ii| ConfRef::new(&confs[ii * all_dim..(ii + 1) * all_dim], layout)),
        )
    }

    /// Copies the masses and probabilities out.
    pub fn to_vecs(&self) -> (Vec<f64>, Vec<f64>) {
        (self.masses().to_vec(), self.probs().to_vec())
    }

    /// Sorts the peaks by increasing mass (configurations along).
    pub fn sort_by_mass(&mut self) {
        sys::envelope_sort_by_mass(self.inner.pin_mut()).expect("sorting an envelope");
    }

    /// Sorts the peaks by increasing probability (configurations along).
    pub fn sort_by_prob(&mut self) {
        sys::envelope_sort_by_prob(self.inner.pin_mut()).expect("sorting an envelope");
    }

    /// Sum of the probabilities (computed once, then remembered).
    pub fn total_prob(&mut self) -> f64 {
        sys::envelope_total_prob_of(self.inner.pin_mut())
    }

    /// Multiplies every probability by `factor` (finite, nonnegative).
    pub fn scale(&mut self, factor: f64) -> Result<()> {
        sys::envelope_scale(self.inner.pin_mut(), factor)?;
        Ok(())
    }

    /// Scales the probabilities to sum to 1.
    pub fn normalize(&mut self) -> Result<()> {
        sys::envelope_normalize(self.inner.pin_mut())?;
        Ok(())
    }

    /// Adds `shift` to every mass.
    pub fn shift_mass(&mut self, shift: f64) -> Result<()> {
        sys::envelope_shift_mass(self.inner.pin_mut(), shift)?;
        Ok(())
    }

    /// Replaces the probabilities by counts of `ionic_current` molecules
    /// sampled from them (a simulated measurement). `beta_bias` tunes the
    /// algorithm's speed only; `seed` makes it reproducible.
    pub fn resample(
        &mut self,
        ionic_current: usize,
        beta_bias: f64,
        seed: Option<u64>,
    ) -> Result<()> {
        sys::envelope_resample(
            self.inner.pin_mut(),
            ionic_current,
            beta_bias,
            seed.is_some(),
            seed.unwrap_or(0),
        )?;
        Ok(())
    }

    /// Probability-weighted mean mass of the peaks.
    pub fn average_mass(&mut self) -> Result<f64> {
        Ok(sys::envelope_average_mass(self.inner.pin_mut())?)
    }

    /// Probability-weighted variance of the peaks' masses.
    pub fn variance(&mut self) -> Result<f64> {
        Ok(sys::envelope_variance(self.inner.pin_mut())?)
    }

    /// Square root of [`variance`](Self::variance).
    pub fn stddev(&mut self) -> Result<f64> {
        Ok(self.variance()?.sqrt())
    }

    /// The Wasserstein (earth mover's) distance to `other`; both must be
    /// normalised. Sorts both by mass.
    pub fn wasserstein(&mut self, other: &mut Envelope) -> Result<f64> {
        Ok(sys::envelope_wasserstein(
            self.inner.pin_mut(),
            other.inner.pin_mut(),
        )?)
    }

    /// The signed Wasserstein distance: positive when `other` lies to the
    /// heavier side. Sorts both by mass.
    pub fn oriented_wasserstein(&mut self, other: &mut Envelope) -> Result<f64> {
        Ok(sys::envelope_oriented_wasserstein(
            self.inner.pin_mut(),
            other.inner.pin_mut(),
        )?)
    }

    /// The Wasserstein distance with an abyss: probability may also be moved
    /// to or from an abyss at cost `abyss_depth` per unit instead of being
    /// transported, so the envelopes need not have equal total probability.
    /// `other`'s probabilities are multiplied by `other_scale` first.
    pub fn abyssal_wasserstein(
        &mut self,
        other: &mut Envelope,
        abyss_depth: f64,
        other_scale: f64,
    ) -> Result<f64> {
        Ok(sys::envelope_abyssal_wasserstein(
            self.inner.pin_mut(),
            other.inner.pin_mut(),
            abyss_depth,
            other_scale,
        )?)
    }

    /// Matches probability between the envelopes, peaks at most
    /// `flow_distance` apart, `other`'s probabilities multiplied by
    /// `other_scale`; returns what was matched and what was left over.
    pub fn wasserstein_match(
        &mut self,
        other: &mut Envelope,
        flow_distance: f64,
        other_scale: f64,
    ) -> Result<Match> {
        let r = sys::envelope_wasserstein_match(
            self.inner.pin_mut(),
            other.inner.pin_mut(),
            flow_distance,
            other_scale,
        )?;
        Ok(Match {
            unmatched_this: r.unmatched_this,
            unmatched_other: r.unmatched_other,
            flow: r.flow,
        })
    }

    /// The peaks summed into bins of width `bin_width`, one bin centred on
    /// `middle`; each bin's mass is its centre.
    pub fn bin(&mut self, bin_width: f64, middle: f64) -> Result<Envelope> {
        Ok(Envelope::adopt(
            sys::envelope_bin(self.inner.pin_mut(), bin_width, middle)?,
            None,
        ))
    }

    /// Both envelopes' peaks together (a mixture, unnormalised).
    pub fn try_add(&self, other: &Envelope) -> Result<Envelope> {
        Ok(Envelope::adopt(
            sys::envelope_add(self.raw(), other.raw())?,
            None,
        ))
    }

    /// The convolution: the distribution of the sum of the two masses, as for
    /// a molecule made of two parts.
    pub fn convolve(&self, other: &Envelope) -> Result<Envelope> {
        Ok(Envelope::adopt(
            sys::envelope_convolve(self.raw(), other.raw())?,
            None,
        ))
    }

    /// A deep copy; as [`Clone::clone`], with allocation failure reported.
    pub fn try_clone(&self) -> Result<Envelope> {
        Ok(Envelope::adopt(
            sys::envelope_clone(self.raw())?,
            self.isotope_numbers.as_ref(),
        ))
    }
}

impl Clone for Envelope {
    /// Panics if allocation fails; see [`Envelope::try_clone`].
    fn clone(&self) -> Envelope {
        self.try_clone().expect("cloning an Envelope")
    }
}

impl Default for Envelope {
    /// An envelope with no peaks.
    fn default() -> Envelope {
        Envelope::adopt(
            sys::envelope_empty().expect("allocating an empty envelope"),
            None,
        )
    }
}

/// [`Envelope::try_add`]; panics where that returns an error.
impl Add for &Envelope {
    type Output = Envelope;
    fn add(self, other: &Envelope) -> Envelope {
        self.try_add(other).expect("adding envelopes")
    }
}

/// [`Envelope::convolve`]; panics where that returns an error.
impl Mul for &Envelope {
    type Output = Envelope;
    fn mul(self, other: &Envelope) -> Envelope {
        self.convolve(other).expect("convolving envelopes")
    }
}

impl fmt::Debug for Envelope {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Envelope")
            .field("len", &self.len())
            .field("has_confs", &self.has_confs())
            .finish_non_exhaustive()
    }
}
