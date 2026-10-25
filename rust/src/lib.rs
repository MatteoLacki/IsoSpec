//! Fine-structure isotopic distributions of molecules, computed by
//! [IsoSpec](https://github.com/MatteoLacki/IsoSpec).
//!
//! The C++ library is compiled from source by the build script (a C++20
//! compiler is needed); this crate is a safe API over it.
//!
//! # Overview
//!
//! A molecule is an [`Iso`]: from a formula, a peptide sequence (with Unimod
//! modifications), a [`Composition`] or raw isotope data. Its isotopic
//! distribution is then computed
//!
//! - whole, as an [`Envelope`] (arrays of masses and probabilities): the
//!   isotopologues above a probability threshold, the smallest set covering a
//!   given total probability, a simulated measurement, or binned on a mass
//!   grid; envelopes can be added, convolved, binned and compared
//!   (Wasserstein distances);
//! - or one isotopologue at a time, by a generator: [`ThresholdGenerator`],
//!   [`LayeredGenerator`], [`OrderedGenerator`] (strictly decreasing
//!   probability) and [`StochasticGenerator`]. Generators iterate over
//!   [`Peak`]s (`for peak in generator`), or step one isotopologue at a time
//!   with `advance()`, which also gives each one's configuration.
//!
//! ```
//! use isospec::{Iso, Threshold, ThresholdGenerator};
//!
//! // The smallest set of isotopologues of a peptide covering 99.9%.
//! let iso = Iso::from_formula("C254H377N65O75S6")?;
//! let mut env = iso.envelope().total_prob(0.999)?;
//! env.sort_by_mass();
//! println!("{} peaks, {:.4} Da to {:.4} Da", env.len(), env.masses()[0], env.masses()[env.len() - 1]);
//!
//! // Or streamed, never stored.
//! let n = ThresholdGenerator::new(iso, 1e-3, Threshold::Relative)?.into_iter().count();
//! assert!(n > 0);
//! # Ok::<(), isospec::Error>(())
//! ```
//!
//! # Threads
//!
//! Every type can be sent to another thread. [`Iso`], [`Envelope`] and
//! [`UnimodTable`] can be shared: their `&self` methods may run concurrently.
//! Generators are stateful and are not `Sync`.
//!
//! # Errors
//!
//! Invalid input is an [`Error`], never undefined behaviour: the C++
//! library's exceptions are caught at the boundary, and the arguments the
//! C++ cannot check itself (slice lengths) are checked before calling it.
//! The [`ErrorKind`] names the exception class.

#![deny(missing_docs)]

mod composition;
mod conf;
mod element;
mod envelope;
mod error;
mod ffi;
mod generator;
pub mod isa;
mod iso;
mod unimod;

pub use composition::{Composition, ExpandedComposition, PeptideOptions};
pub use conf::{Conf, ConfRef};
pub use element::{Element, Isotope};
pub use envelope::{Envelope, EnvelopeBuilder, Match};
pub use error::{Error, ErrorKind, Result};
pub use generator::{
    LayeredGenerator, LayeredOptions, OrderedGenerator, Peak, Peaks, StochasticGenerator,
    StochasticOptions, Threshold, ThresholdGenerator, ThresholdOptions,
};
pub use iso::{ElementSpec, Iso, Masses};
pub use unimod::{UnimodEntry, UnimodTable};

/// The version of the bundled IsoSpec C++ library.
pub const LIBRARY_VERSION: &str = env!("ISOSPEC_LIBRARY_VERSION");

/// What the SIMD kernels run as on this machine: see [`isa::simd_level`].
pub fn simd_level() -> &'static str {
    isa::simd_level()
}
