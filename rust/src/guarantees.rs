//! What the borrow checker enforces, as doctests: each `compile_fail` block
//! has a twin that compiles, so that a block failing for an unrelated reason
//! (a typo) does not pass as the guarantee.
//!
//! A generator is not `Sync` (it is stateful), though it is `Send`:
//!
//! ```compile_fail
//! fn shared<T: Sync>(_: &T) {}
//! let g = isospec::OrderedGenerator::new(isospec::Iso::from_formula("H2O1").unwrap()).unwrap();
//! shared(&g);
//! ```
//! ```
//! fn sent<T: Send>(_: T) {}
//! let g = isospec::OrderedGenerator::new(isospec::Iso::from_formula("H2O1").unwrap()).unwrap();
//! sent(g);
//! ```
//!
//! An envelope's arrays are borrowed: they cannot outlive an operation that
//! may reorder or reallocate them.
//!
//! ```compile_fail
//! let iso = isospec::Iso::from_formula("H2O1").unwrap();
//! let mut env = isospec::Envelope::from_total_prob(&iso, 0.9).unwrap();
//! let masses = env.masses();
//! env.sort_by_mass();
//! println!("{}", masses[0]);
//! ```
//! ```
//! let iso = isospec::Iso::from_formula("H2O1").unwrap();
//! let mut env = isospec::Envelope::from_total_prob(&iso, 0.9).unwrap();
//! env.sort_by_mass();
//! let masses = env.masses();
//! println!("{}", masses[0]);
//! ```
//!
//! A generator takes its molecule: clone it to keep one.
//!
//! ```compile_fail
//! let iso = isospec::Iso::from_formula("H2O1").unwrap();
//! let g = isospec::OrderedGenerator::new(iso).unwrap();
//! println!("{}", iso.monoisotopic_mass());
//! ```
//! ```
//! let iso = isospec::Iso::from_formula("H2O1").unwrap();
//! let g = isospec::OrderedGenerator::new(iso.clone()).unwrap();
//! println!("{}", iso.monoisotopic_mass());
//! ```
//!
//! A generator's configuration is valid only until it moves on.
//!
//! ```compile_fail
//! let mut g = isospec::OrderedGenerator::new(isospec::Iso::from_formula("H2O1").unwrap()).unwrap();
//! g.advance().unwrap();
//! let conf = g.conf().unwrap();
//! g.advance().unwrap();
//! println!("{:?}", conf.counts());
//! ```
//! ```
//! let mut g = isospec::OrderedGenerator::new(isospec::Iso::from_formula("H2O1").unwrap()).unwrap();
//! g.advance().unwrap();
//! let conf = g.conf_owned().unwrap();
//! g.advance().unwrap();
//! println!("{:?}", conf.counts());
//! ```
//!
//! Two envelopes are compared through `&mut` both: the same one cannot be both.
//!
//! ```compile_fail
//! let iso = isospec::Iso::from_formula("H2O1").unwrap();
//! let mut env = isospec::Envelope::from_total_prob(&iso, 0.9).unwrap();
//! env.wasserstein(&mut env).unwrap();
//! ```
//! ```
//! let iso = isospec::Iso::from_formula("H2O1").unwrap();
//! let mut env = isospec::Envelope::from_total_prob(&iso, 0.9).unwrap();
//! let mut copy = env.clone();
//! env.wasserstein(&mut copy).unwrap();
//! ```
