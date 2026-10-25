use std::fmt;
use std::path::Path;
use std::sync::OnceLock;

use crate::Composition;
use crate::error::{Error, Result};
use crate::ffi::bridge as sys;

/// A table of [Unimod](https://www.unimod.org) modifications: id → name,
/// monoisotopic mass delta and composition delta. Used to resolve
/// `[UNIMOD:<id>]` in peptide sequences.
///
/// ```
/// use isospec::UnimodTable;
/// let oxidation = UnimodTable::embedded().lookup(35).unwrap();
/// assert_eq!(oxidation.name, "Oxidation");
/// assert_eq!(oxidation.composition.to_string(), "O1");
/// ```
pub struct UnimodTable {
    inner: Table,
}

enum Table {
    // The embedded table and those loaded by path live for the whole process
    // in the C++ library.
    Static(&'static sys::UnimodTable),
    Owned(cxx::UniquePtr<sys::UnimodTable>),
}

/// One entry of a [`UnimodTable`].
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct UnimodEntry {
    /// The Unimod accession number.
    pub id: u32,
    /// The Unimod name, e.g. `"Oxidation"`.
    pub name: String,
    /// Monoisotopic mass change, in daltons.
    pub mono_mass: f64,
    /// Composition change (counts may be negative).
    pub composition: Composition,
}

impl UnimodTable {
    /// The table embedded in the library at compile time (from Unimod's own
    /// export, isotope-labelled and glycan "brick" entries excluded).
    pub fn embedded() -> &'static UnimodTable {
        static EMBEDDED: OnceLock<UnimodTable> = OnceLock::new();
        EMBEDDED.get_or_init(|| UnimodTable {
            inner: Table::Static(
                sys::unimod_embedded().expect("parsing the embedded Unimod table"),
            ),
        })
    }

    /// Parses a CSV of rows `id,name,mono_mass,composition`, the composition
    /// a formula with signed counts (the shape of the library's
    /// `data/unimod.csv`).
    ///
    /// ```
    /// use isospec::{Composition, PeptideOptions, UnimodTable};
    /// let table = UnimodTable::from_csv("id,name,mono_mass,composition\n1,Acetyl,42.010565,H2C2O1\n")?;
    /// assert_eq!(table.lookup(1).unwrap().name, "Acetyl");
    /// assert_eq!(table.entries().count(), 1);
    /// let options = PeptideOptions::default().with_mods(&table);
    /// let acetylated = Composition::from_peptide("[UNIMOD:1]-PEPTIDE", &options)?;
    /// assert_eq!(acetylated, "C36H55N7O16".parse()?); // PEPTIDE plus C2H2O1
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_csv(csv: &str) -> Result<UnimodTable> {
        Ok(UnimodTable {
            inner: Table::Owned(sys::unimod_from_csv(csv)?),
        })
    }

    /// Loads a CSV file as [`from_csv`](Self::from_csv) does. Each distinct
    /// path is read once and kept for the life of the process: a file changed
    /// after its first load is not read again.
    ///
    /// ```
    /// use isospec::UnimodTable;
    /// assert!(UnimodTable::from_path("/nonexistent/unimod.csv").is_err());
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_path(path: impl AsRef<Path>) -> Result<UnimodTable> {
        let path = path.as_ref();
        let path = path.to_str().ok_or_else(|| {
            Error::invalid(format!(
                "Unimod table path is not UTF-8: {}",
                path.display()
            ))
        })?;
        if path.is_empty() {
            return Err(Error::invalid("empty Unimod table path"));
        }
        Ok(UnimodTable {
            inner: Table::Static(sys::unimod_for_path(path)?),
        })
    }

    pub(crate) fn raw(&self) -> &sys::UnimodTable {
        match &self.inner {
            Table::Static(t) => t,
            Table::Owned(t) => t.as_ref().expect("non-null"),
        }
    }

    /// One more than the largest id the table can hold.
    pub fn capacity(&self) -> usize {
        sys::unimod_size(self.raw())
    }

    /// The entry for `id`, if the table has it.
    ///
    /// ```
    /// use isospec::UnimodTable;
    /// let table = UnimodTable::embedded();
    /// let carbamidomethyl = table.lookup(4).unwrap();
    /// assert_eq!(carbamidomethyl.name, "Carbamidomethyl");
    /// assert!((carbamidomethyl.mono_mass - 57.021464).abs() < 1e-6);
    /// assert!(table.lookup(0).is_none());
    /// assert!(table.entries().count() > 100);
    /// assert!(table.capacity() > 35);
    /// ```
    pub fn lookup(&self, id: u32) -> Option<UnimodEntry> {
        let mut name = String::new();
        let mut mono_mass = 0.0;
        let mut first = Vec::new();
        let mut counts = Vec::new();
        let found = sys::unimod_lookup(
            self.raw(),
            id,
            &mut name,
            &mut mono_mass,
            &mut first,
            &mut counts,
        )
        .expect("allocation failed");
        if !found {
            return None;
        }
        let composition = first
            .iter()
            .zip(counts)
            .map(|(&f, c)| {
                (
                    crate::Element::from_first_index(f).expect("valid element"),
                    c,
                )
            })
            .collect();
        Some(UnimodEntry {
            id,
            name,
            mono_mass,
            composition,
        })
    }

    /// Every entry, by id.
    pub fn entries(&self) -> impl Iterator<Item = UnimodEntry> + '_ {
        (0..self.capacity() as u32).filter_map(|id| self.lookup(id))
    }
}

impl fmt::Debug for UnimodTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("UnimodTable")
            .field("capacity", &self.capacity())
            .finish_non_exhaustive()
    }
}
