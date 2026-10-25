use std::fmt;
use std::str::FromStr;
use std::sync::Arc;

use cxx::UniquePtr;

use crate::composition::PeptideOptions;
use crate::error::{Error, ErrorKind, Result};
use crate::ffi::bridge as sys;
use crate::{Composition, Conf};

/// Which isotope masses a molecule is built with.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Masses {
    /// The isotopes' real masses.
    #[default]
    Exact,
    /// Nucleon numbers instead of masses (so `12.0` and `13.0` for carbon).
    Nominal,
}

/// One element of a molecule given by its isotopes rather than by symbol:
/// for labelled or otherwise non-natural elements. See
/// [`Iso::from_isotopes`] and [`Iso::add_element`]; made with
/// [`ElementSpec::new`].
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct ElementSpec<'a> {
    /// Number of atoms of the element.
    pub atoms: u32,
    /// Mass of each isotope.
    pub masses: &'a [f64],
    /// Abundance of each isotope, in (0, 1]. They must sum to 1 within 1e-3,
    /// and are divided by their sum.
    pub probabilities: &'a [f64],
}

/// A molecule, as IsoSpec computes on it: for each element, its atom count
/// and its isotopes' masses and abundances.
///
/// Molecules come from a formula, a peptide sequence, a [`Composition`] or
/// raw isotope data. Isotopic distributions are then computed by an
/// [`Envelope`](crate::Envelope) (all at once) or a generator (one
/// configuration at a time).
///
/// ```
/// use isospec::Iso;
/// let glucose = Iso::from_formula("C6H12O6")?;
/// assert!((glucose.monoisotopic_mass() - 180.063388).abs() < 1e-6);
/// # Ok::<(), isospec::Error>(())
/// ```
pub struct Iso {
    inner: UniquePtr<sys::Iso>,
}

impl<'a> ElementSpec<'a> {
    /// `atoms` atoms of an element whose isotopes have these masses and
    /// abundances (checked when the molecule is built).
    ///
    /// ```
    /// use isospec::ElementSpec;
    /// // Six carbons, 99% carbon-13.
    /// let carbon13 = ElementSpec::new(6, &[12.0, 13.0033548352], &[0.01, 0.99]);
    /// assert_eq!(carbon13.atoms, 6);
    /// ```
    pub fn new(atoms: u32, masses: &'a [f64], probabilities: &'a [f64]) -> Self {
        ElementSpec {
            atoms,
            masses,
            probabilities,
        }
    }
}

fn atoms(n: u32) -> Result<i32> {
    i32::try_from(n).map_err(|_| Error::new(ErrorKind::OutOfRange, "atom count out of range"))
}

fn check_isotopes(spec: &ElementSpec<'_>) -> Result<i32> {
    if spec.masses.len() != spec.probabilities.len() {
        return Err(Error::invalid(format!(
            "an element has {} isotope masses but {} probabilities",
            spec.masses.len(),
            spec.probabilities.len()
        )));
    }
    i32::try_from(spec.masses.len()).map_err(|_| Error::new(ErrorKind::Length, "too many isotopes"))
}

impl Iso {
    pub(crate) fn from_raw(inner: UniquePtr<sys::Iso>) -> Iso {
        debug_assert!(!inner.is_null());
        Iso { inner }
    }

    pub(crate) fn raw(&self) -> &sys::Iso {
        self.inner.as_ref().expect("non-null")
    }

    pub(crate) fn into_raw(self) -> UniquePtr<sys::Iso> {
        self.inner
    }

    /// Parses a formula such as `"C6H12O6"`: element symbols, each followed
    /// by its count (`"H2O1"`, not `"H2O"`), each element once.
    ///
    /// ```
    /// use isospec::Iso;
    /// let water = Iso::from_formula("H2O1")?;
    /// assert_eq!(water.dim_number(), 2);
    /// assert!(Iso::from_formula("H2O").is_err()); // every count explicit
    /// let parsed: Iso = "H2O1".parse()?; // the same, through FromStr
    /// assert_eq!(parsed.monoisotopic_mass(), water.monoisotopic_mass());
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_formula(formula: &str) -> Result<Iso> {
        Iso::from_formula_with(formula, Masses::Exact)
    }

    /// As [`from_formula`](Self::from_formula), with the choice of masses.
    ///
    /// ```
    /// use isospec::{Iso, Masses};
    /// let nominal = Iso::from_formula_with("C6H12O6", Masses::Nominal)?;
    /// assert_eq!(nominal.monoisotopic_mass(), 180.0);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_formula_with(formula: &str, masses: Masses) -> Result<Iso> {
        Ok(Iso::from_raw(sys::iso_from_formula(
            formula,
            masses == Masses::Nominal,
        )?))
    }

    /// A peptide from a plain amino-acid sequence (FASTA one-letter codes,
    /// any case, including selenocysteine `U` and `J`); any other character is
    /// ignored. With `add_water` the terminal H and OH groups are included.
    ///
    /// ```
    /// use isospec::{Iso, Masses};
    /// let peptide = Iso::from_fasta("PEPTIDE", true, Masses::Exact)?;
    /// assert!((peptide.monoisotopic_mass() - 799.359964).abs() < 1e-6);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_fasta(sequence: &str, add_water: bool, masses: Masses) -> Result<Iso> {
        Ok(Iso::from_raw(sys::iso_from_fasta(
            sequence,
            masses == Masses::Nominal,
            add_water,
        )?))
    }

    /// A peptide from a sequence with `[UNIMOD:<id>]` modifications: see
    /// [`Composition::from_peptide`] for the notation.
    ///
    /// ```
    /// use isospec::{Iso, PeptideOptions};
    /// let oxidised = Iso::from_peptide("PEPM[UNIMOD:35]TIDE", &PeptideOptions::default())?;
    /// let plain = Iso::from_peptide("PEPMTIDE", &PeptideOptions::default())?;
    /// let delta = oxidised.monoisotopic_mass() - plain.monoisotopic_mass();
    /// assert!((delta - 15.994915).abs() < 1e-6);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_peptide(sequence: &str, options: &PeptideOptions<'_>) -> Result<Iso> {
        Ok(Iso::from_raw(sys::iso_from_peptide(
            sequence,
            options.mods.raw(),
            options.masses == Masses::Nominal,
            options.add_water,
        )?))
    }

    /// A molecule of the elements of `composition`, with the element table's
    /// isotopes. Negative counts are an error.
    ///
    /// ```
    /// use isospec::{Composition, Iso, Masses};
    /// let glucose: Composition = "C6H12O6".parse()?;
    /// let water: Composition = "H2O1".parse()?;
    /// let anhydro = Iso::from_composition(&(&glucose - &water), Masses::Exact)?;
    /// assert!((anhydro.monoisotopic_mass() - 162.052823).abs() < 1e-6);
    /// assert!(Iso::from_composition(&(&water - &glucose), Masses::Exact).is_err()); // negative
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_composition(composition: &Composition, masses: Masses) -> Result<Iso> {
        Ok(Iso::from_raw(sys::iso_from_composition(
            &composition.raw_first_indexes(),
            composition.raw_counts(),
            masses == Masses::Nominal,
        )?))
    }

    /// A molecule given element by element by raw isotope data.
    ///
    /// ```
    /// use isospec::{ElementSpec, Iso};
    /// // Six carbons, 99% carbon-13.
    /// let labelled = Iso::from_isotopes(&[ElementSpec::new(6, &[12.0, 13.0033548352], &[0.01, 0.99])])?;
    /// assert_eq!(labelled.dim_number(), 1);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn from_isotopes(elements: &[ElementSpec<'_>]) -> Result<Iso> {
        let mut isotope_numbers = Vec::with_capacity(elements.len());
        let mut atom_counts = Vec::with_capacity(elements.len());
        let mut masses = Vec::new();
        let mut probs = Vec::new();
        for spec in elements {
            isotope_numbers.push(check_isotopes(spec)?);
            atom_counts.push(atoms(spec.atoms)?);
            masses.extend_from_slice(spec.masses);
            probs.extend_from_slice(spec.probabilities);
        }
        if i32::try_from(elements.len()).is_err() || i32::try_from(masses.len()).is_err() {
            return Err(Error::new(ErrorKind::Length, "too many isotopes in total"));
        }
        Ok(Iso::from_raw(sys::iso_from_isotopes(
            &isotope_numbers,
            &atom_counts,
            &masses,
            &probs,
        )?))
    }

    /// Adds one more element, given by raw isotope data (even one the
    /// molecule already has: it becomes a separate dimension, as for a
    /// labelled subset of a molecule's atoms).
    ///
    /// ```
    /// use isospec::{ElementSpec, Iso};
    /// let mut iso = Iso::from_formula("C50H80N12O15")?;
    /// // Six more carbons, 13C-labelled.
    /// iso.add_element(&ElementSpec::new(6, &[12.0, 13.0033548352], &[0.01, 0.99]))?;
    /// assert_eq!(iso.dim_number(), 5);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn add_element(&mut self, element: &ElementSpec<'_>) -> Result<()> {
        check_isotopes(element)?;
        let atom_count = atoms(element.atoms)?;
        sys::iso_add_element(
            self.inner.pin_mut(),
            atom_count,
            element.masses,
            element.probabilities,
        )?;
        Ok(())
    }

    /// A deep copy; as [`Clone::clone`], with allocation failure reported.
    ///
    /// ```
    /// use isospec::Iso;
    /// let iso = Iso::from_formula("C6H12O6")?;
    /// let copy = iso.try_clone()?;
    /// assert_eq!(copy.monoisotopic_mass(), iso.clone().monoisotopic_mass());
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn try_clone(&self) -> Result<Iso> {
        Ok(Iso::from_raw(sys::iso_clone(self.raw())?))
    }

    /// Number of elements (each a separate dimension of the distribution).
    ///
    /// ```
    /// use isospec::Iso;
    /// let glucose = Iso::from_formula("C6H12O6")?;
    /// assert_eq!(glucose.dim_number(), 3);
    /// assert_eq!(glucose.isotope_numbers(), &[2, 2, 3]); // C, H, O
    /// assert_eq!(glucose.all_dim(), 7);
    /// assert_eq!(glucose.atom_counts(), &[6, 12, 6]);
    /// assert_eq!(glucose.element_masses(0), &[12.0, 13.0033548352]);
    /// assert!(glucose.element_lprobs(0)[0] > glucose.element_lprobs(0)[1]); // 12C is the commoner
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn dim_number(&self) -> usize {
        self.raw().getDimNumber() as usize
    }

    /// Total number of isotopes over all elements: the length of a
    /// configuration's flat counts.
    pub fn all_dim(&self) -> usize {
        self.raw().getAllDim() as usize
    }

    /// Number of isotopes of each element.
    pub fn isotope_numbers(&self) -> &[i32] {
        sys::iso_isotope_numbers(self.raw())
    }

    pub(crate) fn isotope_numbers_arc(&self) -> Arc<[i32]> {
        Arc::from(self.isotope_numbers())
    }

    /// Number of atoms of each element.
    pub fn atom_counts(&self) -> &[i32] {
        sys::iso_atom_counts(self.raw())
    }

    /// Isotope masses of element `element`. Panics if `element >= dim_number()`.
    pub fn element_masses(&self, element: usize) -> &[f64] {
        assert!(element < self.dim_number(), "element index out of range");
        sys::iso_element_masses(self.raw(), element as i32)
    }

    /// Natural logarithms of the isotope abundances of element `element`, as
    /// the molecule uses them (renormalised to sum to 1). Panics if
    /// `element >= dim_number()`.
    pub fn element_lprobs(&self, element: usize) -> &[f64] {
        assert!(element < self.dim_number(), "element index out of range");
        sys::iso_element_lprobs(self.raw(), element as i32)
    }

    fn conf_with(&self, fill: fn(&sys::Iso, &mut [i32])) -> Conf {
        let mut counts = vec![0; self.all_dim()];
        fill(self.raw(), &mut counts);
        Conf::from_parts(counts, self.isotope_numbers_arc())
    }

    /// Mass of the lightest isotopologue.
    ///
    /// ```
    /// use isospec::Iso;
    /// let glucose = Iso::from_formula("C6H12O6")?;
    /// // All 12C, 1H and 16O at one end, all 13C, 2H and 18O at the other.
    /// assert_eq!(glucose.lightest_peak_conf().to_nested(), [vec![6, 0], vec![12, 0], vec![6, 0, 0]]);
    /// assert_eq!(glucose.heaviest_peak_conf().to_nested(), [vec![0, 6], vec![0, 12], vec![0, 0, 6]]);
    /// assert!((glucose.heaviest_peak_mass() - glucose.lightest_peak_mass() - 30.121).abs() < 1e-3);
    /// assert!(glucose.heaviest_peak_lprob() < glucose.lightest_peak_lprob());
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn lightest_peak_mass(&self) -> f64 {
        self.raw().getLightestPeakMass()
    }

    /// Log-probability of the lightest isotopologue.
    pub fn lightest_peak_lprob(&self) -> f64 {
        self.raw().getLightestPeakLProb()
    }

    /// Configuration of the lightest isotopologue.
    pub fn lightest_peak_conf(&self) -> Conf {
        self.conf_with(sys::iso_lightest_conf)
    }

    /// Mass of the heaviest isotopologue.
    pub fn heaviest_peak_mass(&self) -> f64 {
        self.raw().getHeaviestPeakMass()
    }

    /// Log-probability of the heaviest isotopologue.
    pub fn heaviest_peak_lprob(&self) -> f64 {
        self.raw().getHeaviestPeakLProb()
    }

    /// Configuration of the heaviest isotopologue.
    pub fn heaviest_peak_conf(&self) -> Conf {
        self.conf_with(sys::iso_heaviest_conf)
    }

    /// Mass of the monoisotopic isotopologue: every atom its element's most
    /// abundant isotope.
    ///
    /// ```
    /// use isospec::Iso;
    /// let glucose = Iso::from_formula("C6H12O6")?;
    /// assert!((glucose.monoisotopic_mass() - 180.063388).abs() < 1e-6);
    /// assert_eq!(glucose.monoisotopic_conf().to_nested(), [vec![6, 0], vec![12, 0], vec![6, 0, 0]]);
    /// assert!((glucose.monoisotopic_lprob().exp() - 0.922).abs() < 1e-3);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn monoisotopic_mass(&self) -> f64 {
        self.raw().getMonoisotopicPeakMass()
    }

    /// Log-probability of the monoisotopic isotopologue.
    pub fn monoisotopic_lprob(&self) -> f64 {
        self.raw().getMonoisotopicPeakLProb()
    }

    /// Configuration of the monoisotopic isotopologue.
    pub fn monoisotopic_conf(&self) -> Conf {
        self.conf_with(sys::iso_monoisotopic_conf)
    }

    /// Mass of the most probable isotopologue.
    ///
    /// ```
    /// use isospec::Iso;
    /// // Small molecules peak at the monoisotopic mass; large ones do not.
    /// let glucose = Iso::from_formula("C6H12O6")?;
    /// assert_eq!(glucose.mode_mass(), glucose.monoisotopic_mass());
    /// let insulin = Iso::from_formula("C254H377N65O75S6")?;
    /// assert!(insulin.mode_mass() > insulin.monoisotopic_mass() + 2.0);
    /// assert!(insulin.mode_lprob() > insulin.monoisotopic_lprob());
    /// assert!(insulin.unlikeliest_peak_lprob() < insulin.mode_lprob());
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn mode_mass(&self) -> f64 {
        self.raw().getModeMass()
    }

    /// Log-probability of the most probable isotopologue.
    pub fn mode_lprob(&self) -> f64 {
        self.raw().getModeLProb()
    }

    /// Log-probability of the least probable isotopologue.
    pub fn unlikeliest_peak_lprob(&self) -> f64 {
        self.raw().getUnlikeliestPeakLProb()
    }

    /// Mean mass of the distribution.
    ///
    /// ```
    /// use isospec::Iso;
    /// let glucose = Iso::from_formula("C6H12O6")?;
    /// assert!((glucose.average_mass() - 180.156695).abs() < 1e-6);
    /// assert!((glucose.stddev() - glucose.variance().sqrt()).abs() < 1e-15);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn average_mass(&self) -> f64 {
        self.raw().getTheoreticalAverageMass()
    }

    /// Variance of the mass.
    pub fn variance(&self) -> f64 {
        self.raw().variance()
    }

    /// Standard deviation of the mass.
    pub fn stddev(&self) -> f64 {
        self.variance().sqrt()
    }

    /// A bound on how far any isotopologue mass the library computes may be
    /// from the exact sum, from rounding.
    ///
    /// ```
    /// use isospec::Iso;
    /// let glucose = Iso::from_formula("C6H12O6")?;
    /// assert!(glucose.mass_rounding_bound() < 1e-9);
    /// assert!((glucose.max_atom_mass_span() - 2.004245).abs() < 1e-6); // 16O to 18O
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn mass_rounding_bound(&self) -> f64 {
        self.raw().getMassRoundingBound()
    }

    /// The largest mass difference between two isotopes of one element of the
    /// molecule.
    pub fn max_atom_mass_span(&self) -> f64 {
        self.raw().getMaxAtomMassSpan()
    }

    /// For each element, an estimate of the log of the number of its
    /// subisotopologues needed to cover `target_total_prob` of its
    /// distribution (what the generators use to order the elements).
    ///
    /// ```
    /// use isospec::Iso;
    /// let estimates = Iso::from_formula("C6H12O6")?.marginal_log_size_estimates(0.99)?;
    /// assert_eq!(estimates.len(), 3);
    /// // Carbon needs the most subisotopologues.
    /// assert!(estimates[0] > estimates[1] && estimates[0] > estimates[2]);
    /// # Ok::<(), isospec::Error>(())
    /// ```
    pub fn marginal_log_size_estimates(&self, target_total_prob: f64) -> Result<Vec<f64>> {
        let mut out = vec![0.0; self.dim_number()];
        sys::iso_marginal_log_size_estimates(self.raw(), &mut out, target_total_prob)?;
        Ok(out)
    }
}

impl Clone for Iso {
    /// Panics if allocation fails; see [`Iso::try_clone`].
    fn clone(&self) -> Iso {
        self.try_clone().expect("cloning an Iso")
    }
}

impl FromStr for Iso {
    type Err = Error;
    fn from_str(formula: &str) -> Result<Iso> {
        Iso::from_formula(formula)
    }
}

impl fmt::Debug for Iso {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Iso")
            .field("isotope_numbers", &self.isotope_numbers())
            .field("atom_counts", &self.atom_counts())
            .field("monoisotopic_mass", &self.monoisotopic_mass())
            .finish()
    }
}
