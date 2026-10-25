use std::fmt;
use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};
use std::str::FromStr;

use crate::error::{Error, Result};
use crate::ffi::bridge as sys;
use crate::unimod::UnimodTable;
use crate::{Element, Masses};

/// Element counts: a molecular formula, or a difference between two (a
/// modification's delta, say), so counts may be negative.
///
/// Each element appears once, in the order first added. A composition turns
/// into a molecule with [`Iso::from_composition`](crate::Iso::from_composition),
/// which rejects negative counts.
///
/// ```
/// use isospec::{Composition, Element};
/// let glucose: Composition = "C6H12O6".parse().unwrap();
/// let water: Composition = "H2O1".parse().unwrap();
/// let anhydro = &glucose - &water;
/// assert_eq!(anhydro.to_string(), "C6H10O5");
/// assert_eq!(anhydro.count(Element::from_symbol("H").unwrap()), 10);
/// ```
///
/// Two compositions are equal when they have the same nonzero counts, in
/// whatever order: `C2H2O1` equals `O1H2C2N0`.
#[derive(Clone, Debug, Default)]
pub struct Composition {
    elements: Vec<Element>,
    counts: Vec<i32>,
}

/// A composition resolved to the flat isotope arrays IsoSpec computes on:
/// one entry per element in `isotope_numbers` and `atom_counts`, and one per
/// isotope, element after element, in `masses` and `probabilities`.
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct ExpandedComposition {
    /// Number of isotopes of each element.
    pub isotope_numbers: Vec<i32>,
    /// Number of atoms of each element.
    pub atom_counts: Vec<i32>,
    /// Isotope masses, element after element.
    pub masses: Vec<f64>,
    /// Isotope abundances, element after element.
    pub probabilities: Vec<f64>,
}

/// How a peptide sequence is read: see [`Composition::from_peptide`] and
/// [`Iso::from_peptide`](crate::Iso::from_peptide). Start from
/// [`default`](Default::default) and change what you need.
///
/// ```
/// use isospec::{Composition, PeptideOptions};
/// let residues = Composition::from_peptide("PEPTIDE", &PeptideOptions::default().with_add_water(false))?;
/// let whole = Composition::from_peptide("PEPTIDE", &PeptideOptions::default())?;
/// assert_eq!(&residues + &Composition::from_formula("H2O1")?, whole);
/// # Ok::<(), isospec::Error>(())
/// ```
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct PeptideOptions<'a> {
    /// Add the terminal H and OH groups (a whole molecule, rather than the
    /// residues alone). Default: `true`.
    pub add_water: bool,
    /// The table `[UNIMOD:<id>]` modifications are looked up in. Default:
    /// the table embedded in the library.
    pub mods: &'a UnimodTable,
    /// Isotope masses, when building an [`Iso`](crate::Iso). Default: exact.
    pub masses: Masses,
}

impl Default for PeptideOptions<'static> {
    fn default() -> Self {
        PeptideOptions {
            add_water: true,
            mods: UnimodTable::embedded(),
            masses: Masses::Exact,
        }
    }
}

impl<'a> PeptideOptions<'a> {
    /// Sets [`add_water`](Self::add_water).
    pub fn with_add_water(mut self, add_water: bool) -> Self {
        self.add_water = add_water;
        self
    }

    /// Sets the Unimod table, [`mods`](Self::mods).
    pub fn with_mods<'b>(self, mods: &'b UnimodTable) -> PeptideOptions<'b> {
        PeptideOptions {
            add_water: self.add_water,
            mods,
            masses: self.masses,
        }
    }

    /// Sets [`masses`](Self::masses).
    pub fn with_masses(mut self, masses: Masses) -> Self {
        self.masses = masses;
        self
    }
}

impl Composition {
    /// An empty composition.
    pub fn new() -> Self {
        Composition::default()
    }

    /// Parses a formula such as `"C6H12O6"`: element symbols each followed by
    /// a count, which may be negative (`"H-2O-1"`), each element at most once.
    /// The formula is exactly the bytes of `formula`.
    pub fn from_formula(formula: &str) -> Result<Self> {
        let mut first = Vec::new();
        let mut counts = Vec::new();
        sys::composition_from_formula(formula, &mut first, &mut counts)?;
        Ok(Composition::from_raw(&first, counts))
    }

    /// Reads a peptide sequence of one-letter amino-acid codes, with Unimod
    /// modifications in the notation `[UNIMOD:<id>]-PEPTIDE` (N-terminal),
    /// `PEPM[UNIMOD:35]TIDE` (after the modified residue) and
    /// `PEPTIDE-[UNIMOD:<id>]` (C-terminal). Characters that are not amino
    /// acids are ignored; a malformed bracket or an id the table lacks is an
    /// error.
    pub fn from_peptide(sequence: &str, options: &PeptideOptions<'_>) -> Result<Self> {
        let mut composition = Composition::new();
        composition.parse_peptide_into(sequence, options)?;
        Ok(composition)
    }

    /// As [`from_peptide`](Self::from_peptide), but reusing this composition's
    /// storage: for loops over many sequences.
    pub fn parse_peptide_into(
        &mut self,
        sequence: &str,
        options: &PeptideOptions<'_>,
    ) -> Result<()> {
        let mut first = Vec::new();
        let mut counts = Vec::new();
        sys::composition_from_peptide(
            sequence,
            options.mods.raw(),
            options.add_water,
            &mut first,
            &mut counts,
        )?;
        self.elements.clear();
        self.elements.extend(first.iter().map(|&f| element_of(f)));
        self.counts = counts;
        Ok(())
    }

    /// Reads a plain amino-acid sequence (FASTA one-letter codes, any case,
    /// including `U` and `J`), ignoring every other character. Modifications
    /// are not understood: in `M[UNIMOD:35]` the brackets, colon and digits
    /// are skipped and the letters read as residues, so use
    /// [`from_peptide`](Self::from_peptide) for those. With `add_water` the
    /// terminal H and OH groups are added.
    pub fn from_fasta(sequence: &str, add_water: bool) -> Self {
        let mut chnosse = [0i32; 6];
        sys::composition_from_fasta(sequence, add_water, &mut chnosse);
        let mut composition = Composition::new();
        for (symbol, count) in ["C", "H", "N", "O", "S", "Se"].iter().zip(chnosse) {
            if count != 0 {
                composition
                    .add_element(Element::from_symbol(symbol).expect("element table"), count);
            }
        }
        composition
    }

    fn from_raw(first: &[i32], counts: Vec<i32>) -> Self {
        Composition {
            elements: first.iter().map(|&f| element_of(f)).collect(),
            counts,
        }
    }

    pub(crate) fn raw_first_indexes(&self) -> Vec<i32> {
        self.elements.iter().map(|e| e.first_index()).collect()
    }

    pub(crate) fn raw_counts(&self) -> &[i32] {
        &self.counts
    }

    /// Adds `count` atoms of `element` (subtracts, if negative).
    pub fn add_element(&mut self, element: Element, count: i32) {
        match self.elements.iter().position(|&e| e == element) {
            Some(ii) => self.counts[ii] += count,
            None => {
                self.elements.push(element);
                self.counts.push(count);
            }
        }
    }

    /// The count of `element`, 0 if absent.
    pub fn count(&self, element: Element) -> i32 {
        self.elements
            .iter()
            .position(|&e| e == element)
            .map_or(0, |ii| self.counts[ii])
    }

    /// `(element, count)` pairs, in order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (Element, i32)> + '_ {
        self.elements
            .iter()
            .copied()
            .zip(self.counts.iter().copied())
    }

    /// Number of elements listed (including any with a count of 0).
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// Whether no element is listed.
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Whether any count is negative.
    pub fn has_negative_counts(&self) -> bool {
        self.counts.iter().any(|&c| c < 0)
    }

    /// Resolves the composition to isotope masses and abundances from the
    /// element table. Negative counts are passed through.
    pub fn expand(&self, masses: Masses) -> Result<ExpandedComposition> {
        let mut out = ExpandedComposition::default();
        sys::composition_expand(
            &self.raw_first_indexes(),
            &self.counts,
            masses == Masses::Nominal,
            &mut out.isotope_numbers,
            &mut out.atom_counts,
            &mut out.masses,
            &mut out.probabilities,
        )?;
        Ok(out)
    }

    fn combine(&self, other: &Composition, sign: i32) -> Composition {
        let mut result = self.clone();
        for (element, count) in other.iter() {
            result.add_element(element, sign * count);
        }
        result
    }
}

impl Composition {
    fn sorted(&self) -> Vec<(Element, i32)> {
        let mut pairs: Vec<_> = self.iter().filter(|&(_, c)| c != 0).collect();
        pairs.sort_unstable();
        pairs
    }
}

impl PartialEq for Composition {
    fn eq(&self, other: &Composition) -> bool {
        self.sorted() == other.sorted()
    }
}

impl Eq for Composition {}

impl std::hash::Hash for Composition {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sorted().hash(state);
    }
}

fn element_of(first: i32) -> Element {
    Element::from_first_index(first).expect("the library returned an invalid element index")
}

impl FromStr for Composition {
    type Err = Error;
    fn from_str(formula: &str) -> Result<Self> {
        Composition::from_formula(formula)
    }
}

/// Writes the formula, e.g. `C6H12O6`, every count explicit (so it parses
/// back), negative ones as `H-2`.
impl fmt::Display for Composition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (element, count) in self.iter() {
            write!(f, "{element}{count}")?;
        }
        Ok(())
    }
}

impl Add for &Composition {
    type Output = Composition;
    fn add(self, other: &Composition) -> Composition {
        self.combine(other, 1)
    }
}

impl Sub for &Composition {
    type Output = Composition;
    fn sub(self, other: &Composition) -> Composition {
        self.combine(other, -1)
    }
}

impl AddAssign<&Composition> for Composition {
    fn add_assign(&mut self, other: &Composition) {
        for (element, count) in other.iter() {
            self.add_element(element, count);
        }
    }
}

impl SubAssign<&Composition> for Composition {
    fn sub_assign(&mut self, other: &Composition) {
        for (element, count) in other.iter() {
            self.add_element(element, -count);
        }
    }
}

impl Mul<i32> for &Composition {
    type Output = Composition;
    fn mul(self, factor: i32) -> Composition {
        Composition {
            elements: self.elements.clone(),
            counts: self.counts.iter().map(|&c| c * factor).collect(),
        }
    }
}

impl Neg for &Composition {
    type Output = Composition;
    fn neg(self) -> Composition {
        self * -1
    }
}

impl FromIterator<(Element, i32)> for Composition {
    fn from_iter<I: IntoIterator<Item = (Element, i32)>>(iter: I) -> Self {
        let mut composition = Composition::new();
        for (element, count) in iter {
            composition.add_element(element, count);
        }
        composition
    }
}
