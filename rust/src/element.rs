use std::fmt;
use std::ops::Range;

use crate::ffi::bridge as sys;

/// A chemical element of IsoSpec's built-in table, with its stable isotopes'
/// masses and natural abundances.
///
/// ```
/// use isospec::Element;
/// let carbon = Element::from_symbol("C").unwrap();
/// assert_eq!(carbon.atomic_number(), 6);
/// assert_eq!(carbon.isotopes().count(), 2);
/// assert!(Element::from_symbol("Xx").is_none());
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Element {
    // The element's first row in the C++ library's flat isotope table; only
    // ever obtained from the library, so always a valid first row.
    first: u32,
}

/// One isotope of an [`Element`].
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct Isotope {
    /// Mass in daltons.
    pub mass: f64,
    /// Natural abundance, between 0 and 1.
    pub abundance: f64,
    /// Nucleon number.
    pub mass_number: u32,
}

impl Element {
    /// The element with this symbol (`"C"`, `"Se"`, ...), if the table has it.
    ///
    /// ```
    /// use isospec::Element;
    /// let se = Element::from_symbol("Se").unwrap();
    /// assert_eq!((se.symbol(), se.name(), se.atomic_number()), ("Se", "selenium", 34));
    /// assert_eq!(se.isotopes().len(), 6);
    /// assert!(Element::from_symbol("Xx").is_none());
    /// ```
    pub fn from_symbol(symbol: &str) -> Option<Element> {
        let first = sys::element_find(symbol);
        Element::from_first_index(first)
    }

    pub(crate) fn from_first_index(first: i32) -> Option<Element> {
        if first < 0 || sys::element_isotope_count(first) == 0 {
            return None;
        }
        Some(Element {
            first: first as u32,
        })
    }

    pub(crate) fn first_index(self) -> i32 {
        self.first as i32
    }

    fn rows(self) -> Range<usize> {
        let first = self.first as usize;
        first..first + sys::element_isotope_count(self.first as i32) as usize
    }

    /// Every element in the table, in table order (by atomic number).
    ///
    /// ```
    /// use isospec::Element;
    /// let first = Element::all().next().unwrap();
    /// assert_eq!(first.symbol(), "H");
    /// assert!(Element::all().any(|e| e.symbol() == "U"));
    /// ```
    pub fn all() -> impl Iterator<Item = Element> {
        let len = sys::element_table_len();
        let mut row = 0usize;
        std::iter::from_fn(move || {
            if row >= len {
                return None;
            }
            let element = Element { first: row as u32 };
            row += sys::element_isotope_count(row as i32).max(1) as usize;
            Some(element)
        })
    }

    /// The element's symbol, e.g. `"Se"`.
    pub fn symbol(self) -> &'static str {
        sys::element_table_symbol(self.first as usize)
    }

    /// The element's English name, e.g. `"selenium"`.
    pub fn name(self) -> &'static str {
        sys::element_table_name(self.first as usize)
    }

    /// The atomic number.
    pub fn atomic_number(self) -> u32 {
        sys::element_table_atomic_number(self.first as usize) as u32
    }

    /// The stable isotopes, lightest first.
    ///
    /// ```
    /// use isospec::Element;
    /// let carbon = Element::from_symbol("C").unwrap();
    /// let isotopes: Vec<_> = carbon.isotopes().collect();
    /// assert_eq!(isotopes[0].mass_number, 12);
    /// assert_eq!(isotopes[0].mass, 12.0);
    /// assert!((isotopes[1].abundance - 0.0107).abs() < 1e-3);
    /// ```
    pub fn isotopes(self) -> impl ExactSizeIterator<Item = Isotope> {
        self.rows().map(|row| Isotope {
            mass: sys::element_table_mass(row),
            abundance: sys::element_table_probability(row),
            mass_number: sys::element_table_mass_number(row) as u32,
        })
    }
}

impl fmt::Debug for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Element({})", self.symbol())
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}
