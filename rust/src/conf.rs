use std::fmt;
use std::sync::Arc;

/// An isotopologue's configuration, borrowed: for each element of the
/// molecule, how many atoms are of each of its isotopes (in the order of the
/// element's isotopes, lightest first for table elements).
///
/// The counts are stored flat, element after element; [`element`](Self::element)
/// and [`elements`](Self::elements) split them.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct ConfRef<'a> {
    counts: &'a [i32],
    isotope_numbers: &'a [i32],
}

/// An isotopologue's configuration, owned. See [`ConfRef`].
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Conf {
    counts: Vec<i32>,
    isotope_numbers: Arc<[i32]>,
}

impl<'a> ConfRef<'a> {
    pub(crate) fn new(counts: &'a [i32], isotope_numbers: &'a [i32]) -> Self {
        debug_assert_eq!(
            counts.len(),
            isotope_numbers.iter().map(|&n| n as usize).sum::<usize>()
        );
        ConfRef {
            counts,
            isotope_numbers,
        }
    }

    /// All counts, element after element.
    pub fn counts(&self) -> &'a [i32] {
        self.counts
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.isotope_numbers.len()
    }

    /// Whether the molecule has no elements.
    pub fn is_empty(&self) -> bool {
        self.isotope_numbers.is_empty()
    }

    /// The isotope counts of element `index` (in the molecule's element
    /// order). Panics if `index >= self.len()`.
    pub fn element(&self, index: usize) -> &'a [i32] {
        let start: usize = self.isotope_numbers[..index]
            .iter()
            .map(|&n| n as usize)
            .sum();
        &self.counts[start..start + self.isotope_numbers[index] as usize]
    }

    /// The isotope counts of each element, in order.
    pub fn elements(self) -> impl ExactSizeIterator<Item = &'a [i32]> + 'a {
        let counts = self.counts;
        let isotope_numbers: &'a [i32] = self.isotope_numbers;
        let mut start = 0usize;
        isotope_numbers.iter().map(move |&n| {
            let part = &counts[start..start + n as usize];
            start += n as usize;
            part
        })
    }

    /// An owned copy.
    pub fn to_owned(&self) -> Conf {
        Conf {
            counts: self.counts.to_vec(),
            isotope_numbers: Arc::from(self.isotope_numbers),
        }
    }
}

impl Conf {
    pub(crate) fn from_parts(counts: Vec<i32>, isotope_numbers: Arc<[i32]>) -> Self {
        Conf {
            counts,
            isotope_numbers,
        }
    }

    /// A borrowed view.
    pub fn as_ref(&self) -> ConfRef<'_> {
        ConfRef::new(&self.counts, &self.isotope_numbers)
    }

    /// All counts, element after element.
    pub fn counts(&self) -> &[i32] {
        &self.counts
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.isotope_numbers.len()
    }

    /// Whether the molecule has no elements.
    pub fn is_empty(&self) -> bool {
        self.isotope_numbers.is_empty()
    }

    /// The isotope counts of element `index`. Panics if `index >= self.len()`.
    pub fn element(&self, index: usize) -> &[i32] {
        self.as_ref().element(index)
    }

    /// The isotope counts of each element, in order.
    pub fn elements(&self) -> impl ExactSizeIterator<Item = &[i32]> {
        self.as_ref().elements()
    }

    /// The counts, as one vector per element.
    pub fn to_nested(&self) -> Vec<Vec<i32>> {
        self.elements().map(<[i32]>::to_vec).collect()
    }
}

impl fmt::Debug for ConfRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.elements()).finish()
    }
}

impl fmt::Debug for Conf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}
