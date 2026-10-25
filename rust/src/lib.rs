//! IsoSpec: fine-structure isotopic distributions.

mod ffi;

use ffi::bridge as sys;

/// A molecule.
pub struct Iso {
    inner: cxx::UniquePtr<sys::Iso>,
}

impl Iso {
    /// Parses a formula such as `"C6H12O6"`.
    pub fn from_formula(formula: &str) -> Result<Self, cxx::Exception> {
        Ok(Iso {
            inner: sys::iso_from_formula(formula, false)?,
        })
    }

    /// Mass of the monoisotopic peak.
    pub fn monoisotopic_mass(&self) -> f64 {
        self.inner.getMonoisotopicPeakMass()
    }
}

/// The SIMD level the kernels run at on this machine.
pub fn simd_level() -> &'static str {
    sys::simd_level()
}

#[cfg(test)]
mod tests {
    #[test]
    fn glucose() {
        let iso = super::Iso::from_formula("C6H12O6").unwrap();
        assert!((iso.monoisotopic_mass() - 180.0633881).abs() < 1e-6);
        assert!(super::Iso::from_formula("C6H12O").is_err());
        eprintln!("simd: {}", super::simd_level());
    }
}
