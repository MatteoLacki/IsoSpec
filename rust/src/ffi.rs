//! The cxx bridge. Everything here is a raw binding; the safe API is built on
//! it elsewhere in the crate.

#[cxx::bridge(namespace = "isospec_rs")]
pub(crate) mod bridge {
    unsafe extern "C++" {
        include!("isospec/rust/shim/shim.h");

        #[namespace = "IsoSpec"]
        type Iso;

        fn iso_from_formula(formula: &str, use_nominal_masses: bool) -> Result<UniquePtr<Iso>>;
        fn getMonoisotopicPeakMass(self: &Iso) -> f64;

        fn simd_level() -> &'static str;
    }
}
