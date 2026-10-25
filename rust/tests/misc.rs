// The smaller surfaces: Debug output, error details, the remaining error paths.

use std::collections::HashSet;

use isospec::{
    Composition, Element, ElementSpec, Envelope, ErrorKind, Iso, LayeredGenerator,
    OrderedGenerator, StochasticGenerator, Threshold, ThresholdGenerator, UnimodTable,
};

#[test]
fn debug_output() {
    let iso = Iso::from_formula("H2O1").unwrap();
    let dbg = format!("{iso:?}");
    assert!(dbg.contains("isotope_numbers: [2, 3]"), "{dbg}");
    let env = iso.envelope().with_confs().total_prob(0.9).unwrap();
    assert_eq!(
        format!("{env:?}"),
        format!("Envelope {{ len: {}, has_confs: true, .. }}", env.len())
    );
    assert_eq!(
        format!("{:?}", Element::from_symbol("O").unwrap()),
        "Element(O)"
    );
    assert!(format!("{:?}", UnimodTable::embedded()).starts_with("UnimodTable"));
    let g = ThresholdGenerator::new(iso.clone(), 0.1, Threshold::Absolute).unwrap();
    assert!(format!("{g:?}").starts_with("ThresholdGenerator"));
    let peaks = LayeredGenerator::new(iso.clone()).unwrap().into_iter();
    assert!(format!("{peaks:?}").starts_with("Peaks { generator: LayeredGenerator"));
    assert!(
        format!("{:?}", OrderedGenerator::new(iso.clone()).unwrap())
            .starts_with("OrderedGenerator")
    );
    assert!(
        format!("{:?}", StochasticGenerator::new(iso, 10).unwrap())
            .starts_with("StochasticGenerator")
    );
}

#[test]
fn errors_carry_kind_and_message() {
    let err = Iso::from_formula("Xx1").unwrap_err();
    assert_eq!(err.kind(), ErrorKind::InvalidArgument);
    assert_eq!(err.message(), err.to_string());
    assert!(!err.message().is_empty());
    let boxed: Box<dyn std::error::Error + Send + Sync> = Box::new(err.clone());
    assert_eq!(boxed.to_string(), err.to_string());
}

#[test]
fn compositions_hash_like_they_compare() {
    let mut set = HashSet::new();
    set.insert(Composition::from_formula("C2H6O1").unwrap());
    assert!(set.contains(&Composition::from_formula("O1H6C2N0").unwrap()));
    assert!(!set.contains(&Composition::from_formula("C2H6O2").unwrap()));
}

#[test]
fn empty_molecules_and_configurations() {
    let empty = Iso::from_isotopes(&[]).unwrap();
    assert_eq!(empty.dim_number(), 0);
    let conf = empty.monoisotopic_conf();
    assert!(conf.is_empty() && conf.as_ref().is_empty());
    let full = Iso::from_formula("H2O1").unwrap().monoisotopic_conf();
    assert!(!full.is_empty() && !full.as_ref().is_empty());
    assert_eq!(full.as_ref().element(1), &[1, 0, 0]);
    assert_eq!(full.as_ref().to_owned(), full);
}

#[test]
fn remaining_error_paths() {
    let iso = Iso::from_formula("C6H12O6").unwrap();
    let mut a = Envelope::from_total_prob(&iso, 0.99).unwrap();
    let mut b = a.clone();
    assert_eq!(
        a.resample(100, -1.0, Some(1)).unwrap_err().kind(),
        ErrorKind::InvalidArgument
    );
    assert!(a.oriented_wasserstein(&mut b).is_ok());
    assert!(a.abyssal_wasserstein(&mut b, f64::NAN, 1.0).is_err());
    assert!(a.wasserstein_match(&mut b, -1.0, 1.0).is_err());
    let mut scaled = a.clone();
    scaled.scale(3.0).unwrap();
    assert!(a.oriented_wasserstein(&mut scaled).is_err());

    let mut iso = Iso::from_formula("C1").unwrap();
    let nan = ElementSpec::new(1, &[f64::NAN], &[1.0]);
    assert_eq!(
        iso.add_element(&nan).unwrap_err().kind(),
        ErrorKind::InvalidArgument
    );
    let big = ElementSpec::new(u32::MAX, &[1.0], &[1.0]);
    assert_eq!(
        iso.add_element(&big).unwrap_err().kind(),
        ErrorKind::OutOfRange
    );
    assert_eq!(iso.dim_number(), 1);

    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        let path = std::ffi::OsStr::from_bytes(b"/tmp/\xff.csv");
        assert_eq!(
            UnimodTable::from_path(path).unwrap_err().kind(),
            ErrorKind::InvalidArgument
        );
    }
}
