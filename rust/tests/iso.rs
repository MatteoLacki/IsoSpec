// Molecules: construction from every source, the per-molecule statistics, and
// the errors. Expected values come from IsoSpecPy 2.5.0.

mod common;
use common::assert_close;

use isospec::{Composition, Element, ElementSpec, ErrorKind, Iso, Masses, PeptideOptions};

#[test]
fn glucose_statistics_match_python() {
    let iso = Iso::from_formula("C6H12O6").unwrap();
    assert_close(iso.monoisotopic_mass(), 180.06338810844, 1e-12);
    assert_close(iso.average_mass(), 180.15669542527945, 1e-12);
    assert_close(iso.variance(), 0.11749298880254992, 1e-9);
    assert_close(iso.stddev(), 0.11749298880254992f64.sqrt(), 1e-9);
    assert_close(iso.mode_mass(), 180.06338810844, 1e-12);
    assert_close(iso.mode_lprob(), -0.08108074545533395, 1e-9);
    assert_close(iso.lightest_peak_mass(), 180.06338810844, 1e-12);
    assert_close(iso.lightest_peak_lprob(), -0.08108074545533395, 1e-9);
    assert_close(iso.heaviest_peak_mass(), 210.18430803168002, 1e-12);
    assert_close(iso.heaviest_peak_lprob(), -173.08441265848222, 1e-9);
    assert_close(iso.monoisotopic_lprob(), -0.08108074545533395, 1e-9);
    assert!(iso.unlikeliest_peak_lprob() <= iso.heaviest_peak_lprob());
    assert!(iso.mass_rounding_bound() > 0.0 && iso.mass_rounding_bound() < 1e-9);
    // Oxygen: 18O - 16O.
    assert_close(
        iso.max_atom_mass_span(),
        17.9991596128 - 15.9949146196,
        1e-6,
    );
}

#[test]
fn configurations_of_the_extreme_peaks() {
    let iso = Iso::from_formula("C6H12O6").unwrap();
    assert_eq!(iso.dim_number(), 3);
    assert_eq!(iso.all_dim(), 7);
    assert_eq!(iso.isotope_numbers(), &[2, 2, 3]);
    assert_eq!(iso.atom_counts(), &[6, 12, 6]);
    assert_eq!(
        iso.heaviest_peak_conf().to_nested(),
        vec![vec![0, 6], vec![0, 12], vec![0, 0, 6]]
    );
    assert_eq!(
        iso.monoisotopic_conf().to_nested(),
        vec![vec![6, 0], vec![12, 0], vec![6, 0, 0]]
    );
    assert_eq!(iso.lightest_peak_conf().counts(), &[6, 0, 12, 0, 6, 0, 0]);
    let conf = iso.heaviest_peak_conf();
    assert_eq!(conf.len(), 3);
    assert_eq!(conf.element(2), &[0, 0, 6]);
    assert_eq!(conf.elements().count(), 3);
    assert_eq!(format!("{conf:?}"), "[[0, 6], [0, 12], [0, 0, 6]]");
}

#[test]
fn element_isotopes_are_readable() {
    let iso = Iso::from_formula("C6H12O6").unwrap();
    assert_close(iso.element_masses(0)[1], 13.0033548352, 1e-12);
    let probs: f64 = iso.element_lprobs(2).iter().map(|l| l.exp()).sum();
    assert_close(probs, 1.0, 1e-12);
}

#[test]
#[should_panic(expected = "element index out of range")]
fn element_index_is_checked() {
    Iso::from_formula("C6H12O6").unwrap().element_masses(3);
}

#[test]
fn every_constructor_agrees_on_glucose() {
    let reference = Iso::from_formula("C6H12O6").unwrap().monoisotopic_mass();
    let comp: Composition = "C6H12O6".parse().unwrap();
    let from_comp = Iso::from_composition(&comp, Masses::Exact).unwrap();
    assert_eq!(from_comp.monoisotopic_mass(), reference);

    let expanded = comp.expand(Masses::Exact).unwrap();
    let mut specs = Vec::new();
    let mut at = 0;
    for (ii, &n) in expanded.isotope_numbers.iter().enumerate() {
        let n = n as usize;
        specs.push(ElementSpec::new(
            expanded.atom_counts[ii] as u32,
            &expanded.masses[at..at + n],
            &expanded.probabilities[at..at + n],
        ));
        at += n;
    }
    let from_isotopes = Iso::from_isotopes(&specs).unwrap();
    assert_close(from_isotopes.monoisotopic_mass(), reference, 1e-14);

    let parsed: Iso = "C6H12O6".parse().unwrap();
    assert_eq!(parsed.monoisotopic_mass(), reference);
    let cloned = parsed.clone();
    assert_eq!(cloned.monoisotopic_mass(), reference);
    assert_eq!(cloned.isotope_numbers(), parsed.isotope_numbers());
}

#[test]
fn nominal_masses_are_nucleon_numbers() {
    let iso = Iso::from_formula_with("C6H12O6", Masses::Nominal).unwrap();
    assert_eq!(iso.monoisotopic_mass(), 180.0);
    assert_eq!(
        iso.heaviest_peak_mass(),
        6.0 * 13.0 + 12.0 * 2.0 + 6.0 * 18.0
    );
}

#[test]
fn peptides_with_and_without_modifications() {
    let options = PeptideOptions::default();
    let oxidised = Iso::from_peptide("PEPM[UNIMOD:35]TIDE", &options).unwrap();
    // IsoSpecPy gives the residues C39H60N8O16S1, 928.3847990671001 Da; plus water.
    assert_close(
        oxidised.monoisotopic_mass(),
        928.3847990671001 + 18.0105646863,
        1e-10,
    );
    assert_eq!(
        oxidised.monoisotopic_mass(),
        Iso::from_formula("C39H62N8O17S1")
            .unwrap()
            .monoisotopic_mass()
    );
    let plain = Iso::from_peptide("PEPMTIDE", &options).unwrap();
    assert_close(
        oxidised.monoisotopic_mass() - plain.monoisotopic_mass(),
        15.9949146196,
        1e-8,
    );
    // The mods-unaware reader ignores what it does not know, and agrees on
    // the plain sequence.
    let fasta = Iso::from_fasta("PEPMTIDE", true, Masses::Exact).unwrap();
    assert_eq!(fasta.monoisotopic_mass(), plain.monoisotopic_mass());
    let no_water = Iso::from_fasta("PEPMTIDE", false, Masses::Exact).unwrap();
    assert_close(
        plain.monoisotopic_mass() - no_water.monoisotopic_mass(),
        18.0105646863,
        1e-8,
    );
    let residues = Iso::from_peptide("PEPMTIDE", &options.with_add_water(false)).unwrap();
    assert_eq!(residues.monoisotopic_mass(), no_water.monoisotopic_mass());
}

#[test]
fn labelled_element_added_to_a_formula() {
    let mut iso = Iso::from_formula("C50H80N12O15").unwrap();
    let before = iso.monoisotopic_mass();
    iso.add_element(&ElementSpec::new(6, &[12.0, 13.0033548352], &[0.01, 0.99]))
        .unwrap();
    assert_eq!(iso.dim_number(), 5);
    assert_eq!(iso.isotope_numbers()[4], 2);
    // The new element's most abundant isotope is 13C.
    assert_close(iso.monoisotopic_mass() - before, 6.0 * 13.0033548352, 1e-10);
}

#[test]
fn bad_input_is_an_error_of_the_right_kind() {
    let kind = |r: isospec::Result<Iso>| r.unwrap_err().kind();
    assert_eq!(kind(Iso::from_formula("")), ErrorKind::InvalidArgument);
    assert_eq!(kind(Iso::from_formula("H2O")), ErrorKind::InvalidArgument);
    assert_eq!(kind(Iso::from_formula("Xx1")), ErrorKind::InvalidArgument);
    assert_eq!(kind(Iso::from_formula("C1C1")), ErrorKind::InvalidArgument);
    assert_eq!(kind(Iso::from_formula("H-2O1")), ErrorKind::InvalidArgument);
    assert_eq!(
        kind(Iso::from_formula("C1\0H1")),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Iso::from_formula("C99999999999")),
        ErrorKind::OutOfRange
    );
    assert_eq!(kind(Iso::from_formula("C100000000")), ErrorKind::Length);
    assert_eq!(
        kind(Iso::from_peptide(
            "PEPM[UNIMOD:999999]",
            &PeptideOptions::default()
        )),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Iso::from_peptide("PEP[UNIMOD:", &PeptideOptions::default())),
        ErrorKind::InvalidArgument
    );
    let spec = |masses: &'static [f64], probabilities: &'static [f64]| {
        ElementSpec::new(3, masses, probabilities)
    };
    assert_eq!(
        kind(Iso::from_isotopes(&[spec(&[1.0, 2.0], &[1.0])])),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Iso::from_isotopes(&[spec(&[], &[])])),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Iso::from_isotopes(&[spec(&[1.0, 2.0], &[0.5, 0.2])])),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Iso::from_isotopes(&[spec(&[f64::NAN], &[1.0])])),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Iso::from_isotopes(&[spec(&[1.0], &[f64::NAN])])),
        ErrorKind::InvalidArgument
    );
    let huge = ElementSpec::new(u32::MAX, &[1.0], &[1.0]);
    assert_eq!(kind(Iso::from_isotopes(&[huge])), ErrorKind::OutOfRange);
    let negative: Composition = "H-2".parse().unwrap();
    assert_eq!(
        kind(Iso::from_composition(&negative, Masses::Exact)),
        ErrorKind::InvalidArgument
    );
    let mut iso = Iso::from_formula("C1").unwrap();
    let err = iso
        .add_element(&ElementSpec::new(1, &[1.0], &[]))
        .unwrap_err();
    assert_eq!(err.kind(), ErrorKind::InvalidArgument);
    assert_eq!(
        iso.dim_number(),
        1,
        "a failed add_element leaves the molecule as it was"
    );
    // Messages come through.
    let err = Iso::from_formula("H2O").unwrap_err();
    assert!(err.to_string().contains("H2O1"), "{err}");
}

#[test]
fn marginal_size_estimates_grow_with_the_target() {
    let iso = Iso::from_formula("C100H202S5").unwrap();
    let small = iso.marginal_log_size_estimates(0.5).unwrap();
    let large = iso.marginal_log_size_estimates(0.999).unwrap();
    assert_eq!(small.len(), 3);
    for (a, b) in small.iter().zip(&large) {
        assert!(b >= a);
    }
}

#[test]
fn elements_and_their_isotopes() {
    let selenium = Element::from_symbol("Se").unwrap();
    assert_eq!(selenium.symbol(), "Se");
    assert_eq!(selenium.atomic_number(), 34);
    assert_eq!(selenium.isotopes().len(), 6);
    let total: f64 = selenium.isotopes().map(|i| i.abundance).sum();
    assert_close(total, 1.0, 1e-3);
    assert_eq!(selenium.to_string(), "Se");
    assert!(!selenium.name().is_empty());
    let carbon = Element::from_symbol("C").unwrap();
    assert_eq!(
        carbon.isotopes().map(|i| i.mass_number).collect::<Vec<_>>(),
        vec![12, 13]
    );
    for bad in ["", "c", "CC", "Cxx", "Q", "\0"] {
        assert!(Element::from_symbol(bad).is_none(), "{bad:?}");
    }
    let all: Vec<Element> = Element::all().collect();
    assert!(all.len() > 80);
    assert_eq!(all[0].symbol(), "H");
    for element in &all {
        assert_eq!(Element::from_symbol(element.symbol()), Some(*element));
    }
}

#[test]
fn library_version_is_the_repository_version() {
    assert_eq!(
        isospec::LIBRARY_VERSION,
        include_str!("../../VERSION").trim()
    );
}

#[test]
fn crate_version_metadata_names_the_bundled_library() {
    // Cargo.toml's version is "<crate semver>+<VERSION>": bumping VERSION
    // without it would publish a crate that misstates what it contains.
    let metadata = env!("CARGO_PKG_VERSION").split_once('+').map(|(_, m)| m);
    assert_eq!(metadata, Some(isospec::LIBRARY_VERSION));
}
