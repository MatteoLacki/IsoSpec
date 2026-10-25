// Compositions and the Unimod table.

use isospec::{Composition, Element, ErrorKind, Masses, PeptideOptions, UnimodTable};

fn el(symbol: &str) -> Element {
    Element::from_symbol(symbol).unwrap()
}

#[test]
fn formulas_parse_and_print() {
    let c: Composition = "C6H12O6".parse().unwrap();
    assert_eq!(c.len(), 3);
    assert_eq!(c.count(el("C")), 6);
    assert_eq!(c.count(el("N")), 0);
    assert_eq!(c.to_string(), "C6H12O6");
    assert_eq!(
        c.iter().collect::<Vec<_>>(),
        vec![(el("C"), 6), (el("H"), 12), (el("O"), 6)]
    );
    let delta = Composition::from_formula("H-2O-1").unwrap();
    assert!(delta.has_negative_counts());
    assert_eq!(delta.to_string(), "H-2O-1");
    assert_eq!(
        Composition::from_formula(&delta.to_string()).unwrap(),
        delta
    );
    for bad in ["", "H2O", "C1C1", "Zz1", "C1 H1"] {
        assert_eq!(
            Composition::from_formula(bad).unwrap_err().kind(),
            ErrorKind::InvalidArgument,
            "{bad:?}"
        );
    }
}

#[test]
fn arithmetic_and_equality() {
    let glucose: Composition = "C6H12O6".parse().unwrap();
    let water: Composition = "H2O1".parse().unwrap();
    assert_eq!((&glucose - &water).to_string(), "C6H10O5");
    assert_eq!(&(&glucose - &water) + &water, glucose);
    assert_eq!((&water * 3).to_string(), "H6O3");
    assert_eq!((-&water).to_string(), "H-2O-1");
    let mut c = Composition::new();
    assert!(c.is_empty());
    c += &water;
    c -= &water;
    assert_eq!(c.to_string(), "H0O0");
    // Order does not matter for equality.
    assert_eq!(Composition::from_formula("O6C6H12").unwrap(), glucose);
    let collected: Composition = [(el("C"), 1), (el("C"), 2)].into_iter().collect();
    assert_eq!(collected.to_string(), "C3");
}

#[test]
fn expansion_resolves_isotopes() {
    let e = Composition::from_formula("H2O1")
        .unwrap()
        .expand(Masses::Exact)
        .unwrap();
    assert_eq!(e.isotope_numbers, vec![2, 3]);
    assert_eq!(e.atom_counts, vec![2, 1]);
    assert_eq!(e.masses.len(), 5);
    assert_eq!(e.probabilities.len(), 5);
    let n = Composition::from_formula("H2O1")
        .unwrap()
        .expand(Masses::Nominal)
        .unwrap();
    assert_eq!(n.masses, vec![1.0, 2.0, 16.0, 17.0, 18.0]);
}

#[test]
fn peptides() {
    let options = PeptideOptions::default();
    let p = Composition::from_peptide("PEPM[UNIMOD:35]TIDE", &options).unwrap();
    // IsoSpecPy's ParsePeptideSequence gives the residues, {'H': 60, 'C': 39,
    // 'N': 8, 'O': 16, 'S': 1}; the default here adds the terminal water.
    assert_eq!(p, "C39H62N8O17S1".parse().unwrap());
    let residues = options.with_add_water(false);
    assert_eq!(
        Composition::from_peptide("PEPM[UNIMOD:35]TIDE", &residues).unwrap(),
        "C39H60N8O16S1".parse().unwrap()
    );
    let mut reused = Composition::new();
    reused.parse_peptide_into("PEPMTIDE", &options).unwrap();
    assert_eq!(reused, "C39H62N8O16S1".parse().unwrap());
    reused.parse_peptide_into("G", &options).unwrap();
    assert_eq!(reused.count(el("C")), 2);
    let fasta = Composition::from_fasta("PEPMTIDE", true);
    assert_eq!(
        fasta,
        Composition::from_peptide("PEPMTIDE", &options).unwrap()
    );
    // The plain reader does not understand brackets: it skips the brackets,
    // the colon and the digits, and reads the letters as residues.
    let bracketed = Composition::from_fasta("PEPM[UNIMOD:35]TIDE", true);
    assert_eq!(bracketed, Composition::from_fasta("PEPMUNIMODTIDE", true));
    assert_eq!(Composition::from_fasta("", false), Composition::new());
    // Selenocysteine brings selenium.
    assert_eq!(Composition::from_fasta("U", false).count(el("Se")), 1);
    // A modification can make a count negative.
    let g = Composition::from_peptide("G[UNIMOD:11]", &options.with_add_water(false));
    assert!(g.is_ok());
}

#[test]
fn embedded_unimod_table() {
    let table = UnimodTable::embedded();
    assert!(std::ptr::eq(table, UnimodTable::embedded()));
    let ox = table.lookup(35).unwrap();
    assert_eq!(ox.id, 35);
    assert_eq!(ox.name, "Oxidation");
    assert!((ox.mono_mass - 15.994915).abs() < 1e-6);
    assert_eq!(ox.composition, "O1".parse().unwrap());
    assert!(table.lookup(0).is_none());
    assert!(table.lookup(u32::MAX).is_none());
    assert!(table.capacity() > 35);
    assert!(table.entries().count() > 100);
    assert!(
        table
            .entries()
            .all(|e| table.lookup(e.id) == Some(e.clone()))
    );
}

#[test]
fn unimod_tables_from_csv_and_path() {
    let csv =
        "id,name,mono_mass,composition\n1,Acetyl,42.010565,H2C2O1\n5,Carbamyl,43.005814,H1C1N1O1\n";
    let table = UnimodTable::from_csv(csv).unwrap();
    assert_eq!(table.lookup(1).unwrap().name, "Acetyl");
    assert!(table.lookup(35).is_none());
    let options = PeptideOptions::default().with_mods(&table);
    let acetylated = Composition::from_peptide("[UNIMOD:1]-PEPTIDE", &options).unwrap();
    let plain = Composition::from_peptide("PEPTIDE", &options).unwrap();
    assert_eq!((&acetylated - &plain), "C2H2O1".parse().unwrap());
    assert!(Composition::from_peptide("PEPM[UNIMOD:35]", &options).is_err());

    let dir = std::env::temp_dir().join(format!("isospec-rs-unimod-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("mods.csv");
    std::fs::write(&path, csv).unwrap();
    let from_path = UnimodTable::from_path(&path).unwrap();
    assert_eq!(from_path.lookup(5).unwrap().name, "Carbamyl");
    std::fs::remove_dir_all(&dir).unwrap();
    // Read once and kept: still there with the file gone.
    assert_eq!(
        UnimodTable::from_path(&path)
            .unwrap()
            .lookup(5)
            .unwrap()
            .name,
        "Carbamyl"
    );

    assert_eq!(
        UnimodTable::from_csv("id,name\nnonsense")
            .unwrap_err()
            .kind(),
        ErrorKind::InvalidArgument
    );
    assert!(UnimodTable::from_path(dir.join("missing.csv")).is_err());
    assert_eq!(
        UnimodTable::from_path("").unwrap_err().kind(),
        ErrorKind::InvalidArgument
    );
}
