// The text entry points of the safe API: formulas, peptide sequences (with
// Unimod brackets), plain FASTA, compositions, Unimod CSVs. Any input gives a
// value or an Err; anything else (a panic, a sanitizer report) is a bug. What
// parses is used a little, to reach the code behind the parsers.
#![no_main]

use isospec::{Composition, Envelope, Iso, Masses, PeptideOptions, UnimodTable};
use libfuzzer_sys::fuzz_target;

fn exercise(iso: &Iso) {
    let _ = (
        iso.monoisotopic_mass(),
        iso.average_mass(),
        iso.variance(),
        iso.mode_lprob(),
    );
    let _ = iso.lightest_peak_conf();
    let _ = iso.marginal_log_size_estimates(0.9);
    // Small molecules only: the work grows with the atom counts.
    if iso.atom_counts().iter().all(|&c| c <= 200) {
        if let Ok(env) = Envelope::from_total_prob(iso, 0.9) {
            assert!(env.probs().iter().all(|p| p.is_finite() && *p >= 0.0));
        }
    }
}

fuzz_target!(|data: &[u8]| {
    let Some((&selector, rest)) = data.split_first() else {
        return;
    };
    let text = String::from_utf8_lossy(rest);
    match selector % 6 {
        0 => {
            if let Ok(iso) = Iso::from_formula(&text) {
                exercise(&iso);
            }
        }
        1 => {
            let options = PeptideOptions {
                add_water: selector & 8 != 0,
                ..PeptideOptions::default()
            };
            if let Ok(iso) = Iso::from_peptide(&text, &options) {
                exercise(&iso);
            }
        }
        2 => {
            let c = Composition::from_fasta(&text, selector & 8 != 0);
            let _ = c.expand(Masses::Exact);
        }
        3 => {
            if let Ok(c) = Composition::from_formula(&text) {
                assert_eq!(
                    Composition::from_formula(&c.to_string()).ok(),
                    Some(c.clone())
                );
                if let Ok(iso) = Iso::from_composition(&c, Masses::Nominal) {
                    exercise(&iso);
                }
            }
        }
        4 => {
            if let Ok(table) = UnimodTable::from_csv(&text) {
                for entry in table.entries() {
                    let _ = entry.composition.to_string();
                }
                let options = PeptideOptions {
                    mods: &table,
                    ..PeptideOptions::default()
                };
                let _ = Composition::from_peptide("PEPM[UNIMOD:35]TIDE[UNIMOD:1]", &options);
            }
        }
        _ => {
            let mut c = Composition::new();
            let options = PeptideOptions::default();
            for part in text.split('\n') {
                let _ = c.parse_peptide_into(part, &options);
            }
        }
    }
});
