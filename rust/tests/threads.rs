// What the Send/Sync claims rest on: shared molecules and envelopes read from
// many threads at once, generators moved between threads.

use std::thread;

use isospec::{Envelope, Iso, Threshold, ThresholdGenerator, UnimodTable};

fn assert_send_sync<T: Send + Sync>() {}
fn assert_send<T: Send>() {}

#[test]
fn marker_traits() {
    assert_send_sync::<Iso>();
    assert_send_sync::<Envelope>();
    assert_send_sync::<UnimodTable>();
    assert_send_sync::<isospec::Error>();
    assert_send::<ThresholdGenerator>();
    assert_send::<isospec::LayeredGenerator>();
    assert_send::<isospec::OrderedGenerator>();
    assert_send::<isospec::StochasticGenerator>();
}

#[test]
fn shared_iso_and_envelope() {
    let iso = Iso::from_formula("C254H377N65O75S6").unwrap();
    let env = Envelope::from_total_prob(&iso, 0.999).unwrap();
    let expected_len = env.len();
    let expected_mass: f64 = env.iter().map(|(m, p)| m * p).sum();
    thread::scope(|s| {
        for _ in 0..8 {
            s.spawn(|| {
                for _ in 0..20 {
                    let mine = Envelope::from_total_prob(&iso, 0.999).unwrap();
                    assert_eq!(mine.len(), expected_len);
                    let mass: f64 = env.iter().map(|(m, p)| m * p).sum();
                    assert_eq!(mass, expected_mass);
                    assert!(iso.monoisotopic_mass() > 0.0);
                    let _copy = env.clone();
                    let _iso = iso.clone();
                }
            });
        }
    });
}

#[test]
fn generators_move_between_threads() {
    let iso = Iso::from_formula("C100H202").unwrap();
    let mut g = ThresholdGenerator::new(iso, 1e-6, Threshold::Absolute).unwrap();
    let first = g.advance().unwrap().unwrap();
    let (g, rest) = thread::spawn(move || {
        let mut rest = 0;
        while g.advance().unwrap().is_some() {
            rest += 1;
        }
        (g, rest)
    })
    .join()
    .unwrap();
    assert!(first.prob > 0.0 && rest > 0);
    drop(thread::spawn(move || drop(g)).join());
}

#[test]
fn independent_work_on_many_threads() {
    let peptides = ["PEPTIDE", "PEPM[UNIMOD:35]TIDE", "ACDEFGHIK", "LMNPQRSTVWY"];
    let handles: Vec<_> = (0..16)
        .map(|ii| {
            thread::spawn(move || {
                let options = isospec::PeptideOptions::default();
                let iso = Iso::from_peptide(peptides[ii % 4], &options).unwrap();
                let mut env = iso
                    .envelope()
                    .stochastic(
                        10000,
                        isospec::StochasticOptions {
                            seed: Some(ii as u64),
                            ..Default::default()
                        },
                    )
                    .unwrap();
                env.total_prob()
            })
        })
        .collect();
    for h in handles {
        assert!(h.join().unwrap() <= 10000.0);
    }
}
