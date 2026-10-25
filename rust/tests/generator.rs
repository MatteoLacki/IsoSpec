// Generators: what they enumerate, the cursor API, and agreement with the
// envelopes.

mod common;
use common::assert_close;

use isospec::{
    Envelope, Iso, LayeredGenerator, LayeredOptions, OrderedGenerator, StochasticGenerator,
    StochasticOptions, Threshold, ThresholdGenerator, ThresholdOptions,
};

fn iso(formula: &str) -> Iso {
    Iso::from_formula(formula).unwrap()
}

fn sorted_masses(mut v: Vec<f64>) -> Vec<f64> {
    v.sort_by(f64::total_cmp);
    v
}

#[test]
fn threshold_generator_agrees_with_the_threshold_envelope() {
    let env = Envelope::from_threshold(&iso("C6H12O6"), 1e-6, Threshold::Absolute).unwrap();
    let peaks: Vec<_> = ThresholdGenerator::new(iso("C6H12O6"), 1e-6, Threshold::Absolute)
        .unwrap()
        .collect();
    assert_eq!(peaks.len(), 21);
    assert_eq!(
        sorted_masses(peaks.iter().map(|p| p.mass).collect()),
        sorted_masses(env.masses().to_vec())
    );
    for p in &peaks {
        assert!(p.prob >= 1e-6);
        assert_close(p.prob.ln(), p.lprob, 1e-9);
    }
    let unordered = ThresholdGenerator::with_options(
        iso("C6H12O6"),
        1e-6,
        Threshold::Absolute,
        ThresholdOptions {
            reorder_marginals: false,
        },
    )
    .unwrap();
    assert_eq!(unordered.count(), 21);
}

#[test]
fn threshold_cursor_count_and_reset() {
    let mut g = ThresholdGenerator::new(iso("C100H202"), 1e-3, Threshold::Relative).unwrap();
    assert_eq!(g.current(), None);
    assert!(g.conf().is_none());
    let n = g.count_confs().unwrap();
    assert!(n > 1);
    let first = g.advance().unwrap().unwrap();
    assert_eq!(g.current(), Some(first));
    let conf = g.conf().unwrap();
    assert_eq!(conf.len(), 2);
    let atoms: Vec<i32> = conf.elements().map(|e| e.iter().sum()).collect();
    assert_eq!(atoms, vec![100, 202]);
    let owned = g.conf_owned().unwrap();
    let mut flat = vec![0; g.all_dim()];
    assert!(g.write_conf(&mut flat));
    assert_eq!(flat, owned.counts());
    assert_eq!(g.isotope_numbers(), &[2, 2]);
    let rest = g.by_ref().count();
    assert_eq!(rest + 1, n);
    // Exhausted, and stays so.
    assert_eq!(g.advance().unwrap(), None);
    assert!(g.conf().is_none());
    assert!(!g.write_conf(&mut flat));
    g.reset().unwrap();
    assert_eq!(g.count(), n);
    // Counting mid-way starts over too.
    let mut g = ThresholdGenerator::new(iso("C100"), 1e-3, Threshold::Relative).unwrap();
    g.next();
    let single = g.count_confs().unwrap();
    assert_eq!(g.count(), single);
}

#[test]
#[should_panic(expected = "wrong length")]
fn write_conf_checks_the_buffer() {
    let mut g = ThresholdGenerator::new(iso("H2O1"), 1e-3, Threshold::Absolute).unwrap();
    g.next();
    g.write_conf(&mut [0; 2]);
}

#[test]
fn layered_generator_covers_like_total_prob() {
    let mut covered = 0.0;
    let mut n = 0;
    for p in LayeredGenerator::new(iso("C100H202")).unwrap() {
        covered += p.prob;
        n += 1;
        if covered >= 0.99 {
            break;
        }
    }
    // The optimal set for 0.99 has 8 peaks; a layer may hold more.
    assert!(n >= 8);
    let options = LayeredOptions {
        reorder_marginals: false,
        prob_hint: 0.5,
    };
    let all: f64 = LayeredGenerator::with_options(iso("H2O1"), options)
        .unwrap()
        .map(|p| p.prob)
        .sum();
    assert_close(all, 1.0, 1e-12);
}

#[test]
fn ordered_generator_is_ordered_and_complete() {
    let peaks: Vec<_> = OrderedGenerator::new(iso("C6H12O6")).unwrap().collect();
    // C6H12O6 has 7 * 13 * 28 configurations.
    assert_eq!(peaks.len(), 7 * 13 * 28);
    assert!(peaks.windows(2).all(|w| w[0].prob >= w[1].prob));
    let total: f64 = peaks.iter().map(|p| p.prob).sum();
    assert_close(total, 1.0, 1e-12);
    assert_close(peaks[0].mass, iso("C6H12O6").mode_mass(), 1e-12);
}

#[test]
fn stochastic_generator_is_reproducible_when_seeded() {
    let options = StochasticOptions {
        seed: Some(42),
        ..Default::default()
    };
    let collect = || {
        let mut v: Vec<(f64, f64)> = StochasticGenerator::with_options(iso("C100"), 1000, options)
            .unwrap()
            .map(|p| (p.mass, p.prob))
            .collect();
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        v
    };
    let a = collect();
    assert_eq!(a, collect());
    let total: f64 = a.iter().map(|p| p.1).sum();
    assert!(total <= 1000.0 && total > 990.0);
    let mut g = StochasticGenerator::new(iso("C100"), 100).unwrap();
    assert_eq!(g.molecule_count(), None);
    g.next();
    assert!(g.molecule_count().unwrap() >= 1);
}

#[test]
fn generators_reject_bad_parameters() {
    assert!(ThresholdGenerator::new(iso("C1"), f64::NAN, Threshold::Absolute).is_err());
    assert!(
        StochasticGenerator::with_options(
            iso("C1"),
            1,
            StochasticOptions {
                precision: 2.0,
                ..Default::default()
            }
        )
        .is_err()
    );
    let empty = Iso::from_isotopes(&[]).unwrap();
    assert!(ThresholdGenerator::new(empty.clone(), 0.1, Threshold::Absolute).is_err());
    assert!(LayeredGenerator::new(empty.clone()).is_err());
    assert!(OrderedGenerator::new(empty).is_err());
}
