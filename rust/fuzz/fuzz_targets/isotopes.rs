// Molecules from raw isotope data -- arbitrary doubles and counts -- through
// every Iso accessor, the generators (cursor, configurations and batched
// iteration) and the envelopes built from them.
#![no_main]

use arbitrary::Arbitrary;
use isospec::{
    ElementSpec, Envelope, Iso, LayeredGenerator, OrderedGenerator, StochasticGenerator,
    StochasticOptions, Threshold, ThresholdGenerator,
};
use libfuzzer_sys::fuzz_target;

#[derive(Arbitrary, Debug)]
struct Element {
    atoms: u8,
    isotopes: Vec<(f64, f64)>,
}

#[derive(Arbitrary, Debug)]
struct Input {
    elements: Vec<Element>,
    extra: Option<Element>,
    threshold: f64,
    which: u8,
    steps: u16,
    seed: u64,
}

fuzz_target!(|input: Input| {
    if input.elements.len() > 4 || input.elements.iter().any(|e| e.isotopes.len() > 5) {
        return;
    }
    let masses: Vec<Vec<f64>> = input
        .elements
        .iter()
        .map(|e| e.isotopes.iter().map(|i| i.0).collect())
        .collect();
    let probs: Vec<Vec<f64>> = input
        .elements
        .iter()
        .map(|e| e.isotopes.iter().map(|i| i.1).collect())
        .collect();
    let specs: Vec<ElementSpec> = input
        .elements
        .iter()
        .zip(masses.iter().zip(&probs))
        .map(|(e, (m, p))| ElementSpec::new(e.atoms as u32 % 64, m, p))
        .collect();
    let Ok(mut iso) = Iso::from_isotopes(&specs) else {
        return;
    };
    if let Some(extra) = &input.extra {
        let m: Vec<f64> = extra.isotopes.iter().map(|i| i.0).collect();
        let p: Vec<f64> = extra.isotopes.iter().map(|i| i.1).collect();
        let dims = iso.dim_number();
        if iso
            .add_element(&ElementSpec::new(extra.atoms as u32 % 64, &m, &p))
            .is_err()
        {
            assert_eq!(iso.dim_number(), dims);
        }
    }
    let _ = (
        iso.monoisotopic_mass(),
        iso.lightest_peak_mass(),
        iso.heaviest_peak_mass(),
        iso.stddev(),
    );
    for c in [
        iso.lightest_peak_conf(),
        iso.heaviest_peak_conf(),
        iso.monoisotopic_conf(),
    ] {
        assert_eq!(c.counts().len(), iso.all_dim());
        for (ii, part) in c.elements().enumerate() {
            assert_eq!(part.iter().sum::<i32>(), iso.atom_counts()[ii]);
        }
    }
    for ii in 0..iso.dim_number() {
        assert_eq!(
            iso.element_masses(ii).len(),
            iso.isotope_numbers()[ii] as usize
        );
        assert_eq!(
            iso.element_lprobs(ii).len(),
            iso.isotope_numbers()[ii] as usize
        );
    }
    let steps = input.steps as usize % 2000;
    let threshold = if input.threshold.is_finite() {
        input.threshold.abs().max(1e-12)
    } else {
        input.threshold
    };
    match input.which % 6 {
        0 => {
            if let Ok(mut g) = ThresholdGenerator::new(iso.clone(), threshold, Threshold::Relative)
            {
                for _ in 0..steps {
                    if g.advance().unwrap().is_none() {
                        break;
                    }
                    let conf = g.conf().unwrap();
                    assert_eq!(conf.counts().len(), iso.all_dim());
                }
                let _ = g.count_confs();
                let n = g.into_iter().take(steps).count();
                assert!(n <= steps);
            }
        }
        1 => {
            if let Ok(g) = LayeredGenerator::new(iso.clone()) {
                let total: f64 = g.into_iter().take(steps).map(|p| p.prob).sum();
                assert!(total <= 1.0 + 1e-9);
            }
        }
        2 => {
            if let Ok(mut g) = OrderedGenerator::new(iso.clone()) {
                let mut last = f64::INFINITY;
                for _ in 0..steps {
                    match g.advance().unwrap() {
                        Some(p) => {
                            assert!(p.prob <= last * (1.0 + 1e-9));
                            last = p.prob;
                            let _ = g.conf_owned();
                        }
                        None => break,
                    }
                }
            }
        }
        3 => {
            let options = StochasticOptions::default().with_seed(input.seed);
            if let Ok(g) = StochasticGenerator::with_options(iso.clone(), steps, options) {
                let total: f64 = g.into_iter().map(|p| p.prob).sum();
                assert!(total <= steps as f64);
            }
        }
        4 => {
            if let Ok(env) = iso
                .envelope()
                .with_confs()
                .threshold(threshold, Threshold::Relative)
            {
                if let Some(confs) = env.confs() {
                    for c in confs.take(steps) {
                        assert_eq!(c.counts().len(), iso.all_dim());
                    }
                }
            }
        }
        _ => {
            let width = if threshold.is_finite() {
                threshold
            } else {
                1.0
            };
            if let Ok(env) = Envelope::binned(&iso, 0.9, width, input.seed as f64) {
                assert!(env.probs().iter().all(|p| p.is_finite() && *p > 0.0));
            }
        }
    }
});
