// Envelopes: every constructor, the accessors, and the operations. Expected
// values come from IsoSpecPy 2.5.0.

mod common;
use common::{assert_close, close};

use isospec::{Envelope, ErrorKind, Iso, StochasticOptions, Threshold};

fn iso(formula: &str) -> Iso {
    Iso::from_formula(formula).unwrap()
}

#[test]
fn total_prob_matches_python() {
    let mut env = Envelope::from_total_prob(&iso("C100H202"), 0.99).unwrap();
    assert_eq!(env.len(), 8);
    assert_close(env.total_prob(), 0.9929694469968894, 1e-12);
    let mut layered = iso("C100H202").envelope().total_prob_layered(0.99).unwrap();
    assert_eq!(layered.len(), 8);
    assert_close(layered.total_prob(), 0.9922602266897552, 1e-12);
}

#[test]
fn threshold_matches_python() {
    let mut env = Envelope::from_threshold(&iso("C6H12O6"), 1e-6, Threshold::Absolute).unwrap();
    assert_eq!(env.len(), 21);
    assert_close(env.total_prob(), 0.9999976192461025, 1e-12);
    let rel = Envelope::from_threshold(&iso("C6H12O6"), 1e-2, Threshold::Relative).unwrap();
    assert_eq!(rel.len(), 3);
}

#[test]
fn seeded_stochastic_envelope_is_the_same_as_python_s() {
    // IsoSpecPy.IsoStochastic(1000, formula="C100", seed=42)
    let expected = [
        (1200.0, 329.0),
        (1201.0033548352, 362.0),
        (1202.0067096704, 202.0),
        (1203.0100645056, 84.0),
        (1204.0134193408, 17.0),
        (1205.016774176, 3.0),
        (1206.0201290112, 3.0),
    ];
    let options = StochasticOptions {
        seed: Some(42),
        ..Default::default()
    };
    let mut env = iso("C100").envelope().stochastic(1000, options).unwrap();
    env.sort_by_mass();
    let got: Vec<(f64, f64)> = env.iter().collect();
    assert_eq!(got.len(), expected.len());
    for ((m, p), (em, ep)) in got.iter().zip(expected) {
        assert_close(*m, em, 1e-12);
        assert_eq!(*p, ep);
    }
    // And the same twice.
    let again = iso("C100").envelope().stochastic(1000, options).unwrap();
    let mut again = again;
    again.sort_by_mass();
    assert_eq!(again.probs(), env.probs());
}

#[test]
fn seeded_resample_is_the_same_as_python_s() {
    // x = IsoTotalProb(0.999, formula="C6H12O6"); x.resample(1000, seed=7)
    let expected = [
        (180.06338810844, 925.0),
        (181.06674294364, 59.0),
        (181.06760524584, 1.0),
        (181.06966485435998, 1.0),
        (182.06763310194, 8.0),
        (182.07009777884, 2.0),
        (183.07098793714, 4.0),
    ];
    let mut env = Envelope::from_total_prob(&iso("C6H12O6"), 0.999).unwrap();
    env.resample(1000, 1.0, Some(7)).unwrap();
    env.sort_by_mass();
    let got: Vec<(f64, f64)> = env.iter().filter(|&(_, p)| p > 0.0).collect();
    assert_eq!(got.len(), expected.len());
    for ((m, p), (em, ep)) in got.iter().zip(expected) {
        assert_close(*m, em, 1e-12);
        assert_eq!(*p, ep);
    }
}

#[test]
fn binned_matches_python() {
    let mut env = Envelope::binned(&iso("C100H202"), 0.999, 1.0, 0.0).unwrap();
    env.sort_by_mass();
    assert_eq!(env.len(), 7);
    assert_eq!(&env.masses()[..3], &[1404.0, 1405.0, 1406.0]);
    assert_close(env.probs()[0], 0.33020491691200565, 1e-12);
    assert_close(env.probs()[1], 0.367830792929928, 1e-12);
    assert_close(env.total_prob(), 0.9990856895148157, 1e-12);
    assert_eq!(
        Envelope::binned(&iso("C1"), 0.9, 0.0, 0.0)
            .unwrap_err()
            .kind(),
        ErrorKind::InvalidArgument
    );
}

#[test]
fn wasserstein_distances_match_python() {
    let mut a = Envelope::from_total_prob(&iso("C6H12O6"), 0.999).unwrap();
    let mut b = Envelope::from_total_prob(&iso("C6H12O5N1"), 0.999).unwrap();
    a.normalize().unwrap();
    b.normalize().unwrap();
    assert_close(a.wasserstein(&mut b).unwrap(), 1.9930987487082414, 1e-10);
    assert_close(
        a.oriented_wasserstein(&mut b).unwrap(),
        -1.9930987487082414,
        1e-10,
    );
    let m = a.wasserstein_match(&mut b, 0.5, 1.0).unwrap();
    assert_close(m.unmatched_this, 0.988260559968297, 1e-10);
    assert_close(m.unmatched_other, 0.9882605599682968, 1e-10);
    assert_close(m.flow, 0.011739440031703186, 1e-10);
    // With an abyss deeper than any transport, nothing is condemned.
    let d = a.abyssal_wasserstein(&mut b, 100.0, 1.0).unwrap();
    assert_close(d, 1.9930987487082414, 1e-10);
    // Unnormalised envelopes are rejected by the plain distance.
    a.scale(2.0).unwrap();
    assert!(a.wasserstein(&mut b).is_err());
}

#[test]
fn configurations_ride_along() {
    let iso = iso("C6H12O6");
    let mut env = iso.envelope().with_confs().total_prob(0.999).unwrap();
    assert!(env.has_confs());
    env.sort_by_prob();
    let n = env.len();
    // The most probable is the monoisotopic one.
    assert_eq!(
        env.conf(n - 1).unwrap().counts(),
        iso.monoisotopic_conf().counts()
    );
    let confs: Vec<_> = env.confs().unwrap().collect();
    assert_eq!(confs.len(), n);
    for (ii, conf) in confs.iter().enumerate() {
        assert_eq!(conf.len(), 3);
        // Every configuration has the molecule's atoms.
        let atoms: Vec<i32> = conf.elements().map(|e| e.iter().sum()).collect();
        assert_eq!(atoms, vec![6, 12, 6]);
        assert_eq!(conf.to_owned().counts(), env.conf(ii).unwrap().counts());
    }
    // Operations keep them where the library does, and drop them where it
    // does not; a clone keeps them.
    assert!(env.clone().has_confs());
    assert!(!env.convolve(&env).unwrap().has_confs());
    let plain = Envelope::from_total_prob(&iso, 0.9).unwrap();
    assert!(!plain.has_confs());
    assert!(plain.conf(0).is_none());
    assert!(plain.confs().is_none());
}

#[test]
#[should_panic(expected = "peak index out of range")]
fn conf_index_is_checked() {
    let env = iso("H2O1").envelope().with_confs().total_prob(0.9).unwrap();
    env.conf(env.len());
}

#[test]
fn arrays_in_and_out() {
    let mut env = Envelope::from_arrays(&[3.0, 1.0, 2.0], &[0.2, 0.5, 0.3]).unwrap();
    assert_eq!(env.len(), 3);
    assert!(!env.is_empty());
    env.sort_by_mass();
    assert_eq!(env.masses(), &[1.0, 2.0, 3.0]);
    assert_eq!(env.probs(), &[0.5, 0.3, 0.2]);
    env.sort_by_prob();
    assert_eq!(env.probs(), &[0.2, 0.3, 0.5]);
    let (m, p) = env.to_vecs();
    assert_eq!(m, vec![3.0, 2.0, 1.0]);
    assert_eq!(p, vec![0.2, 0.3, 0.5]);
    assert_close(
        env.average_mass().unwrap(),
        0.2 * 3.0 + 0.3 * 2.0 + 0.5,
        1e-12,
    );
    let var = env.variance().unwrap();
    assert_close(env.stddev().unwrap(), var.sqrt(), 1e-12);
    env.shift_mass(10.0).unwrap();
    assert!(env.masses().iter().all(|&m| m > 10.0));
    env.scale(2.0).unwrap();
    assert_close(env.total_prob(), 2.0, 1e-12);
    env.normalize().unwrap();
    assert_close(env.total_prob(), 1.0, 1e-12);

    let empty = Envelope::from_arrays(&[], &[]).unwrap();
    assert!(empty.is_empty());
    assert!(empty.masses().is_empty());
    assert_eq!(Envelope::default().len(), 0);

    let kind = |r: isospec::Result<Envelope>| r.unwrap_err().kind();
    assert_eq!(
        kind(Envelope::from_arrays(&[1.0], &[])),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Envelope::from_arrays(&[f64::NAN], &[1.0])),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Envelope::from_arrays(&[1.0], &[-1.0])),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        kind(Envelope::from_arrays(&[1.0], &[f64::INFINITY])),
        ErrorKind::InvalidArgument
    );
}

#[test]
fn sums_convolutions_and_combinations() {
    let a = Envelope::from_arrays(&[1.0, 2.0], &[0.5, 0.5]).unwrap();
    let b = Envelope::from_arrays(&[10.0], &[1.0]).unwrap();
    let mut sum = &a + &b;
    assert_eq!(sum.len(), 3);
    assert_close(sum.total_prob(), 2.0, 1e-12);
    let mut conv = &a * &b;
    conv.sort_by_mass();
    assert_eq!(conv.masses(), &[11.0, 12.0]);
    assert_eq!(conv.probs(), &[0.5, 0.5]);
    let mut lc = Envelope::linear_combination(&[(&a, 2.0), (&b, 3.0)]).unwrap();
    assert_close(lc.total_prob(), 5.0, 1e-12);
    assert!(Envelope::linear_combination(&[]).unwrap().is_empty());
    // A convolution whose masses overflow is an error, not inf.
    let mut far = Envelope::from_arrays(&[1e100], &[1.0]).unwrap();
    far.shift_mass(1.7e308).unwrap();
    assert_eq!(far.convolve(&far).unwrap_err().kind(), ErrorKind::Overflow);
    // Convolving two glucoses gives the distribution of two glucoses.
    let g = Envelope::from_threshold(&iso("C6H12O6"), 0.0, Threshold::Absolute).unwrap();
    let mut gg = g.convolve(&g).unwrap();
    let two = iso("C12H24O12");
    assert_close(gg.average_mass().unwrap(), two.average_mass(), 1e-9);
}

#[test]
fn binning_an_envelope() {
    let mut env = Envelope::from_arrays(&[1.1, 1.2, 2.4, 3.9], &[0.1, 0.2, 0.3, 0.4]).unwrap();
    let mut binned = env.bin(1.0, 0.0).unwrap();
    binned.sort_by_mass();
    assert_eq!(binned.masses(), &[1.0, 2.0, 4.0]);
    assert!(close(binned.probs()[0], 0.3, 1e-12));
    assert!(env.bin(f64::NAN, 0.0).is_err());
}

#[test]
fn scalar_arguments_are_checked() {
    let mut env = Envelope::from_arrays(&[1.0], &[1.0]).unwrap();
    assert_eq!(
        env.scale(f64::NAN).unwrap_err().kind(),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        env.scale(-1.0).unwrap_err().kind(),
        ErrorKind::InvalidArgument
    );
    assert_eq!(
        env.scale(f64::MAX)
            .and_then(|_| env.scale(f64::MAX))
            .unwrap_err()
            .kind(),
        ErrorKind::Overflow
    );
    let iso = iso("C6H12O6");
    assert!(Envelope::from_threshold(&iso, f64::NAN, Threshold::Absolute).is_err());
    assert!(Envelope::from_total_prob(&iso, f64::NAN).is_err());
    assert!(
        iso.envelope()
            .stochastic(
                10,
                StochasticOptions {
                    precision: 0.0,
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        iso.envelope()
            .stochastic(
                10,
                StochasticOptions {
                    beta_bias: -1.0,
                    ..Default::default()
                }
            )
            .is_err()
    );
}
