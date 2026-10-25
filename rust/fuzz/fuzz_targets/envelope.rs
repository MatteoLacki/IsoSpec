// Envelopes from arbitrary arrays through a sequence of operations, as
// tests/fuzz/fuzz_envelope.cpp does for the C++ API, here through the safe
// Rust one. Every operation either succeeds or returns an Err, and leaves the
// envelope's arrays finite.
#![no_main]

use arbitrary::Arbitrary;
use isospec::Envelope;
use libfuzzer_sys::fuzz_target;

#[derive(Arbitrary, Debug)]
enum Op {
    Scale(f64),
    Shift(f64),
    Normalize,
    Resample(u16, f64, u64),
    Bin(f64, f64),
    Add,
    Convolve,
    Combine(f64, f64),
    Wasserstein,
    Oriented,
    Abyssal(f64, f64),
    Match(f64, f64),
    SortMass,
    SortProb,
    Stats,
    Clone,
    Swap,
}

#[derive(Arbitrary, Debug)]
struct Input {
    a: Vec<(f64, f64)>,
    b: Vec<(f64, f64)>,
    ops: Vec<Op>,
}

fn build(peaks: &[(f64, f64)]) -> Option<Envelope> {
    let masses: Vec<f64> = peaks.iter().map(|p| p.0).collect();
    let probs: Vec<f64> = peaks.iter().map(|p| p.1).collect();
    Envelope::from_arrays(&masses, &probs).ok()
}

fn check(env: &Envelope) {
    assert_eq!(env.masses().len(), env.len());
    assert_eq!(env.probs().len(), env.len());
    assert!(env.masses().iter().all(|m| m.is_finite()));
    assert!(env.probs().iter().all(|p| p.is_finite() && *p >= 0.0));
}

fuzz_target!(|input: Input| {
    if input.a.len() > 64 || input.b.len() > 64 || input.ops.len() > 16 {
        return;
    }
    let (Some(mut a), Some(mut b)) = (build(&input.a), build(&input.b)) else {
        return;
    };
    for op in &input.ops {
        let r: isospec::Result<()> = match *op {
            Op::Scale(f) => a.scale(f),
            Op::Shift(s) => a.shift_mass(s),
            Op::Normalize => a.normalize(),
            Op::Resample(n, beta, seed) => a.resample(n as usize, beta, Some(seed)),
            Op::Bin(w, m) => a.bin(w, m).map(|e| {
                if e.len() <= 4096 {
                    a = e;
                }
            }),
            Op::Add => a.try_add(&b).map(|e| a = e),
            Op::Convolve => {
                if a.len() * b.len() <= 4096 {
                    a.convolve(&b).map(|e| a = e)
                } else {
                    Ok(())
                }
            }
            Op::Combine(x, y) => Envelope::linear_combination(&[(&a, x), (&b, y)]).map(|e| a = e),
            Op::Wasserstein => a.wasserstein(&mut b).map(drop),
            Op::Oriented => a.oriented_wasserstein(&mut b).map(drop),
            Op::Abyssal(d, s) => a.abyssal_wasserstein(&mut b, d, s).map(drop),
            Op::Match(d, s) => a.wasserstein_match(&mut b, d, s).map(|m| {
                assert!(m.unmatched_this >= 0.0 && m.unmatched_other >= 0.0 && m.flow >= 0.0);
            }),
            Op::SortMass => {
                a.sort_by_mass();
                assert!(a.masses().windows(2).all(|w| w[0] <= w[1]));
                Ok(())
            }
            Op::SortProb => {
                a.sort_by_prob();
                Ok(())
            }
            Op::Stats => {
                let _ = a.total_prob();
                a.average_mass().and_then(|_| a.variance()).map(drop)
            }
            Op::Clone => {
                a = a.clone();
                Ok(())
            }
            Op::Swap => {
                std::mem::swap(&mut a, &mut b);
                Ok(())
            }
        };
        let _ = r;
        check(&a);
        check(&b);
    }
});
