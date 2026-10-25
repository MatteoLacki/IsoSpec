// The runtime ISA selection. Its own test binary: the level is process-wide,
// and the other tests run in parallel threads.

use isospec::isa::{self, IsaLevel};
use isospec::{Envelope, Iso};

#[test]
fn levels_switch_and_agree() {
    let start = isa::current();
    assert!(start <= isa::max_supported());
    assert!(IsaLevel::Baseline.is_available());
    assert_eq!(
        isa::current(),
        isa::max_supported(),
        "the default is the best available"
    );
    assert!(!isa::simd_level().is_empty());
    assert_eq!(isospec::simd_level(), isa::simd_level());
    let names: Vec<_> = IsaLevel::ALL.iter().map(|l| l.name()).collect();
    assert_eq!(names, vec!["baseline", "avx", "v3", "v4"]);

    // Every available level computes the same envelope.
    let iso = Iso::from_formula("C520H817N139O147S8").unwrap();
    let mut reference: Option<Vec<(f64, f64)>> = None;
    for level in IsaLevel::ALL {
        if !level.is_available() {
            assert!(!isa::set(level));
            assert_eq!(isa::current(), start, "a refused level changes nothing");
            continue;
        }
        assert!(isa::set(level));
        assert_eq!(isa::current(), level);
        let mut env = Envelope::from_total_prob(&iso, 0.999).unwrap();
        env.sort_by_mass();
        let peaks: Vec<_> = env.iter().collect();
        match &reference {
            None => reference = Some(peaks),
            Some(r) => {
                assert_eq!(r.len(), peaks.len(), "{level:?}");
                for (a, b) in r.iter().zip(&peaks) {
                    assert!((a.0 - b.0).abs() <= 1e-9 * a.0, "{level:?}");
                    assert!((a.1 - b.1).abs() <= 1e-12, "{level:?}");
                }
            }
        }
    }
    assert!(isa::set(start));
    assert_eq!(isa::current(), start);
}
