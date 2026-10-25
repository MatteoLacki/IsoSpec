// Shared by the integration tests (each includes it with `mod common;`).
#![allow(dead_code)]

pub fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-300)
}

#[track_caller]
pub fn assert_close(a: f64, b: f64, rel: f64) {
    assert!(close(a, b, rel), "{a} != {b} (relative tolerance {rel})");
}
