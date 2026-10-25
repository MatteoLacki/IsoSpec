//! Streaming isotopologues instead of storing them: the most probable ones in
//! order, and everything above a threshold with configurations.
//!
//!     cargo run --release --example generators

use isospec::{Iso, OrderedGenerator, Threshold, ThresholdGenerator};

fn main() -> Result<(), isospec::Error> {
    let iso = Iso::from_formula("C100H202S2")?;

    println!("The five most probable isotopologues:");
    for peak in OrderedGenerator::new(iso.clone())?.into_iter().take(5) {
        println!("  {:.6} Da  p = {:.6}", peak.mass, peak.prob);
    }

    // The cursor gives each isotopologue's configuration as well.
    println!("Isotopologues within 1% of the most probable, with their isotopes:");
    let mut generator = ThresholdGenerator::new(iso, 0.01, Threshold::Relative)?;
    while let Some(peak) = generator.advance()? {
        let conf = generator.conf().expect("positioned on an isotopologue");
        println!("  {:.6} Da  p = {:.6}  {:?}", peak.mass, peak.prob, conf);
    }
    Ok(())
}
