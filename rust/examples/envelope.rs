//! The isotopologues of a molecule covering 99.9% of its probability, by mass.
//!
//!     cargo run --release --example envelope -- C254H377N65O75S6 0.999

use isospec::Iso;

fn main() -> Result<(), isospec::Error> {
    let mut args = std::env::args().skip(1);
    let formula = args.next().unwrap_or_else(|| "C254H377N65O75S6".into());
    let coverage: f64 = args
        .next()
        .map_or(0.999, |p| p.parse().expect("a probability"));

    let iso = Iso::from_formula(&formula)?;
    println!(
        "{formula}: monoisotopic {:.6} Da, average {:.6} Da",
        iso.monoisotopic_mass(),
        iso.average_mass()
    );

    let mut env = iso.envelope().total_prob(coverage)?;
    env.sort_by_mass();
    println!("{} isotopologues cover {:.6}", env.len(), env.total_prob());
    for (mass, prob) in env.iter() {
        println!("{mass:.6}\t{prob:.6e}");
    }
    Ok(())
}
