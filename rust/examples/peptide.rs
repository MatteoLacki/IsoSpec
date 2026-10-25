//! A modified peptide: its composition, and its spectrum binned to 1 Da.
//!
//!     cargo run --release --example peptide -- 'PEPTC[UNIMOD:4]DEK'

use isospec::{Composition, Envelope, Iso, PeptideOptions};

fn main() -> Result<(), isospec::Error> {
    let sequence = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "PEPTC[UNIMOD:4]DEK".into());
    let options = PeptideOptions::default(); // with the terminal water, embedded Unimod table

    let composition = Composition::from_peptide(&sequence, &options)?;
    println!("{sequence}: {composition}");

    let iso = Iso::from_peptide(&sequence, &options)?;
    println!("monoisotopic mass {:.6} Da", iso.monoisotopic_mass());

    // Binned directly from the molecule: no fine structure kept.
    let mut binned = Envelope::binned(&iso, 0.9999, 1.0, 0.0)?;
    binned.sort_by_mass();
    for (mass, prob) in binned.iter() {
        println!("{mass:.1}\t{prob:.6}");
    }
    Ok(())
}
