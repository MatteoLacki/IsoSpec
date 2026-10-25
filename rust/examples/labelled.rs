//! A molecule with a labelled element: glucose with five of its six carbons
//! 99% carbon-13, compared with natural glucose by Wasserstein distance.
//!
//!     cargo run --release --example labelled

use isospec::{Element, ElementSpec, Envelope, Iso};

fn main() -> Result<(), isospec::Error> {
    let carbon = Element::from_symbol("C").expect("carbon");
    let masses: Vec<f64> = carbon.isotopes().map(|i| i.mass).collect();

    // One natural carbon and the rest of glucose from the formula ...
    let mut labelled = Iso::from_formula("C1H12O6")?;
    // ... and five labelled carbons as an element of their own.
    labelled.add_element(&ElementSpec {
        atoms: 5,
        masses: &masses,
        probabilities: &[0.01, 0.99],
    })?;

    let natural = Iso::from_formula("C6H12O6")?;
    let mut a = Envelope::from_total_prob(&natural, 0.9999)?;
    let mut b = Envelope::from_total_prob(&labelled, 0.9999)?;
    a.normalize()?;
    b.normalize()?;
    println!("natural  glucose: {:.4} Da on average", a.average_mass()?);
    println!("labelled glucose: {:.4} Da on average", b.average_mass()?);
    println!("Wasserstein distance: {:.4} Da", a.wasserstein(&mut b)?);
    Ok(())
}
