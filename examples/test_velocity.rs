//! Quick verification of wavefunction-gradient velocity extraction.
//!
//! Run with:
//!   cargo run --example test_velocity --release

use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric;
use ndarray::IxDyn;

use hermes_rs::core::schrodinger_dynamics::extract_velocity;
use hermes_rs::physics::cosmology::planck_2018;

fn main() {
    let n = 64;
    let box_length = 10000.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let dx = box_length / n as f64;

    let cosmology = planck_2018();
    let rho_bar = cosmology.density_matter();
    let ell_over_m = 7500.0;
    let mass = 1e10;
    let ell = ell_over_m * mass;

    let k0 = 2.0 * PI * 3.0 / box_length;
    let v0 = 100.0; // kpc/Gyr
    let phase_amp = mass * v0 / (ell * k0); // ~7 radians

    println!("Phase amplitude: {phase_amp:.2} radians (wraps through 2pi)");
    println!();

    // alpha(x) = sqrt(rho_bar/m) * exp(-I * phase_amp * cos(k0 x))
    // Analytic velocity: v_x = v0 * sin(k0 x)
    let psi = EvenField::from_fn(&grid, g, |x| {
        let amp = (rho_bar / mass).sqrt();
        let phase = -phase_amp * (k0 * x[0]).cos();
        (amp * phase.cos(), amp * phase.sin())
    });

    let velocity = extract_velocity(&psi, &grid, ell, mass);

    let mut max_error = 0.0_f64;
    println!(
        "{:>5}  {:>12}  {:>12}  {:>12}",
        "index", "v_analytic", "v_extracted", "error"
    );

    for m0 in 0..n {
        let x = m0 as f64 * dx; // morphis uses cell corners
        let v_analytic = v0 * (k0 * x).sin();
        let v_extracted = velocity[0][IxDyn(&[m0, 0, 0])];
        let error = (v_extracted - v_analytic).abs();
        if error > max_error {
            max_error = error;
        }
        if m0 % 8 == 0 {
            println!("{m0:5}  {v_analytic:12.4e}  {v_extracted:12.4e}  {error:12.4e}");
        }
    }

    println!();
    println!("Max absolute error: {max_error:.4e}");
    println!("Relative to v0={v0}: {:.4e}", max_error / v0);
}
