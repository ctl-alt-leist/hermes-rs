//! Constrained Zel'dovich initialization for galaxy group formation.
//!
//! Uses the same Zel'dovich approximation as the cosmic web, but adds
//! a uniform overdensity bias to the initial density field. This
//! ensures the entire box is a Lagrangian region that will collapse
//! into a group-mass structure by z = 0.
//!
//! The bias δ_0 is chosen so that the linear overdensity at z = 0
//! exceeds the spherical collapse threshold δ_c ≈ 1.686. For a
//! starting redshift z_init, the required initial overdensity is:
//!
//! ```text
//! δ_init = δ_c / D_+(a = 1) × D_+(a_init)
//! ```
//!
//! In practice we use a slightly lower value so the collapse happens
//! around z ~ 0.5-1, giving the group time to virialize.

use crate::scenes::cosmic_web::init::power_spectrum;

use ndarray::Array3;
use ndrustfft::{FftHandler, R2cFftHandler, ndifft, ndifft_r2c};
use num_complex::Complex64;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::algebra::vector_from_components;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// Generate constrained Zel'dovich initial conditions for a galaxy group.
///
/// Same pipeline as the standard Zel'dovich approximation, but with a
/// uniform overdensity bias added to the density field. The bias
/// ensures the box collapses into a group-mass halo.
pub fn constrained_zeldovich_init(
    n_per_side: usize,
    grid: &Grid,
    cosmology: &Cosmology,
    scale_factor_initial: f64,
    seed: u64,
) -> Result<Particles, HermesError> {
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;
    let box_length = grid.box_length;

    // Overdensity bias: target δ_c / D+(1) at the initial redshift.
    // δ_c ≈ 1.686 for spherical collapse. We use 0.8 × δ_c so the
    // collapse completes around z ~ 0.5, giving time to virialize.
    let collapse_fraction = 0.8;
    let delta_c = 1.686;
    let growth_initial = cosmology.growth_factor(scale_factor_initial);
    let overdensity_bias = collapse_fraction * delta_c * growth_initial;

    // Generate Gaussian random overdensity in Fourier space.
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let mut delta_hat = Array3::<Complex64>::zeros((n, n, n_complex));

    let k_to_hmpc = 1000.0 / cosmology.hubble;
    let volume_box = box_length * box_length * box_length;

    for m0 in 0..n {
        let kx = grid.wavevector_component(m0);
        for m1 in 0..n {
            let ky = grid.wavevector_component(m1);
            for m2 in 0..n_complex {
                let kz = grid.wavevector_component(m2);
                let k2 = kx * kx + ky * ky + kz * kz;

                if k2 < 1e-30 {
                    // Zero mode: set to the overdensity bias.
                    // In Fourier convention, δ̂(k=0) = δ_0 × V_box^{1/2}
                    // (the mean overdensity times the box volume normalization).
                    delta_hat[[m0, m1, m2]] = Complex64::new(
                        overdensity_bias * volume_box.sqrt() / k_to_hmpc.powf(1.5),
                        0.0,
                    );
                    continue;
                }

                let k_mag = k2.sqrt();
                let k_hmpc = k_mag * k_to_hmpc;
                let power = power_spectrum(k_hmpc, cosmology);

                let sigma = (power * volume_box / k_to_hmpc.powi(3)).sqrt();

                let re: f64 = rng.random_range(-1.0..1.0);
                let im: f64 = rng.random_range(-1.0..1.0);
                delta_hat[[m0, m1, m2]] = Complex64::new(re * sigma, im * sigma);
            }
        }
    }

    // Compute displacement field: Ψ̂_d(k) = i k_d / k² × δ̂(k)
    let growth = cosmology.growth_factor(scale_factor_initial);

    let mut displacement: [Array3<f64>; 3] = std::array::from_fn(|_| Array3::zeros((n, n, n)));

    for (d, displacement_d) in displacement.iter_mut().enumerate() {
        let mut psi_hat = Array3::<Complex64>::zeros((n, n, n_complex));

        for m0 in 0..n {
            let kx = grid.wavevector_component(m0);
            for m1 in 0..n {
                let ky = grid.wavevector_component(m1);
                for m2 in 0..n_complex {
                    let kz = grid.wavevector_component(m2);
                    let k2 = kx * kx + ky * ky + kz * kz;

                    if k2 < 1e-30 {
                        continue;
                    }

                    let kd = match d {
                        0 => kx,
                        1 => ky,
                        _ => kz,
                    };

                    psi_hat[[m0, m1, m2]] = Complex64::new(0.0, kd / k2) * delta_hat[[m0, m1, m2]];
                }
            }
        }

        let handler_c2c_0 = FftHandler::new(n);
        let handler_c2c_1 = FftHandler::new(n);
        let handler_r2c = R2cFftHandler::new(n);
        let mut scratch = psi_hat.clone();

        ndifft(&psi_hat, &mut scratch, &handler_c2c_0, 0);
        psi_hat.assign(&scratch);
        ndifft(&psi_hat, &mut scratch, &handler_c2c_1, 1);
        psi_hat.assign(&scratch);

        let mut real_out = Array3::<f64>::zeros((n, n, n));
        ndifft_r2c(&psi_hat, &mut real_out, &handler_r2c, 2);

        *displacement_d = real_out;
    }

    // Place particles on Lagrangian lattice and displace.
    let density_mean = cosmology.density_matter();
    let mut particles = Particles::on_lattice(n_per_side, grid, density_mean);

    let growth_rate = cosmology.growth_rate(scale_factor_initial);
    let hubble_a = cosmology.hubble_parameter(scale_factor_initial);
    let velocity_factor = scale_factor_initial
        * scale_factor_initial
        * growth_rate
        * hubble_a
        * particles.mass_particle;

    for p in 0..particles.count() {
        let pos = particles.position_components(p);
        let psi = interpolate_displacement(&displacement, &pos, grid);

        let position_displaced = vector_from_components(
            pos[0] + growth * psi[0],
            pos[1] + growth * psi[1],
            pos[2] + growth * psi[2],
        );
        particles.set_position(p, &position_displaced);

        let momentum = vector_from_components(
            velocity_factor * psi[0],
            velocity_factor * psi[1],
            velocity_factor * psi[2],
        );
        particles.set_momentum(p, &momentum);
    }

    particles.wrap_positions(grid);

    Ok(particles)
}

/// Interpolate the displacement field at a position using nearest grid point.
fn interpolate_displacement(
    displacement: &[Array3<f64>; 3],
    pos: &[f64; 3],
    grid: &Grid,
) -> [f64; 3] {
    let h_inv = 1.0 / grid.cell_length;
    let n = grid.n_cells;

    let m0 = ((pos[0] * h_inv) as usize).min(n - 1);
    let m1 = ((pos[1] * h_inv) as usize).min(n - 1);
    let m2 = ((pos[2] * h_inv) as usize).min(n - 1);

    [
        displacement[0][[m0, m1, m2]],
        displacement[1][[m0, m1, m2]],
        displacement[2][[m0, m1, m2]],
    ]
}
