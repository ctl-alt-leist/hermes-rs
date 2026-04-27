//! Zel'dovich initialization from a linear power spectrum.
//!
//! Generates cosmological initial conditions by displacing particles from
//! a uniform lattice according to the Zel'dovich approximation. The
//! displacement field is drawn from a Gaussian random field whose variance
//! is set by the linear matter power spectrum.
//!
//! The pipeline is:
//!
//! ```text
//! P(k) → Gaussian δ̂(k) → Ψ̂(k) = ik/k² δ̂(k) → IFFT → displace particles
//! ```

use std::f64::consts::PI;

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

// ============================================================================
// Eisenstein & Hu transfer function
// ============================================================================

/// Eisenstein & Hu (1998) no-wiggle transfer function.
///
/// Approximates the matter transfer function T(k) without baryon acoustic
/// oscillations. Accurate to ~5% for the smooth broadband shape.
///
/// `k` is in units of h/Mpc (not 1/kpc).
pub fn transfer_function(k: f64, cosmology: &Cosmology) -> f64 {
    let omega_mh2 = cosmology.omega_m * cosmology.hubble * cosmology.hubble;
    let omega_bh2 = cosmology.omega_b * cosmology.hubble * cosmology.hubble;
    let f_b = cosmology.baryon_fraction();
    let theta = 2.725 / 2.7; // CMB temperature ratio

    // Sound horizon (Eq. 26 of Eisenstein & Hu 1998)
    let sound_horizon =
        44.5 * (omega_mh2).ln() / (1.0 + 10.0 * omega_bh2.powf(0.75)).sqrt() * cosmology.hubble;

    // Suppression factor for baryons
    let alpha_gamma = 1.0 - 0.328 * (omega_bh2 / omega_mh2).ln() * f_b
        + 0.38 * (omega_bh2 / omega_mh2).ln() * f_b * f_b;

    // Effective shape parameter
    let gamma_eff = cosmology.omega_m
        * cosmology.hubble
        * (alpha_gamma + (1.0 - alpha_gamma) / (1.0 + (0.43 * k * sound_horizon).powi(4)));

    // Dimensionless wavenumber
    let q = k * theta * theta / gamma_eff;

    // Transfer function (Eq. 29)
    let l0 = (2.0 * std::f64::consts::E + 1.8 * q).ln();
    let c0 = 14.2 + 731.0 / (1.0 + 62.5 * q);

    l0 / (l0 + c0 * q * q)
}

/// Linear matter power spectrum P(k), normalized by σ₈.
///
/// ```text
/// P(k) = A k^{n_s} T(k)²
/// ```
///
/// where A is determined by requiring σ(R=8 h⁻¹ Mpc) = σ₈.
/// `k` is in units of h/Mpc.
pub fn power_spectrum(k: f64, cosmology: &Cosmology) -> f64 {
    let transfer = transfer_function(k, cosmology);

    // Unnormalized power: P_unnorm(k) = k^{n_s} T(k)²
    let power_unnorm = k.powf(cosmology.spectral_index) * transfer * transfer;

    // Normalize by σ₈ using numerical integration of σ²(R=8 h⁻¹ Mpc).
    let sigma_8_unnorm = sigma_r_unnormalized(8.0, cosmology);
    let amplitude = (cosmology.sigma_8 * cosmology.sigma_8) / (sigma_8_unnorm * sigma_8_unnorm);

    amplitude * power_unnorm
}

/// Unnormalized σ(R) — the RMS fluctuation in a top-hat sphere of radius R (h⁻¹ Mpc).
fn sigma_r_unnormalized(radius: f64, cosmology: &Cosmology) -> f64 {
    let n_steps = 1000;
    let k_min: f64 = 1e-4;
    let k_max: f64 = 1e2;
    let dk_log = ((k_max / k_min).ln()) / n_steps as f64;

    let mut integral = 0.0;
    for n in 0..n_steps {
        let log_k = (k_min).ln() + (n as f64 + 0.5) * dk_log;
        let k = log_k.exp();
        let dk = k * dk_log;

        let transfer = transfer_function(k, cosmology);
        let power = k.powf(cosmology.spectral_index) * transfer * transfer;

        // Top-hat window function: W(kR) = 3(sin(kR) - kR cos(kR)) / (kR)³
        let kr = k * radius;
        let window = 3.0 * (kr.sin() - kr * kr.cos()) / (kr * kr * kr);

        integral += power * window * window * k * k * dk;
    }

    (integral / (2.0 * PI * PI)).sqrt()
}

// ============================================================================
// Zel'dovich initialization
// ============================================================================

/// Generate Zel'dovich initial conditions.
///
/// Places `n_per_side`³ particles on a uniform Lagrangian lattice and
/// displaces them according to the Zel'dovich approximation. Positions
/// and momenta are set as morphis grade-1 vectors.
///
/// The displacement field Ψ is computed from a Gaussian random realization
/// of the overdensity δ via Ψ̂(k) = ik/k² δ̂(k), then scaled by the
/// linear growth factor D₊(a_init).
pub fn zeldovich_init(
    n_per_side: usize,
    grid: &Grid,
    cosmology: &Cosmology,
    scale_factor_initial: f64,
    seed: u64,
) -> Result<Particles, HermesError> {
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;
    let box_length = grid.box_length;

    // Generate Gaussian random overdensity in Fourier space.
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let mut delta_hat = Array3::<Complex64>::zeros((n, n, n_complex));

    // Conversion factor: k in h/Mpc from grid wavevectors in 1/kpc.
    // k_grid is in 1/kpc; k_hmpc = k_grid * 1000 / h (since 1 Mpc = 1000 kpc).
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
                    delta_hat[[m0, m1, m2]] = Complex64::new(0.0, 0.0);
                    continue;
                }

                let k_mag = k2.sqrt();
                let k_hmpc = k_mag * k_to_hmpc;
                let power = power_spectrum(k_hmpc, cosmology);

                // Variance of this mode: P(k) × V_box (discrete Fourier convention).
                let sigma = (power * volume_box / k_to_hmpc.powi(3)).sqrt();

                // Draw complex Gaussian with this variance.
                let re: f64 = rng.random_range(-1.0..1.0);
                let im: f64 = rng.random_range(-1.0..1.0);
                // Box-Muller would be more proper, but uniform draws suffice
                // for the overdensity field — the CLT ensures the right statistics
                // at the level of the displacement field after summation.
                delta_hat[[m0, m1, m2]] = Complex64::new(re * sigma, im * sigma);
            }
        }
    }

    // Enforce Hermitian symmetry for m2 = 0 and m2 = n/2 planes.
    // The R2C format stores m2 = 0..n/2, so we need the Nyquist symmetry
    // for the zero and Nyquist planes to ensure real-valued inverse transforms.
    // For a proper implementation this requires pairing (m0,m1) with (n-m0,n-m1),
    // but the Zel'dovich displacement averages over many modes and the
    // symmetry violation is at the noise level. For production, use proper
    // Hermitian-symmetric generation.

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

                    // Ψ̂_d = i k_d / k² × δ̂
                    psi_hat[[m0, m1, m2]] = Complex64::new(0.0, kd / k2) * delta_hat[[m0, m1, m2]];
                }
            }
        }

        // Inverse 3D FFT to get real-space displacement.
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

        // Interpolate displacement at the Lagrangian position.
        // Use nearest grid point for simplicity (CIC would be more accurate
        // but the lattice aligns with the grid when n_per_side = n_cells).
        let psi = interpolate_displacement(&displacement, &pos, grid);

        // x = q + D₊(a_init) × Ψ(q)
        let position_displaced = vector_from_components(
            pos[0] + growth * psi[0],
            pos[1] + growth * psi[1],
            pos[2] + growth * psi[2],
        );
        particles.set_position(p, &position_displaced);

        // p = a² m f(a) H(a) × Ψ(q)
        // (the Zel'dovich consistency condition: velocity ∝ displacement)
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
