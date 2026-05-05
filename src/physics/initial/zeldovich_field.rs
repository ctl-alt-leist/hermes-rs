//! Cosmic web field initial conditions.
//!
//! Density and velocity fields converted to the even-subalgebra
//! representation via the inverse Madelung transformation:
//!
//!   α(x) = sqrt(ρ(x) / m) * exp(I * m * φ_v(x) / ℓ)

use morphis::even_field::EvenField;
use morphis::field::Field;
use std::f64::consts::PI;

use morphis::metric;
use ndarray::Array3;
use num_complex::Complex64;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use super::zeldovich::power_spectrum;
use crate::core::content::FieldParams;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid as HermesGrid;
use crate::physics::spectral::{fft_3d, ifft_3d};

/// Generate Zel'dovich initial conditions as a wavefunction.
pub fn zeldovich_wavefunction(
    grid: &HermesGrid,
    cosmology: &Cosmology,
    params: &FieldParams,
    scale_factor_initial: f64,
    perturbation_amplitude: f64,
    seed: u64,
) -> Result<EvenField<3>, HermesError> {
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;
    let box_length = grid.box_length;
    let cell_length = grid.cell_length;

    let ell = params.smoothing_length;
    let mass = params.mass_alpha;

    let density_mean = cosmology.density_matter();
    let growth = cosmology.growth_factor(scale_factor_initial);
    let growth_rate = cosmology.growth_rate(scale_factor_initial);
    let hubble_a = cosmology.hubble_parameter(scale_factor_initial);

    // Generate overdensity in Fourier space.
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
                    delta_hat[[m0, m1, m2]] = Complex64::new(0.0, 0.0);
                    continue;
                }

                let k_hmpc = k2.sqrt() * k_to_hmpc;
                let power = power_spectrum(k_hmpc, cosmology);
                let volume_hmpc3 = volume_box / k_to_hmpc.powi(3);
                let sigma = (power * volume_hmpc3).sqrt();

                let re: f64 = rng.random_range(-1.0..1.0);
                let im: f64 = rng.random_range(-1.0..1.0);
                delta_hat[[m0, m1, m2]] = Complex64::new(re * sigma, im * sigma);
            }
        }
    }

    // Compute displacement field: Psi_d(k) = i k_d / k^2 * delta(k).
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

        *displacement_d = ifft_3d(&psi_hat, n);
    }

    // Compute overdensity: delta = -D+(a) * div(Psi), then normalize
    // to the target perturbation amplitude. The natural linear overdensity
    // at high z is tiny at grid scales; normalizing to a target RMS makes
    // the initial perturbation amplitude an explicit input rather than an
    // artifact of the starting redshift.
    let mut div_psi = Array3::<f64>::zeros((n, n, n));
    #[allow(clippy::needless_range_loop)]
    for d in 0..3 {
        let psi_d_hat = fft_3d(&displacement[d], n);
        let mut deriv_hat = Array3::<Complex64>::zeros((n, n, n_complex));

        for m0 in 0..n {
            let kx = grid.wavevector_component(m0);
            for m1 in 0..n {
                let ky = grid.wavevector_component(m1);
                for m2 in 0..n_complex {
                    let kz = grid.wavevector_component(m2);
                    let kd = match d {
                        0 => kx,
                        1 => ky,
                        _ => kz,
                    };
                    deriv_hat[[m0, m1, m2]] = psi_d_hat[[m0, m1, m2]] * Complex64::new(0.0, kd);
                }
            }
        }

        let deriv_real = ifft_3d(&deriv_hat, n);
        div_psi += &deriv_real;
    }

    // Normalize the overdensity to the target perturbation amplitude.
    // The natural Zel'dovich delta = -D+(a) * div(Psi) is tiny at high z;
    // normalizing to a target RMS replaces the old hardcoded 50000x boost.
    let delta_raw = div_psi.mapv(|d| -growth * d);
    let delta_rms = (delta_raw.iter().map(|d| d * d).sum::<f64>() / delta_raw.len() as f64).sqrt();
    let norm = perturbation_amplitude / delta_rms.max(1e-30);

    // Compute velocity potential: phi_v_hat = v_scale * delta_hat / k^2.
    //
    // The Madelung phase is S = (m / ell) * phi_v, NOT (m / ell) * v . x.
    // Using v . x produces an unbounded, wildly oscillating phase that the
    // kinetic step immediately disperses. The velocity potential is the
    // scalar whose gradient gives the velocity field, and it is smooth
    // and periodic on the box.
    let v_scale = growth * growth_rate * hubble_a;
    let mut phi_v_hat = Array3::<Complex64>::zeros((n, n, n_complex));
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
                phi_v_hat[[m0, m1, m2]] = delta_hat[[m0, m1, m2]] * (norm * v_scale / k2);
            }
        }
    }
    let phi_v = ifft_3d(&phi_v_hat, n);

    // Build morphis fields for the inverse Madelung transform.
    let morphis_grid = morphis::grid::Grid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let nu = ell / mass;

    let rho_field = Field::scalar_field(&morphis_grid, g, |x| {
        let m0 = ((x[0] / cell_length) as usize).min(n - 1);
        let m1 = ((x[1] / cell_length) as usize).min(n - 1);
        let m2 = ((x[2] / cell_length) as usize).min(n - 1);
        let delta = norm * delta_raw[[m0, m1, m2]];
        density_mean * (1.0 + delta).max(1e-10)
    });

    let phi_v_field = Field::scalar_field(&morphis_grid, g, |x| {
        let m0 = ((x[0] / cell_length) as usize).min(n - 1);
        let m1 = ((x[1] / cell_length) as usize).min(n - 1);
        let m2 = ((x[2] / cell_length) as usize).min(n - 1);
        phi_v[[m0, m1, m2]]
    });

    let result = EvenField::madelung_inverse(&rho_field, &phi_v_field, mass, nu);

    Ok(result)
}

/// Random multi-scale density field converted to wavefunction.
///
/// Generates Fourier modes with Gaussian random amplitudes and a
/// band-pass spectrum: suppressed below k_min (removes the box-scale
/// mode that produces a trivial half-and-half split) and above k_max
/// (removes grid-scale noise). The density amplitude is normalized
/// to a target delta_rms, producing visible structure immediately.
///
/// Velocity is derived from the density gradient (Zel'dovich-like
/// relation v ~ -grad(delta) / k^2), scaled to be dynamically active.
pub fn random_density_field(
    grid: &HermesGrid,
    cosmology: &Cosmology,
    params: &FieldParams,
    scale_factor_initial: f64,
    perturbation_amplitude: f64,
    band_pass: [f64; 2],
    seed: u64,
) -> EvenField<3> {
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;
    let box_length = grid.box_length;
    let cell_length = grid.cell_length;
    let ell = params.smoothing_length;
    let mass = params.mass_alpha;
    let density_mean = cosmology.density_matter();

    let mut rng = ChaCha20Rng::seed_from_u64(seed);

    let k_fundamental = 2.0 * PI / box_length;

    let k_min = band_pass[0] * k_fundamental;
    let k_max = band_pass[1] * PI * n as f64 / box_length;

    // Generate overdensity in Fourier space with Gaussian amplitudes
    // and a red spectrum P(k) ~ k^-2 that concentrates power at large
    // scales, producing visible clumps rather than fine-grained noise.
    let mut delta_hat = Array3::<Complex64>::zeros((n, n, n_complex));

    for m0 in 0..n {
        let kx = grid.wavevector_component(m0);
        for m1 in 0..n {
            let ky = grid.wavevector_component(m1);
            for m2 in 0..n_complex {
                let kz = grid.wavevector_component(m2);
                let k2 = kx * kx + ky * ky + kz * kz;
                let k = k2.sqrt();

                if k < k_min || k > k_max {
                    continue;
                }

                // Red spectrum: large-scale modes dominate.
                let amplitude = k_fundamental / k;

                // Gaussian random amplitudes (Box-Muller).
                let u1: f64 = rng.random_range(1e-10..1.0);
                let u2: f64 = rng.random_range(0.0..2.0 * PI);
                let gauss_r = (-2.0 * u1.ln()).sqrt() * u2.cos();
                let gauss_i = (-2.0 * u1.ln()).sqrt() * u2.sin();

                delta_hat[[m0, m1, m2]] = Complex64::new(gauss_r * amplitude, gauss_i * amplitude);
            }
        }
    }

    // IFFT to real space, then normalize to target amplitude.
    let delta_raw = ifft_3d(&delta_hat, n);
    let delta_rms = (delta_raw.iter().map(|d| d * d).sum::<f64>() / delta_raw.len() as f64).sqrt();
    let target_rms = perturbation_amplitude;
    let norm = target_rms / delta_rms.max(1e-30);

    // Compute the velocity potential for the Zel'dovich growing mode.
    //
    // The velocity is v = grad(phi_v), where the velocity potential is:
    //   phi_v_hat(k) = a H(a) f(a) * delta_hat(k) / k^2
    //
    // The Madelung phase is S = (m / ell) * phi_v, NOT (m / ell) * v . x.
    // Using v . x instead of the velocity potential produces an unbounded,
    // wildly oscillating phase that the kinetic step immediately disperses.
    let delta_hat_normalized = delta_hat.mapv(|c| c * norm);

    let a = scale_factor_initial;
    let v_scale = a * cosmology.hubble_parameter(a) * cosmology.growth_rate(a);

    let mut phi_v_hat = Array3::<Complex64>::zeros((n, n, n_complex));
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
                // phi_v_hat = v_scale * delta_hat / k^2
                phi_v_hat[[m0, m1, m2]] = delta_hat_normalized[[m0, m1, m2]] * (v_scale / k2);
            }
        }
    }

    let phi_v = ifft_3d(&phi_v_hat, n);

    // Build morphis fields for the inverse Madelung transform.
    let morphis_grid = morphis::grid::Grid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let nu = ell / mass;

    let rho_field = Field::scalar_field(&morphis_grid, g, |x| {
        let m0 = ((x[0] / cell_length) as usize).min(n - 1);
        let m1 = ((x[1] / cell_length) as usize).min(n - 1);
        let m2 = ((x[2] / cell_length) as usize).min(n - 1);
        let delta = norm * delta_raw[[m0, m1, m2]];
        density_mean * (1.0 + delta).max(1e-10)
    });

    let phi_v_field = Field::scalar_field(&morphis_grid, g, |x| {
        let m0 = ((x[0] / cell_length) as usize).min(n - 1);
        let m1 = ((x[1] / cell_length) as usize).min(n - 1);
        let m2 = ((x[2] / cell_length) as usize).min(n - 1);
        phi_v[[m0, m1, m2]]
    });

    EvenField::madelung_inverse(&rho_field, &phi_v_field, mass, nu)
}
