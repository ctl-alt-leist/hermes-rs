//! Schrodinger-Poisson dynamics: split-step spectral integrator.
//!
//! Evolves the dark matter field α in G^+ (even subalgebra) under
//! self-gravity via the symmetric split-step method:
//!
//!   1. Kinetic half-step (Fourier space): phase rotation by -ℓ |k|² dt / (2m a²)
//!   2. Potential full step (real space): Poisson solve, phase rotation by -m Φ dt / ℓ
//!   3. Kinetic half-step (repeat step 1)
//!
//! The kinetic step uses FFT to transform the field components
//! to Fourier space, applies the k-dependent phase rotation, and
//! transforms back. The potential step uses morphis's laplacian_inverse
//! for the Poisson solve.

use std::f64::consts::PI;

use morphis::field::Field;
use morphis::metric;
use ndarray::Array3;
use ndrustfft::{FftHandler, R2cFftHandler, ndfft, ndfft_r2c, ndifft, ndifft_r2c};
use num_complex::Complex64;

use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::error::HermesError;
use crate::physics::constants::G as GRAV;
use crate::physics::cosmology::Cosmology;

/// Schrodinger-Poisson dynamics for wavefunction dark matter.
pub struct SchrodingerPoissonDynamics;

impl SchrodingerPoissonDynamics {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SchrodingerPoissonDynamics {
    fn default() -> Self {
        Self::new()
    }
}

impl Dynamics for SchrodingerPoissonDynamics {
    fn step(
        &mut self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor_prev: f64,
        scale_factor_next: f64,
    ) -> Result<(), HermesError> {
        let fields = content.fields_mut().ok_or_else(|| {
            HermesError::Config("Schrodinger dynamics requires field content".to_string())
        })?;

        let alpha = fields.alpha.as_mut().ok_or_else(|| {
            HermesError::Config("Schrodinger dynamics requires alpha field".to_string())
        })?;

        let scale_factor = (scale_factor_prev + scale_factor_next) / 2.0;
        let dt = (scale_factor_next - scale_factor_prev)
            / (scale_factor * cosmology.hubble_parameter(scale_factor));

        let ell = fields.params.smoothing_length;
        let mass = fields.params.mass_alpha;
        let density_mean = cosmology.density_matter();

        // 1. Kinetic half-step in Fourier space.
        kinetic_step(alpha, &fields.grid, ell, mass, scale_factor, dt / 2.0);

        // 2. Potential full step in real space.
        potential_step(
            alpha,
            &fields.grid,
            ell,
            mass,
            density_mean,
            scale_factor,
            dt,
        );

        // 3. Kinetic half-step in Fourier space.
        kinetic_step(alpha, &fields.grid, ell, mass, scale_factor, dt / 2.0);

        Ok(())
    }
}

/// Kinetic half-step: FFT α, rotate by exp(I θ(k)), IFFT.
///
/// θ(k) = -ℓ |k|² dt / (2m a²)
pub fn kinetic_step(
    alpha: &mut morphis::even_field::EvenField<3>,
    grid: &morphis::grid::Grid<3>,
    ell: f64,
    mass: f64,
    scale_factor: f64,
    dt: f64,
) {
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;

    // FFT both components to Fourier space.
    let scalar_hat = fft_3d(&alpha.scalar, n);
    let pseudo_hat = fft_3d(&alpha.pseudoscalar, n);

    // Apply phase rotation in k-space: (a + bI) * (cos theta + sin theta I)
    let mut result_scalar_hat = Array3::<Complex64>::zeros((n, n, n_complex));
    let mut result_pseudo_hat = Array3::<Complex64>::zeros((n, n, n_complex));

    for m0 in 0..n {
        let kx = grid.wavenumber(m0);
        for m1 in 0..n {
            let ky = grid.wavenumber(m1);
            for m2 in 0..n_complex {
                let kz = grid.wavenumber(m2);
                let k2 = kx * kx + ky * ky + kz * kz;

                let theta = -ell * k2 * dt / (2.0 * mass * scale_factor * scale_factor);
                let cos_t = theta.cos();
                let sin_t = theta.sin();

                let a = scalar_hat[[m0, m1, m2]];
                let b = pseudo_hat[[m0, m1, m2]];

                // (a + bI)(cos + sin I) = (a cos - b sin) + (a sin + b cos) I
                result_scalar_hat[[m0, m1, m2]] = a * cos_t - b * sin_t;
                result_pseudo_hat[[m0, m1, m2]] = a * sin_t + b * cos_t;
            }
        }
    }

    // IFFT back to real space.
    alpha.scalar = ifft_3d(&result_scalar_hat, n);
    alpha.pseudoscalar = ifft_3d(&result_pseudo_hat, n);
}

/// Potential full step: compute density, Poisson solve, phase rotation.
pub fn potential_step(
    alpha: &mut morphis::even_field::EvenField<3>,
    grid: &morphis::grid::Grid<3>,
    ell: f64,
    mass: f64,
    density_mean: f64,
    scale_factor: f64,
    dt: f64,
) {
    // Density: rho = m * (a^2 + b^2)
    let rho = alpha.density(mass);

    // Overdensity: delta = rho / rho_bar - 1
    let mut overdensity = &rho * (1.0 / density_mean);
    // Subtract 1 pointwise.
    let _n = grid.n_cells;
    let ones = Field::scalar_field(grid, metric::euclidean::<3>(), |_| 1.0);
    overdensity = &overdensity - &ones;

    // Poisson solve: Phi = (4 pi G rho_bar a^2 * delta).laplacian_inverse()
    let prefactor = 4.0 * PI * GRAV * density_mean * scale_factor * scale_factor;
    let source = &overdensity * prefactor;
    let phi = source.laplacian_inverse();

    // Phase rotation: theta = -m * Phi * dt / ell
    let angle = &phi * (-mass * dt / ell);
    *alpha = alpha.rotate_phase(&angle);
}

// ============================================================================
// Madelung velocity extraction
// ============================================================================

/// Extract the velocity field from the dark matter field via the field-gradient form.
///
/// Computes v_d = (ℓ / (m |α|²)) Im(ᾱ ∇_d α) for each spatial direction,
/// where ᾱ is the even-subalgebra conjugate (s - pI) and ∇_d α is
/// the spectral gradient along direction d.
///
/// This avoids the phase branch cut that makes the naive
/// v = (ℓ/m) ∇ arg(α) form unreliable when the phase exceeds 2π.
pub fn extract_velocity(
    alpha: &morphis::even_field::EvenField<3>,
    grid: &morphis::grid::Grid<3>,
    ell: f64,
    mass: f64,
) -> [ndarray::ArrayD<f64>; 3] {
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;

    let scalar_hat = fft_3d(&alpha.scalar, n);
    let pseudo_hat = fft_3d(&alpha.pseudoscalar, n);

    let mut velocity: [ndarray::ArrayD<f64>; 3] =
        std::array::from_fn(|_| ndarray::ArrayD::zeros(ndarray::IxDyn(&[n, n, n])));

    #[allow(clippy::needless_range_loop)]
    for d in 0..3 {
        // Spectral derivative: multiply by i * k_d.
        let mut ds_hat = Array3::<Complex64>::zeros((n, n, n_complex));
        let mut dp_hat = Array3::<Complex64>::zeros((n, n, n_complex));

        for m0 in 0..n {
            let kx = grid.wavenumber(m0);
            for m1 in 0..n {
                let ky = grid.wavenumber(m1);
                for m2 in 0..n_complex {
                    let kz = grid.wavenumber(m2);
                    let kd = match d {
                        0 => kx,
                        1 => ky,
                        _ => kz,
                    };
                    let ik = Complex64::new(0.0, kd);
                    ds_hat[[m0, m1, m2]] = scalar_hat[[m0, m1, m2]] * ik;
                    dp_hat[[m0, m1, m2]] = pseudo_hat[[m0, m1, m2]] * ik;
                }
            }
        }

        let ds = ifft_3d(&ds_hat, n);
        let dp = ifft_3d(&dp_hat, n);

        // v_d = (ell / (m |alpha|^2)) * (scalar * dp - pseudo * ds)
        let v_d = velocity[d]
            .as_slice_mut()
            .expect("velocity array not contiguous");
        let s = alpha.scalar.as_slice().expect("scalar not contiguous");
        let p = alpha
            .pseudoscalar
            .as_slice()
            .expect("pseudo not contiguous");
        let ds_slice = ds.as_slice().expect("ds not contiguous");
        let dp_slice = dp.as_slice().expect("dp not contiguous");

        for k in 0..v_d.len() {
            let mod_sq = s[k] * s[k] + p[k] * p[k];
            let im_part = s[k] * dp_slice[k] - p[k] * ds_slice[k];
            v_d[k] = (ell / mass) * im_part / mod_sq;
        }
    }

    velocity
}

// ============================================================================
// FFT helpers for EvenField components
// ============================================================================

/// Forward 3D R2C FFT on an ndarray.
pub fn fft_3d(data: &ndarray::ArrayD<f64>, n: usize) -> Array3<Complex64> {
    let n_complex = n / 2 + 1;

    // Reshape to Array3 for ndrustfft.
    let data_3d = data
        .view()
        .into_shape_with_order((n, n, n))
        .expect("data shape mismatch");

    let handler_r2c = R2cFftHandler::new(n);
    let handler_c2c_1 = FftHandler::new(n);
    let handler_c2c_0 = FftHandler::new(n);

    let mut complex = Array3::<Complex64>::zeros((n, n, n_complex));
    ndfft_r2c(&data_3d, &mut complex, &handler_r2c, 2);

    let mut scratch = complex.clone();
    ndfft(&complex, &mut scratch, &handler_c2c_1, 1);
    complex.assign(&scratch);
    ndfft(&complex, &mut scratch, &handler_c2c_0, 0);
    complex.assign(&scratch);

    complex
}

/// Inverse 3D C2R FFT, returning an ndarray::ArrayD.
pub fn ifft_3d(complex: &Array3<Complex64>, n: usize) -> ndarray::ArrayD<f64> {
    let handler_c2c_0 = FftHandler::new(n);
    let handler_c2c_1 = FftHandler::new(n);
    let handler_r2c = R2cFftHandler::new(n);

    let mut work = complex.clone();
    let mut scratch = work.clone();

    ndifft(&work, &mut scratch, &handler_c2c_0, 0);
    work.assign(&scratch);
    ndifft(&work, &mut scratch, &handler_c2c_1, 1);
    work.assign(&scratch);

    let mut real = Array3::<f64>::zeros((n, n, n));
    ndifft_r2c(&work, &mut real, &handler_r2c, 2);

    real.into_dyn()
}
