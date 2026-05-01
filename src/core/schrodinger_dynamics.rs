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
use num_complex::Complex64;

use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::error::HermesError;
use crate::physics::constants::G as GRAV;
use crate::physics::cosmology::Cosmology;
use crate::physics::spectral::{fft_3d_dyn, ifft_3d_dyn};

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
    let scalar_hat = fft_3d_dyn(&alpha.scalar, n);
    let pseudo_hat = fft_3d_dyn(&alpha.pseudoscalar, n);

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
    alpha.scalar = ifft_3d_dyn(&result_scalar_hat, n);
    alpha.pseudoscalar = ifft_3d_dyn(&result_pseudo_hat, n);
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
    // Density: ρ = m (a² + b��)
    let rho = alpha.density(mass);

    // Poisson source: ∇²Φ = 4πG a² (ρ - ρ̄)
    let rho_bar_field = Field::scalar_field(grid, metric::euclidean::<3>(), |_| density_mean);
    let poisson_coupling = 4.0 * PI * GRAV * scale_factor * scale_factor;
    let source = &(&rho - &rho_bar_field) * poisson_coupling;
    let phi = source.laplacian_inverse();

    // Phase rotation: theta = -m * Phi * dt / ell
    let angle = &phi * (-mass * dt / ell);
    *alpha = alpha.rotate_phase(&angle);
}

// ============================================================================
// Madelung velocity extraction
// ============================================================================

/// Extract the velocity field from the dark matter field via the Madelung form.
///
/// Delegates to morphis's `madelung_velocity`, which computes
/// v_d = (ν / |α|²) (a ∂_d b - b ∂_d a) spectrally for each direction,
/// avoiding the phase branch cut of ∇ arg(α).
///
/// Returns raw arrays for compatibility with the snapshot and
/// visualization pipeline. The morphis Field<3> is extracted into
/// three scalar component arrays.
pub fn extract_velocity(
    alpha: &morphis::even_field::EvenField<3>,
    _grid: &morphis::grid::Grid<3>,
    ell: f64,
    mass: f64,
) -> [ndarray::ArrayD<f64>; 3] {
    let nu = ell / mass;
    let v_field = alpha.madelung_velocity(nu);
    let n = alpha.grid.n_cells;

    std::array::from_fn(|d| {
        let component = v_field.component_field(&[d]);
        let mut out = ndarray::ArrayD::zeros(ndarray::IxDyn(&[n, n, n]));
        out.assign(&component.data);
        out
    })
}
