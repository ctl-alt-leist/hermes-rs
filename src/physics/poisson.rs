//! FFT-based Poisson solver for periodic gravity.
//!
//! Solves the Poisson equation on a periodic cubic grid using real-to-complex
//! FFTs via `ndrustfft`. The solver precomputes the discrete Green's function
//! and reuses FFT plans and workspace arrays across calls.
//!
//! The force chain for one step is:
//!
//! ```text
//! overdensity δ(x) → FFT → δ̂(k) × G(k) → IFFT → ϕ(x) → ∇ϕ → F(x)
//! ```

use std::f64::consts::PI;

use ndarray::Array3;
use ndrustfft::{FftHandler, R2cFftHandler, ndfft, ndfft_r2c, ndifft, ndifft_r2c};
use num_complex::Complex64;

use crate::physics::constants::G as GRAV;
use crate::physics::field::{ScalarField, VectorField};
use crate::physics::grid::Grid;

/// FFT-based Poisson solver with precomputed Green's function.
pub struct PoissonSolver {
    /// Discrete Green's function in Fourier space, shape [N, N, N/2+1].
    green: Array3<f64>,
    /// Grid geometry.
    grid: Grid,
    /// R2C handler for axis 2 (the real-to-complex axis).
    handler_r2c: R2cFftHandler<f64>,
    /// C2C handler for axis 1.
    handler_c2c_1: FftHandler<f64>,
    /// C2C handler for axis 0.
    handler_c2c_0: FftHandler<f64>,
}

impl PoissonSolver {
    /// Create a new solver for the given grid.
    ///
    /// Precomputes the discrete Green's function using the finite-difference
    /// Laplacian so that the force is consistent with a second-order
    /// centered gradient stencil.
    pub fn new(grid: &Grid) -> Self {
        let n = grid.n_cells;
        let n_complex = n / 2 + 1;
        let h = grid.cell_length;

        let mut green = Array3::zeros((n, n, n_complex));
        let factor = (2.0 / h) * (2.0 / h);

        for m0 in 0..n {
            let kx_h = PI * freq(m0, n) / n as f64;
            let sx2 = kx_h.sin().powi(2);

            for m1 in 0..n {
                let ky_h = PI * freq(m1, n) / n as f64;
                let sy2 = ky_h.sin().powi(2);

                for m2 in 0..n_complex {
                    let kz_h = PI * freq(m2, n) / n as f64;
                    let sz2 = kz_h.sin().powi(2);

                    let k2_discrete = factor * (sx2 + sy2 + sz2);

                    green[[m0, m1, m2]] = if k2_discrete.abs() < 1e-30 {
                        0.0
                    } else {
                        -1.0 / k2_discrete
                    };
                }
            }
        }

        Self {
            green,
            grid: grid.clone(),
            handler_r2c: R2cFftHandler::new(n),
            handler_c2c_1: FftHandler::new(n),
            handler_c2c_0: FftHandler::new(n),
        }
    }

    /// Reference to the grid geometry.
    pub fn grid_ref(&self) -> &Grid {
        &self.grid
    }

    /// Reference to the precomputed Green's function.
    pub fn green_function(&self) -> &Array3<f64> {
        &self.green
    }

    /// Solve for the gravitational force field from an overdensity.
    ///
    /// Given the overdensity δ = ρ/ρ̄ - 1 on the grid, returns the
    /// gravitational force field F = -∇ϕ with the cosmological prefactor
    /// 4πG ρ̄ a² included, so that the force can be used directly in the
    /// kick step: dp/dt = F/a, with kick_factor absorbing the 1/a.
    ///
    /// The `density_mean` parameter is ρ̄_m in M_☉/kpc³ (comoving),
    /// and `scale_factor` is the current value of a.
    pub fn solve(
        &mut self,
        overdensity: &ScalarField,
        density_mean: f64,
        scale_factor: f64,
    ) -> VectorField {
        let n = self.grid.n_cells;
        let n_complex = n / 2 + 1;

        // Forward 3D R2C FFT: real [N,N,N] → complex [N,N,N/2+1].
        //   1) R2C along axis 2 (innermost)
        //   2) C2C along axis 1
        //   3) C2C along axis 0
        let mut overdensity_hat = Array3::<Complex64>::zeros((n, n, n_complex));
        ndfft_r2c(
            &overdensity.data,
            &mut overdensity_hat,
            &self.handler_r2c,
            2,
        );

        let mut scratch = overdensity_hat.clone();
        ndfft(&overdensity_hat, &mut scratch, &self.handler_c2c_1, 1);
        overdensity_hat.assign(&scratch);
        ndfft(&overdensity_hat, &mut scratch, &self.handler_c2c_0, 0);
        overdensity_hat.assign(&scratch);

        // Multiply by Green's function with physical prefactor.
        // ϕ̂(k) = 4πG ρ̄ a² × G(k) × δ̂(k)
        let prefactor = 4.0 * PI * GRAV * density_mean * scale_factor * scale_factor;
        let mut potential_hat = overdensity_hat;

        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n_complex {
                    potential_hat[[m0, m1, m2]] *= prefactor * self.green[[m0, m1, m2]];
                }
            }
        }

        // Compute force components: F_d = -ik_d × ϕ̂(k), then IFFT.
        let mut force = VectorField::zeros(&self.grid);

        for d in 0..3 {
            // F̂_d(k) = -i k_d × ϕ̂(k)
            let mut force_hat = potential_hat.clone();

            for m0 in 0..n {
                let k0 = self.grid.wavevector_component(m0);
                for m1 in 0..n {
                    let k1 = self.grid.wavevector_component(m1);
                    for m2 in 0..n_complex {
                        let k2 = self.grid.wavevector_component(m2);
                        let kd = match d {
                            0 => k0,
                            1 => k1,
                            _ => k2,
                        };
                        force_hat[[m0, m1, m2]] *= Complex64::new(0.0, -kd);
                    }
                }
            }

            // Inverse 3D C2R FFT: reverse order.
            //   1) C2C inverse along axis 0
            //   2) C2C inverse along axis 1
            //   3) C2R inverse along axis 2
            let mut scratch = force_hat.clone();
            ndifft(&force_hat, &mut scratch, &self.handler_c2c_0, 0);
            force_hat.assign(&scratch);
            ndifft(&force_hat, &mut scratch, &self.handler_c2c_1, 1);
            force_hat.assign(&scratch);

            let mut force_component = Array3::<f64>::zeros((n, n, n));
            ndifft_r2c(&force_hat, &mut force_component, &self.handler_r2c, 2);

            force.data[d] = force_component;
        }

        force
    }
}

/// Convert a grid index to a frequency index with negative wrapping.
fn freq(m: usize, n: usize) -> f64 {
    if m <= n / 2 {
        m as f64
    } else {
        m as f64 - n as f64
    }
}
