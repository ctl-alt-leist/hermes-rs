/// Gravity solver: Poisson equation from aggregated sector densities.
///
/// Takes the total mass density deposited by all sectors, subtracts
/// the cosmological mean, and solves the Poisson equation for the
/// gravitational potential. The solver uses morphis `laplacian_inverse`
/// for the spectral inversion.
use std::f64::consts::PI;

use morphis::field::Field;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric;

use crate::engine::sector::Potential;
use crate::error::HermesError;
use crate::physics::constants::G as GRAV;

/// FFT-based Poisson solver for gravitational coupling.
///
/// Owns the morphis grid needed to construct the background density
/// field. Reused across timesteps.
pub struct GravitySolver {
    /// Morphis grid for spectral operations.
    grid: MorphisGrid<3>,
}

impl GravitySolver {
    /// Create a gravity solver for the given grid.
    pub fn new(grid: MorphisGrid<3>) -> Self {
        Self { grid }
    }

    /// Solve the Poisson equation for the gravitational potential.
    ///
    /// Given the total mass density from all sectors, computes:
    ///   source = (rho - rho_mean) * 4 pi G a^2
    ///   phi = laplacian_inverse(source)
    pub fn solve(
        &self,
        total_density: &Field<3>,
        density_mean: f64,
        scale_factor: f64,
    ) -> Result<Potential, HermesError> {
        let rho_bar = Field::scalar_field(&self.grid, metric::euclidean::<3>(), |_| density_mean);

        let poisson_coupling = 4.0 * PI * GRAV * scale_factor * scale_factor;
        let source = &(total_density - &rho_bar) * poisson_coupling;
        let phi = source.laplacian_inverse();

        Ok(Potential { phi })
    }
}
