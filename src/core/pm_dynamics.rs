//! Particle-mesh dynamics: KDK leapfrog with FFT-Poisson gravity.
//!
//! Delegates to `PoissonGravity` for the gravitational force chain
//! and to the KDK integrator structure for time stepping.

use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::engine::coupling::poisson::PoissonGravity;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::integrator::midpoint;
use crate::physics::poisson::PoissonSolver;

/// Particle-mesh N-body dynamics.
///
/// Delegates to `PoissonGravity` for the full KDK step.
pub struct ParticleMeshDynamics {
    gravity: PoissonGravity,
}

impl ParticleMeshDynamics {
    /// Create PM dynamics for a given grid.
    pub fn new(grid: Grid) -> Self {
        Self {
            gravity: PoissonGravity::new(grid),
        }
    }

    /// Access the grid.
    pub fn grid(&self) -> &Grid {
        self.gravity.grid()
    }

    /// Access the solver (for diagnostics that need it).
    pub fn solver_mut(&mut self) -> &mut PoissonSolver {
        self.gravity.solver_mut()
    }
}

impl Dynamics for ParticleMeshDynamics {
    fn step(
        &mut self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor_prev: f64,
        scale_factor_next: f64,
    ) -> Result<(), HermesError> {
        let scale_factor_mid = midpoint(scale_factor_prev, scale_factor_next);

        self.gravity.step_particles_kdk(
            content,
            cosmology,
            scale_factor_prev,
            scale_factor_mid,
            scale_factor_next,
        )
    }
}
