//! Particle-mesh dynamics: KDK leapfrog with FFT-Poisson gravity.
//!
//! Wraps the existing integrator and Poisson solver into a `Dynamics`
//! implementation for particle content.

use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::integrator::{midpoint, step_kdk};
use crate::physics::poisson::PoissonSolver;

/// Particle-mesh N-body dynamics.
///
/// Owns the PM grid and Poisson solver. Delegates to the existing
/// KDK leapfrog integrator for the actual step.
pub struct ParticleMeshDynamics {
    pub grid: Grid,
    pub solver: PoissonSolver,
    forces_prev: Option<crate::physics::cic::ParticleForces>,
}

impl ParticleMeshDynamics {
    /// Create PM dynamics for a given grid.
    pub fn new(grid: Grid) -> Self {
        let solver = PoissonSolver::new(&grid);

        Self {
            grid,
            solver,
            forces_prev: None,
        }
    }

    /// Access the grid.
    pub fn grid(&self) -> &Grid {
        &self.grid
    }

    /// Access the solver (for diagnostics that need it).
    pub fn solver_mut(&mut self) -> &mut PoissonSolver {
        &mut self.solver
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
        let particles = content.particles_mut().ok_or_else(|| {
            HermesError::Config("PM dynamics requires particle content".to_string())
        })?;

        let scale_factor_mid = midpoint(scale_factor_prev, scale_factor_next);

        let forces = step_kdk(
            particles,
            &mut self.solver,
            &self.grid,
            cosmology,
            scale_factor_prev,
            scale_factor_mid,
            scale_factor_next,
            self.forces_prev.as_ref(),
        )?;

        self.forces_prev = Some(forces);

        Ok(())
    }
}
