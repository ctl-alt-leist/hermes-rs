/// Poisson gravity coupling: instantaneous gravitational interaction.
///
/// Collects density from all participating species onto a shared grid,
/// solves the Poisson equation once, and applies the gravitational
/// potential back to each species in its own representation:
///
///   - Particles: force interpolation + momentum kick
///   - Fields: potential phase rotation
///
/// This module replaces the gravity code previously embedded in both
/// ParticleMeshDynamics and SchrodingerPoissonDynamics.
use std::f64::consts::PI;

use morphis::field::Field;
use morphis::metric;

use crate::core::content::Content;
use crate::engine::coupling::Coupling;
use crate::engine::state::SimulationState;
use crate::error::HermesError;
use crate::physics::cic::{ParticleForces, assign_density, interpolate_force};
use crate::physics::constants::G as GRAV;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::poisson::PoissonSolver;

/// Poisson gravity module.
///
/// Owns the PM grid and Poisson solver. Applies gravitational forces
/// to any combination of particles and fields.
pub struct PoissonGravity {
    /// Grid for CIC assignment and Poisson solving.
    grid: Grid,
    /// FFT-based Poisson solver (reused across steps).
    solver: PoissonSolver,
    /// Cached particle forces from the previous half-kick.
    forces_prev: Option<ParticleForces>,
}

impl PoissonGravity {
    /// Create a Poisson gravity module for a given grid.
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

    /// Access the solver mutably (for diagnostics).
    pub fn solver_mut(&mut self) -> &mut PoissonSolver {
        &mut self.solver
    }

    /// Apply a gravity half-step to particles: compute forces and kick.
    ///
    /// On the first call (or when forces are not cached), computes forces
    /// from the current particle positions. On subsequent calls, reuses
    /// cached forces from the previous closing half-kick.
    pub fn kick_particles(
        &mut self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor_from: f64,
        scale_factor_to: f64,
    ) -> Result<(), HermesError> {
        let particles = content.particles_mut().ok_or_else(|| {
            HermesError::Config("gravity kick requires particle content".to_string())
        })?;

        let kick_factor = cosmology.kick_factor(scale_factor_from, scale_factor_to);

        let forces = match self.forces_prev.take() {
            Some(f) => f,
            None => compute_particle_forces(
                particles,
                &mut self.solver,
                &self.grid,
                cosmology,
                scale_factor_from,
            ),
        };

        crate::physics::integrator::kick(particles, &forces, kick_factor);
        self.forces_prev = Some(forces);

        Ok(())
    }

    /// Recompute forces at current particle positions and kick.
    ///
    /// Used after a drift when particle positions have changed and
    /// the cached forces are stale.
    pub fn recompute_and_kick_particles(
        &mut self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor: f64,
        scale_factor_from: f64,
        scale_factor_to: f64,
    ) -> Result<(), HermesError> {
        let particles = content.particles_mut().ok_or_else(|| {
            HermesError::Config("gravity kick requires particle content".to_string())
        })?;

        let forces = compute_particle_forces(
            particles,
            &mut self.solver,
            &self.grid,
            cosmology,
            scale_factor,
        );

        let kick_factor = cosmology.kick_factor(scale_factor_from, scale_factor_to);
        crate::physics::integrator::kick(particles, &forces, kick_factor);
        self.forces_prev = Some(forces);

        Ok(())
    }

    /// Drift particles using their current momenta.
    pub fn drift_particles(
        &self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor_from: f64,
        scale_factor_to: f64,
    ) -> Result<(), HermesError> {
        let particles = content.particles_mut().ok_or_else(|| {
            HermesError::Config("gravity drift requires particle content".to_string())
        })?;

        let drift_factor = cosmology.drift_factor(scale_factor_from, scale_factor_to);
        crate::physics::integrator::drift(particles, drift_factor, &self.grid);

        Ok(())
    }

    /// Apply gravity to field content: Poisson solve and phase rotation.
    ///
    /// Computes the gravitational potential from the field density and
    /// rotates the wavefunction phase by -m Φ dt / l.
    pub fn potential_step_field(
        &self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError> {
        let fields = content.fields_mut().ok_or_else(|| {
            HermesError::Config("gravity potential step requires field content".to_string())
        })?;

        let alpha = fields.alpha.as_mut().ok_or_else(|| {
            HermesError::Config("gravity potential step requires alpha field".to_string())
        })?;

        let ell = fields.params.smoothing_length;
        let mass = fields.params.mass_alpha;
        let density_mean = cosmology.density_matter();

        field_potential_step(
            alpha,
            &fields.grid,
            ell,
            mass,
            density_mean,
            scale_factor,
            dt,
        );

        Ok(())
    }

    /// Execute a full KDK step for particle content.
    ///
    /// Convenience method that composes half-kick, drift, recompute, half-kick.
    /// This preserves the exact behavior of the old ParticleMeshDynamics.
    pub fn step_particles_kdk(
        &mut self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor_prev: f64,
        scale_factor_mid: f64,
        scale_factor_next: f64,
    ) -> Result<(), HermesError> {
        // Opening half-kick: a_prev → a_mid
        self.kick_particles(content, cosmology, scale_factor_prev, scale_factor_mid)?;

        // Full drift: a_prev → a_next
        self.drift_particles(content, cosmology, scale_factor_prev, scale_factor_next)?;

        // Recompute forces at new positions and closing half-kick: a_mid → a_next
        self.recompute_and_kick_particles(
            content,
            cosmology,
            scale_factor_next,
            scale_factor_mid,
            scale_factor_next,
        )?;

        Ok(())
    }
}

// ============================================================================
// Internal helpers
// ============================================================================

/// Compute gravitational forces on particles via CIC → Poisson → CIC.
fn compute_particle_forces(
    particles: &crate::physics::particles::Particles,
    solver: &mut PoissonSolver,
    grid: &Grid,
    cosmology: &Cosmology,
    scale_factor: f64,
) -> ParticleForces {
    let density = assign_density(particles, grid);
    let density_mean = cosmology.density_matter();

    let mut overdensity = density;
    overdensity.data /= density_mean;
    overdensity.data -= 1.0;

    let force_field = solver.solve(&overdensity, density_mean, scale_factor);

    interpolate_force(&force_field, particles, grid)
}

/// Apply gravitational potential to an even-subalgebra field.
///
/// Computes density from |α|², solves Poisson for Φ, and rotates
/// the wavefunction phase by -m Φ dt / l.
pub fn field_potential_step(
    alpha: &mut morphis::even_field::EvenField<3>,
    grid: &morphis::grid::Grid<3>,
    ell: f64,
    mass: f64,
    density_mean: f64,
    scale_factor: f64,
    dt: f64,
) {
    let rho = alpha.density(mass);

    let rho_bar_field = Field::scalar_field(grid, metric::euclidean::<3>(), |_| density_mean);
    let poisson_coupling = 4.0 * PI * GRAV * scale_factor * scale_factor;
    let source = &(&rho - &rho_bar_field) * poisson_coupling;
    let phi = source.laplacian_inverse();

    let angle = &phi * (-mass * dt / ell);
    *alpha = alpha.rotate_phase(&angle);
}

// ============================================================================
// Coupling trait implementation
// ============================================================================

impl Coupling for PoissonGravity {
    fn opening_half_step(
        &mut self,
        state: &mut SimulationState,
        cosmology: &Cosmology,
        scale_factor: f64,
        dt_half: f64,
    ) -> Result<(), HermesError> {
        let density_mean = cosmology.density_matter();

        // Apply gravity to each field species.
        let morphis_grid = state.morphis_grid;
        for field in state.fields.values_mut() {
            field_potential_step(
                &mut field.data,
                &morphis_grid,
                field.smoothing_length,
                field.mass,
                density_mean,
                scale_factor,
                dt_half,
            );
        }

        // Particle gravity: kick with cached forces from previous closing half.
        for particles in state.particles.values_mut() {
            let forces = match self.forces_prev.take() {
                Some(f) => f,
                None => compute_particle_forces(
                    particles,
                    &mut self.solver,
                    &self.grid,
                    cosmology,
                    scale_factor,
                ),
            };

            let kick_factor = cosmology.kick_factor(scale_factor, scale_factor);
            crate::physics::integrator::kick(particles, &forces, kick_factor * dt_half);
            self.forces_prev = Some(forces);
        }

        Ok(())
    }

    fn closing_half_step(
        &mut self,
        state: &mut SimulationState,
        cosmology: &Cosmology,
        scale_factor: f64,
        dt_half: f64,
    ) -> Result<(), HermesError> {
        let density_mean = cosmology.density_matter();

        // Apply gravity to each field species (same as opening).
        let morphis_grid = state.morphis_grid;
        for field in state.fields.values_mut() {
            field_potential_step(
                &mut field.data,
                &morphis_grid,
                field.smoothing_length,
                field.mass,
                density_mean,
                scale_factor,
                dt_half,
            );
        }

        // Particle gravity: recompute forces at new positions and kick.
        for particles in state.particles.values_mut() {
            let forces = compute_particle_forces(
                particles,
                &mut self.solver,
                &self.grid,
                cosmology,
                scale_factor,
            );

            let kick_factor = cosmology.kick_factor(scale_factor, scale_factor);
            crate::physics::integrator::kick(particles, &forces, kick_factor * dt_half);
            self.forces_prev = Some(forces);
        }

        Ok(())
    }
}
