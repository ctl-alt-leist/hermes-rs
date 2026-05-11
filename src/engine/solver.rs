/// Gravity solver: Poisson equation from aggregated sector densities.
///
/// Takes the total mass density deposited by all sectors, subtracts
/// the cosmological mean, and solves the Poisson equation for the
/// gravitational potential. Produces both a morphis Field<3> potential
/// (for field sectors) and an optional hermes VectorField force
/// (for particle sectors).
use std::f64::consts::PI;

use morphis::field::Field;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric;

use crate::engine::sector::Potential;
use crate::error::HermesError;
use crate::physics::cic::assign_density;
use crate::physics::constants::G as GRAV;
use crate::physics::field::ScalarField;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;
use crate::physics::poisson::PoissonSolver;

/// FFT-based Poisson solver for gravitational coupling.
///
/// Owns the morphis grid for field potential and optionally the
/// hermes grid + solver for particle force computation.
pub struct GravitySolver {
    /// Morphis grid for spectral operations.
    morphis_grid: MorphisGrid<3>,
    /// Hermes grid and Poisson solver for particle forces.
    /// Present only when the simulation has particle species.
    particle_solver: Option<ParticleSolverState>,
}

/// Hermes-side Poisson solver state for particle force computation.
struct ParticleSolverState {
    grid: Grid,
    solver: PoissonSolver,
}

impl GravitySolver {
    /// Create a gravity solver for field-only simulations.
    pub fn new(grid: MorphisGrid<3>) -> Self {
        Self {
            morphis_grid: grid,
            particle_solver: None,
        }
    }

    /// Create a gravity solver that also computes particle forces.
    pub fn with_particles(morphis_grid: MorphisGrid<3>, hermes_grid: Grid) -> Self {
        let solver = PoissonSolver::new(&hermes_grid);
        Self {
            morphis_grid,
            particle_solver: Some(ParticleSolverState {
                grid: hermes_grid,
                solver,
            }),
        }
    }

    /// Access the hermes grid (for particle drift wrapping).
    pub fn hermes_grid(&self) -> Option<&Grid> {
        self.particle_solver.as_ref().map(|s| &s.grid)
    }

    /// Solve the Poisson equation for the gravitational potential.
    ///
    /// Produces:
    ///   - `phi`: morphis Field<3> potential for field phase rotation
    ///   - `force`: hermes VectorField for particle kicks (if particles present)
    pub fn solve(
        &mut self,
        total_density: &Field<3>,
        density_mean: f64,
        scale_factor: f64,
        particle_species: &[(&str, &Particles)],
    ) -> Result<Potential, HermesError> {
        // Morphis path: field potential via laplacian_inverse.
        let rho_bar = Field::scalar_field(&self.morphis_grid, metric::euclidean::<3>(), |_| {
            density_mean
        });

        let poisson_coupling = 4.0 * PI * GRAV * scale_factor * scale_factor;
        let source = &(total_density - &rho_bar) * poisson_coupling;
        let phi = source.laplacian_inverse();

        // Hermes path: particle force via CIC → PoissonSolver → VectorField.
        let force = if let Some(ref mut ps) = self.particle_solver {
            if !particle_species.is_empty() {
                // Deposit total particle density via CIC.
                // For now, deposit all particle species into one density field.
                let mut total_particle_density: Option<ScalarField> = None;
                for &(_, particles) in particle_species {
                    let rho = assign_density(particles, &ps.grid);
                    total_particle_density = Some(match total_particle_density {
                        None => rho,
                        Some(acc) => {
                            let mut sum = acc;
                            sum.data += &rho.data;
                            sum
                        }
                    });
                }

                if let Some(particle_density) = total_particle_density {
                    let mut overdensity = particle_density;
                    overdensity.data /= density_mean;
                    overdensity.data -= 1.0;
                    Some(ps.solver.solve(&overdensity, density_mean, scale_factor))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };

        Ok(Potential { phi, force })
    }
}
