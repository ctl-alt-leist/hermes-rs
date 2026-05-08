/// Physics engine: composable sector dynamics via Strang splitting.
///
/// The engine holds the simulation state and a collection of sectors,
/// each exposing a kinetic flow (T) and potential flow (V). A timestep
/// is a merged Strang splitting composition that saves one kinetic FFT
/// pair per field sector on interior steps:
///
///   First step:    T(dt/2) → V(dt)
///   Interior:      T(dt)   → V(dt)       (merged half-steps)
///   Finalize:      T(dt/2)                (closing bookend)
///
/// This is equivalent to the symmetric form T/2 → V → T/2 at every
/// step, but avoids redundant kinetic evaluations where adjacent
/// half-steps merge.
pub mod coupling;
pub mod free;
pub mod sector;
pub mod solver;
pub mod state;

use morphis::field::Field;

use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::particles::Particles;

use sector::Sector;
use solver::GravitySolver;
use state::SimulationState;

/// The composable physics engine.
///
/// Owns the simulation state and the physics sectors. Advances the
/// state forward by composing kinetic and potential flows via merged
/// Strang splitting.
pub struct Engine {
    /// The simulation state: species on a grid at a moment in time.
    pub state: SimulationState,
    /// Dynamical sectors (one per species).
    pub sectors: Vec<Box<dyn Sector>>,
    /// Gravity solver (None if gravity is disabled).
    pub solver: Option<GravitySolver>,
    /// Cosmological background (or None for static spacetime).
    pub cosmology: Option<Cosmology>,
    /// Whether the next step needs a half-width opening kinetic step.
    /// True at the start of a run; false after the first step, when
    /// the trailing half from the previous step merges with the
    /// leading half of the current step into a full kinetic step.
    needs_opening_half: bool,
}

impl Engine {
    /// Create a new engine.
    pub fn new(
        state: SimulationState,
        sectors: Vec<Box<dyn Sector>>,
        solver: Option<GravitySolver>,
        cosmology: Option<Cosmology>,
    ) -> Self {
        Self {
            state,
            sectors,
            solver,
            cosmology,
            needs_opening_half: true,
        }
    }

    /// Execute one timestep via merged Strang splitting.
    ///
    /// The first call applies a half-width kinetic step (opening
    /// bookend). Subsequent calls apply a full-width kinetic step
    /// (merged adjacent halves). Every call applies a full-width
    /// potential step. Call `finalize()` after the last step to
    /// apply the closing bookend.
    pub fn step(&mut self, scale_factor: f64, dt: f64) -> Result<(), HermesError> {
        let Engine {
            state,
            sectors,
            solver,
            cosmology,
            needs_opening_half,
        } = self;

        // 1. Kinetic step: T(dt/2) on first call, T(dt) on interior.
        let dt_kinetic = if *needs_opening_half { dt / 2.0 } else { dt };

        for sector in sectors.iter_mut() {
            sector.kinetic_flow(state, scale_factor, dt_kinetic)?;
        }

        *needs_opening_half = false;

        // 2. Potential full step V(dt): deposit → solve → apply.
        if let Some(solver) = solver {
            let cosmology = cosmology
                .as_ref()
                .ok_or_else(|| HermesError::Config("gravity requires cosmology".to_string()))?;
            let density_mean = cosmology.density_matter();

            let total_density = aggregate_density(sectors, state)?;

            // Collect particle species references for the solver.
            let particle_species: Vec<(&str, &Particles)> = state
                .particles
                .iter()
                .map(|(name, p)| (name.as_str(), p))
                .collect();

            let potential = solver.solve(
                &total_density,
                density_mean,
                scale_factor,
                &particle_species,
            )?;

            for sector in sectors.iter_mut() {
                sector.potential_flow(state, &potential, scale_factor, dt)?;
            }
        }

        state.step += 1;

        Ok(())
    }

    /// Apply the closing kinetic half-step after the last timestep.
    ///
    /// Call this once after the stepping loop completes to close the
    /// Strang composition. The state is then at a fully symmetric
    /// composition point, suitable for final snapshots and diagnostics.
    pub fn finalize(&mut self, scale_factor: f64, dt: f64) -> Result<(), HermesError> {
        let dt_half = dt / 2.0;

        for sector in self.sectors.iter_mut() {
            sector.kinetic_flow(&mut self.state, scale_factor, dt_half)?;
        }

        Ok(())
    }
}

// ============================================================================
// Internal helpers
// ============================================================================

/// Sum density contributions from all sectors.
fn aggregate_density(
    sectors: &[Box<dyn Sector>],
    state: &SimulationState,
) -> Result<Field<3>, HermesError> {
    let mut total: Option<Field<3>> = None;

    for sector in sectors {
        let rho = sector.deposit_density(state)?;
        total = Some(match total {
            None => rho,
            Some(acc) => &acc + &rho,
        });
    }

    total.ok_or_else(|| HermesError::Config("no sectors to deposit density from".to_string()))
}
