/// Physics engine: composable dynamics from free evolution and couplings.
///
/// The engine holds the simulation state and a set of physics modules
/// (free evolution per species, couplings between species). A timestep
/// is a Strang splitting composition:
///
///   1. Coupling half-steps (dt/2)
///   2. Free evolution full steps (dt)
///   3. Coupling half-steps (dt/2)
pub mod coupling;
pub mod free;
pub mod state;

use std::collections::BTreeMap;

use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;

use coupling::Coupling;
use free::FreeEvolution;
use state::SimulationState;

/// The composable physics engine.
///
/// Owns the simulation state and the physics modules. Advances the
/// state forward by composing free evolution and coupling steps.
pub struct Engine {
    /// The simulation state: species on a grid at a moment in time.
    pub state: SimulationState,
    /// Free evolution modules, keyed by field species name.
    pub free_modules: BTreeMap<String, Box<dyn FreeEvolution>>,
    /// Coupling modules (gravity, electromagnetic, etc.).
    pub couplings: Vec<Box<dyn Coupling>>,
    /// Cosmological background (or None for static spacetime).
    pub cosmology: Option<Cosmology>,
}

impl Engine {
    /// Execute one full timestep via Strang splitting.
    ///
    /// For FLRW spacetimes, dt is computed from the scale factor step
    /// and the Hubble parameter. For static spacetimes, dt is passed
    /// directly from the time schedule.
    pub fn step(&mut self, scale_factor: f64, dt: f64) -> Result<(), HermesError> {
        let cosmology = self.cosmology.clone();

        // 1. Opening coupling half-steps (dt/2).
        for coupling in &mut self.couplings {
            if let Some(ref cosmo) = cosmology {
                coupling.opening_half_step(&mut self.state, cosmo, scale_factor, dt / 2.0)?;
            }
        }

        // 2. Free evolution full steps (dt) for each field species.
        for (name, module) in &mut self.free_modules {
            if let Some(field) = self.state.fields.get_mut(name) {
                let grid = self.state.morphis_grid;
                module.step(field, &grid, scale_factor, dt)?;
            }
        }

        // 3. Closing coupling half-steps (dt/2).
        for coupling in &mut self.couplings {
            if let Some(ref cosmo) = cosmology {
                coupling.closing_half_step(&mut self.state, cosmo, scale_factor, dt / 2.0)?;
            }
        }

        self.state.step += 1;

        Ok(())
    }
}
