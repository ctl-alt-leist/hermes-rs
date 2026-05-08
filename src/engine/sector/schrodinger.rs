/// Schrodinger sector: first-degree field with dispersive kinetics.
///
/// Advances an even-subalgebra field under:
///   T-flow: Fourier-space phase rotation (free kinetic term)
///   V-flow: real-space phase rotation (gravitational potential)
///
/// The kinetic step delegates to `kinetic_step()` in the core module.
/// The potential step applies the gravitational potential produced by
/// the engine's gravity solver.
use crate::engine::sector::{Potential, Sector};
use crate::engine::state::SimulationState;
use crate::error::HermesError;

/// Schrodinger sector for a named field species.
///
/// Stateless — the field parameters (smoothing_length, mass) are read
/// from the FieldEntry in SimulationState at each call.
pub struct SchrodingerSector {
    name: String,
}

impl SchrodingerSector {
    /// Create a Schrodinger sector for the named field species.
    pub fn new(name: String) -> Self {
        Self { name }
    }
}

impl Sector for SchrodingerSector {
    fn name(&self) -> &str {
        &self.name
    }

    fn kinetic_flow(
        &mut self,
        state: &mut SimulationState,
        scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError> {
        let field = state.fields.get_mut(&self.name).ok_or_else(|| {
            HermesError::Config(format!("field '{}' not found in state", self.name))
        })?;

        let grid = state.morphis_grid;

        crate::core::schrodinger_dynamics::kinetic_step(
            &mut field.data,
            &grid,
            field.smoothing_length,
            field.mass,
            scale_factor,
            dt,
        );

        Ok(())
    }

    fn potential_flow(
        &mut self,
        state: &mut SimulationState,
        potential: &Potential,
        _scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError> {
        let field = state.fields.get_mut(&self.name).ok_or_else(|| {
            HermesError::Config(format!("field '{}' not found in state", self.name))
        })?;

        let angle = &potential.phi * (-field.mass * dt / field.smoothing_length);
        field.data = field.data.rotate_phase(&angle);

        Ok(())
    }

    fn deposit_density(
        &self,
        state: &SimulationState,
    ) -> Result<morphis::field::Field<3>, HermesError> {
        let field = state.fields.get(&self.name).ok_or_else(|| {
            HermesError::Config(format!("field '{}' not found in state", self.name))
        })?;

        Ok(field.data.density(field.mass))
    }
}
