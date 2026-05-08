/// Gross-Pitaevskii sector: Schrodinger kinetics with self-interaction.
///
/// Advances an even-subalgebra field under:
///   T-flow: Fourier-space phase rotation (free kinetic term, same as Schrodinger)
///   V-flow: real-space phase rotation from gravitational potential + self-interaction
///
/// The self-interaction enters the potential flow as an additional
/// phase angle proportional to the local density:
///   angle_total = -(m * phi + g * |alpha|^2) * dt / ell
///
/// This gives baryons an effective pressure (sound speed) that
/// prevents collapse below the de Broglie scale.
use crate::engine::sector::{Potential, Sector};
use crate::engine::state::SimulationState;
use crate::error::HermesError;

/// Gross-Pitaevskii sector for a named field species.
///
/// Stateless — the field parameters (smoothing_length, mass,
/// self_interaction) are read from the FieldEntry at each call.
pub struct GrossPitaevskiiSector {
    name: String,
}

impl GrossPitaevskiiSector {
    /// Create a Gross-Pitaevskii sector for the named field species.
    pub fn new(name: String) -> Self {
        Self { name }
    }
}

impl Sector for GrossPitaevskiiSector {
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

        let mass = field.mass;
        let ell = field.smoothing_length;
        let coupling = field.self_interaction.unwrap_or(0.0);

        // Gravitational phase angle: -m * phi * dt / ell
        let mut angle = &potential.phi * (-mass * dt / ell);

        // Self-interaction phase angle: -g * |alpha|^2 * dt / ell
        if coupling != 0.0 {
            let density_sq = field.data.norm_squared();
            let angle_self = &density_sq * (-coupling * dt / ell);
            angle = &angle + &angle_self;
        }

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
