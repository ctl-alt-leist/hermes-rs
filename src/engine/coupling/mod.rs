/// Coupling modules: cross-species interactions.
///
/// Each coupling knows how to apply an interaction between species.
/// The engine composes couplings with free evolution in a Strang
/// splitting loop.
pub mod poisson;

use crate::engine::state::SimulationState;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;

/// A coupling module that applies an interaction across species.
///
/// The engine calls `half_step` before and after the free evolution
/// steps to implement Strang splitting.
pub trait Coupling: Send {
    /// Apply a half-step of the coupling interaction.
    fn half_step(
        &mut self,
        state: &mut SimulationState,
        cosmology: &Cosmology,
        scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError>;
}
