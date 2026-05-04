/// Coupling modules: cross-species interactions.
///
/// Each coupling knows how to apply an interaction between species.
/// The engine composes couplings with free evolution in a Strang
/// splitting loop: opening_half → free → closing_half.
pub mod poisson;

use crate::engine::state::SimulationState;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;

/// A coupling module that applies an interaction across species.
///
/// The engine calls `opening_half_step` before free evolution and
/// `closing_half_step` after. Both receive dt/2. The distinction
/// matters for particles: the opening half can reuse cached forces
/// from the previous closing half, while the closing half must
/// recompute forces at the new (post-drift) positions.
///
/// For fields, both halves are identical: Poisson solve from
/// current density, phase rotation by dt/2.
pub trait Coupling: Send {
    /// Apply the opening coupling half-step (before free evolution).
    fn opening_half_step(
        &mut self,
        state: &mut SimulationState,
        cosmology: &Cosmology,
        scale_factor: f64,
        dt_half: f64,
    ) -> Result<(), HermesError>;

    /// Apply the closing coupling half-step (after free evolution).
    fn closing_half_step(
        &mut self,
        state: &mut SimulationState,
        cosmology: &Cosmology,
        scale_factor: f64,
        dt_half: f64,
    ) -> Result<(), HermesError>;
}
