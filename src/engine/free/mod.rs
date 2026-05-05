/// Free evolution modules: per-species dynamics from free Lagrangian terms.
///
/// Each module knows how to advance a single field species by one timestep
/// under its free (non-interacting) dynamics. Particles have no free
/// evolution — they are inertial between coupling kicks.
pub mod schrodinger;

use crate::engine::state::FieldEntry;
use crate::error::HermesError;

/// A free evolution module for a single field species.
///
/// Advances the field by dt under the free Lagrangian (no couplings).
/// The engine calls this between coupling half-steps in the Strang loop.
pub trait FreeEvolution: Send {
    /// Advance the field by dt at the given scale factor.
    fn step(
        &mut self,
        field: &mut FieldEntry,
        grid: &morphis::grid::Grid<3>,
        scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError>;
}
