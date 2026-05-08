/// Sector: the unit of dynamical content in the engine.
///
/// A sector is a single species — a field or a particle population —
/// that exposes its kinetic flow (T) and potential flow (V) to the
/// integrator. The engine composes these flows across all active
/// sectors through Strang splitting, without knowing what physics
/// any individual sector represents.
///
/// First-degree sectors (Schrodinger-form fields) have:
///   T-flow = Fourier-space phase rotation (dispersive kinetic term)
///   V-flow = real-space phase rotation (gravitational + self-interaction)
///
/// Particle sectors have:
///   T-flow = position drift
///   V-flow = momentum kick from interpolated force
pub mod gross_pitaevskii;
pub mod schrodinger;

use morphis::field::Field;

use crate::engine::state::SimulationState;
use crate::error::HermesError;

// ============================================================================
// Potential: the result of solving the coupling equations
// ============================================================================

/// Gravitational potential produced by the gravity solver.
///
/// Field sectors use `phi` for phase rotation. Particle sectors
/// (when added) will use an interpolated force derived from `phi`.
pub struct Potential {
    /// Gravitational potential on the grid (morphis Field<3>).
    pub phi: Field<3>,
}

// ============================================================================
// Sector trait
// ============================================================================

/// A dynamical sector: one species with kinetic and potential flows.
///
/// The engine calls these methods in a Strang composition:
/// T(dt/2) → V(dt) → T(dt/2) per step, iterating over all sectors
/// at each stage.
pub trait Sector: Send {
    /// The species name this sector governs.
    ///
    /// Must match a key in `SimulationState.fields` (for field sectors)
    /// or `SimulationState.particles` (for particle sectors).
    fn name(&self) -> &str;

    /// Kinetic flow (T-step): advance under the free Lagrangian.
    fn kinetic_flow(
        &mut self,
        state: &mut SimulationState,
        scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError>;

    /// Potential flow (V-step): advance under the gravitational potential.
    fn potential_flow(
        &mut self,
        state: &mut SimulationState,
        potential: &Potential,
        scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError>;

    /// Deposit this sector's mass density onto a shared field.
    ///
    /// Returns the density as a morphis Field<3>. The engine sums
    /// contributions from all sectors before passing the total to
    /// the gravity solver.
    fn deposit_density(&self, state: &SimulationState) -> Result<Field<3>, HermesError>;
}
