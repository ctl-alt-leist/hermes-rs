//! Schrodinger-Poisson dynamics stub.
//!
//! Placeholder dynamics for the field-theoretic dark matter simulation.
//! Currently a no-op — the wavefunction is not evolved. This validates
//! that the Content::Fields path works through the full pipeline.
//! The actual split-step integrator will be implemented here.

use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;

/// Schrodinger-Poisson dynamics for wavefunction dark matter.
///
/// Will implement symmetric split-step: kinetic/2 -> potential -> kinetic/2.
/// Currently a no-op stub.
pub struct SchrodingerPoissonDynamics;

impl SchrodingerPoissonDynamics {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SchrodingerPoissonDynamics {
    fn default() -> Self {
        Self::new()
    }
}

impl Dynamics for SchrodingerPoissonDynamics {
    fn step(
        &mut self,
        _content: &mut Content,
        _cosmology: &Cosmology,
        _scale_factor_prev: f64,
        _scale_factor_next: f64,
    ) -> Result<(), HermesError> {
        // Stub: no evolution. The wavefunction stays at its initial state.
        // TODO: implement split-step kinetic/potential integration.
        Ok(())
    }
}
