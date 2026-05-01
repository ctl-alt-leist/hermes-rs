//! Dynamics trait: how content evolves in time.
//!
//! A dynamics module knows how to advance a simulation's content by
//! one timestep. Different content kinds require different dynamics:
//! particles use KDK leapfrog, fields use split-step, mixed content
//! couples both through a shared Poisson solve.

use crate::core::content::Content;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;

/// A dynamics module that advances content by one step.
///
/// The simulation driver calls `step()` without knowing which
/// dynamics implementation is running. The dynamics module owns
/// whatever internal state it needs (FFT plans, solver workspace, etc.)
/// and mutates the content in place.
pub trait Dynamics: Send {
    /// Advance the content from scale_factor_prev to scale_factor_next.
    fn step(
        &mut self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor_prev: f64,
        scale_factor_next: f64,
    ) -> Result<(), HermesError>;
}
