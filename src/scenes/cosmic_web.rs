//! Cosmic web scene: Zel'dovich PM simulation in a periodic box.
//!
//! Dark-matter-only particle-mesh simulation with Zel'dovich initial
//! conditions from a linear power spectrum. This is the default scene
//! and the foundation onto which finer scales, baryonic matter, and
//! learned closures are eventually attached.

use crate::config::Configuration;
use crate::error::HermesError;
use crate::physics::simulation::Simulation;
use crate::scenes::Scene;

/// Zel'dovich PM simulation in a periodic cosmological box.
pub struct CosmicWeb;

impl Scene for CosmicWeb {
    fn name(&self) -> &str {
        "cosmic-web"
    }

    fn initialize(&self, config: &Configuration, seed: u64) -> Result<Simulation, HermesError> {
        Simulation::from_config(config.clone(), seed)
    }
}
