//! Cosmic web scene: Zel'dovich PM simulation in a periodic box.
//!
//! Dark-matter-only particle-mesh simulation with Zel'dovich initial
//! conditions from a linear power spectrum. This is the default scene
//! and the foundation onto which finer scales, baryonic matter, and
//! learned closures are eventually attached.

pub mod init;

use crate::config::Configuration;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;
use crate::scenes::Scene;

const SCENE_DEFAULTS: &str = include_str!("defaults.toml");

/// Zel'dovich PM simulation in a periodic cosmological box.
pub struct CosmicWeb;

impl Scene for CosmicWeb {
    fn name(&self) -> &str {
        "cosmic-web"
    }

    fn default_overrides(&self) -> Option<toml::Value> {
        toml::from_str(SCENE_DEFAULTS).ok()
    }

    fn initialize_particles(
        &self,
        grid: &Grid,
        cosmology: &Cosmology,
        config: &Configuration,
        seed: u64,
    ) -> Result<Particles, HermesError> {
        init::zeldovich_init(
            config.simulation.n_particles,
            grid,
            cosmology,
            config.time.scale_factor_initial,
            seed,
        )
    }
}
