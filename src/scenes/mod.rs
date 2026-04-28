//! Simulation scenes.
//!
//! A `Scene` defines how to initialize a simulation for a particular
//! physical setup — what initial conditions to generate, what config
//! defaults to apply, and what validation to enforce. The physics engine
//! (Poisson, integrator, CIC) is shared across all scenes.

pub mod cosmic_web;
pub mod galaxy_group;

use crate::config::Configuration;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// Trait for simulation scenarios.
///
/// Each scene provides:
/// - A name for CLI selection
/// - Optional default config overrides (merged above global defaults)
/// - Config validation for scene-specific constraints
/// - Particle initialization for the scene's physics
pub trait Scene {
    /// Human-readable name (used in CLI `--scene` flag).
    fn name(&self) -> &str;

    /// Scene-specific config defaults, merged on top of global defaults
    /// but below user config file and CLI overrides.
    fn default_overrides(&self) -> Option<toml::Value> {
        None
    }

    /// Validate that the final config is compatible with this scene.
    fn validate(&self, _config: &Configuration) -> Result<(), HermesError> {
        Ok(())
    }

    /// Initialize particles for this scene.
    fn initialize_particles(
        &self,
        grid: &Grid,
        cosmology: &Cosmology,
        config: &Configuration,
        seed: u64,
    ) -> Result<Particles, HermesError>;
}

/// Look up a scene by name.
pub fn scene_by_name(name: &str) -> Result<Box<dyn Scene>, HermesError> {
    match name {
        "cosmic-web" => Ok(Box::new(cosmic_web::CosmicWeb)),
        "galaxy-group" => Ok(Box::new(galaxy_group::GalaxyGroup)),
        _ => Err(HermesError::Config(format!(
            "unknown scene: {name}. available: {}",
            available_scenes().join(", ")
        ))),
    }
}

/// List available scene names.
pub fn available_scenes() -> Vec<&'static str> {
    vec!["cosmic-web", "galaxy-group"]
}
