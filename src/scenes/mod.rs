//! Simulation scenes.
//!
//! A `Scene` defines how to initialize a simulation — what content to
//! create (particles, fields, or both), what dynamics to attach, and
//! what config defaults to apply. The scene returns both the initial
//! content and the dynamics module that evolves it.

pub mod cosmic_web;
pub mod galaxy_group;

use crate::config::Configuration;
use crate::error::HermesError;
use crate::physics::content::Content;
use crate::physics::cosmology::Cosmology;
use crate::physics::dynamics::Dynamics;

/// Trait for simulation scenarios.
///
/// Each scene provides:
/// - A name for CLI selection
/// - Optional default config overrides
/// - Config validation for scene-specific constraints
/// - Content and dynamics initialization
pub trait Scene {
    /// Human-readable name (used in CLI `--scene` flag).
    fn name(&self) -> &str;

    /// Scene-specific config defaults.
    fn default_overrides(&self) -> Option<toml::Value> {
        None
    }

    /// Validate that the final config is compatible with this scene.
    fn validate(&self, _config: &Configuration) -> Result<(), HermesError> {
        Ok(())
    }

    /// Initialize content and dynamics for this scene.
    fn initialize(
        &self,
        config: &Configuration,
        cosmology: &Cosmology,
        seed: u64,
    ) -> Result<(Content, Box<dyn Dynamics>), HermesError>;
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
