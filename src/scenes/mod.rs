//! Simulation scenarios.
//!
//! A `Scene` defines how to initialize a simulation — what initial conditions
//! to generate, what physics to activate, and what default configuration
//! overrides to apply. The runner looks up a scene by name and delegates
//! initialization to it.

mod cosmic_web;

pub use cosmic_web::CosmicWeb;

use crate::config::Configuration;
use crate::error::HermesError;
use crate::physics::simulation::Simulation;

/// Trait for simulation scenarios.
///
/// Each scene knows how to initialize a `Simulation` from a configuration.
/// Different scenes produce different initial conditions (Zel'dovich cosmic
/// web, isolated halo collapse, Zeldovich pancake test, etc.).
pub trait Scene {
    /// Human-readable name of the scene.
    fn name(&self) -> &str;

    /// Initialize a simulation for this scene.
    fn initialize(&self, config: &Configuration, seed: u64) -> Result<Simulation, HermesError>;
}

/// Look up a scene by name.
pub fn scene_by_name(name: &str) -> Result<Box<dyn Scene>, HermesError> {
    match name {
        "cosmic-web" => Ok(Box::new(CosmicWeb)),
        _ => Err(HermesError::Config(format!("unknown scene: {name}"))),
    }
}

/// List available scene names.
pub fn available_scenes() -> Vec<&'static str> {
    vec!["cosmic-web"]
}
