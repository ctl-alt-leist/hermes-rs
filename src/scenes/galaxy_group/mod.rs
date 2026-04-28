//! Galaxy group scene: constrained Zel'dovich in a smaller periodic box.
//!
//! A ~3 Mpc periodic box focused on the formation of a galaxy group.
//! Uses the same Zel'dovich approximation as the cosmic web, but with
//! a smaller box, later start redshift, and an initial long-wavelength
//! overdensity bias that ensures the region collapses into a group-mass
//! structure by z = 0.

pub mod init;

use crate::config::Configuration;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;
use crate::scenes::Scene;

const SCENE_DEFAULTS: &str = include_str!("defaults.toml");

/// Galaxy group formation in a smaller periodic box.
pub struct GalaxyGroup;

impl Scene for GalaxyGroup {
    fn name(&self) -> &str {
        "galaxy-group"
    }

    fn default_overrides(&self) -> Option<toml::Value> {
        toml::from_str(SCENE_DEFAULTS).ok()
    }

    fn validate(&self, config: &Configuration) -> Result<(), HermesError> {
        if config.simulation.box_length > 10_000.0 {
            return Err(HermesError::Config(
                "galaxy-group scene expects box_length <= 10 Mpc (10000 kpc)".to_string(),
            ));
        }

        Ok(())
    }

    fn initialize_particles(
        &self,
        grid: &Grid,
        cosmology: &Cosmology,
        config: &Configuration,
        seed: u64,
    ) -> Result<Particles, HermesError> {
        init::colliding_halos_init(
            config.simulation.n_particles,
            grid,
            cosmology,
            config.time.scale_factor_initial,
            seed,
        )
    }
}
