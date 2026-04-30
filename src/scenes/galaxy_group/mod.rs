//! Galaxy group scene: multiple colliding NFW halos.

pub mod init;

use crate::config::Configuration;
use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::core::pm_dynamics::ParticleMeshDynamics;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::scenes::Scene;

const SCENE_DEFAULTS: &str = include_str!("defaults.toml");

/// Galaxy group formation: colliding NFW halos.
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

    fn initialize(
        &self,
        config: &Configuration,
        cosmology: &Cosmology,
        seed: u64,
    ) -> Result<(Content, Box<dyn Dynamics>), HermesError> {
        let grid = Grid::new(config.simulation.n_cells(), config.simulation.box_length);

        let particles = init::colliding_halos_init(
            config.simulation.n_particles,
            &grid,
            cosmology,
            config.time.scale_factor_initial,
            seed,
        )?;

        let dynamics = ParticleMeshDynamics::new(grid);

        Ok((Content::Particles(particles), Box::new(dynamics)))
    }
}
