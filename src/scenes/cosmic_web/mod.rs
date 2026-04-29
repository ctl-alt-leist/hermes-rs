//! Cosmic web scene: Zel'dovich PM simulation in a periodic box.

pub mod init;

use crate::config::Configuration;
use crate::error::HermesError;
use crate::physics::content::Content;
use crate::physics::cosmology::Cosmology;
use crate::physics::dynamics::Dynamics;
use crate::physics::grid::Grid;
use crate::physics::pm_dynamics::ParticleMeshDynamics;
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

    fn initialize(
        &self,
        config: &Configuration,
        cosmology: &Cosmology,
        seed: u64,
    ) -> Result<(Content, Box<dyn Dynamics>), HermesError> {
        let grid = Grid::new(config.simulation.n_cells, config.simulation.box_length);

        let particles = init::zeldovich_init(
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
