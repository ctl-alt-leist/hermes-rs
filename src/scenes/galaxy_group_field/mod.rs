//! Galaxy group field scene: colliding NFW halos as a Schrodinger-Poisson field.
//!
//! The same halo geometry as galaxy-group-pm, but the dark matter is
//! represented as a wavefunction α in the even subalgebra. The NFW
//! density profiles and infall velocities are encoded via the inverse
//! Madelung transform.

pub mod init;

use morphis::grid::Grid as MorphisGrid;

use crate::config::Configuration;
use crate::core::content::{Content, FieldParams, FieldState};
use crate::core::dynamics::Dynamics;
use crate::core::schrodinger_dynamics::SchrodingerPoissonDynamics;
use crate::engine::coupling::poisson::PoissonGravity;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid as HermesGrid;
use crate::scenes::Scene;

const SCENE_DEFAULTS: &str = include_str!("defaults.toml");

/// Galaxy group formation via Schrodinger-Poisson field theory.
pub struct GalaxyGroupField;

impl Scene for GalaxyGroupField {
    fn name(&self) -> &str {
        "galaxy-group-field"
    }

    fn default_overrides(&self) -> Option<toml::Value> {
        toml::from_str(SCENE_DEFAULTS).ok()
    }

    fn validate(&self, config: &Configuration) -> Result<(), HermesError> {
        if config.simulation.box_length > 10_000.0 {
            return Err(HermesError::Config(
                "galaxy-group-field scene expects box_length <= 10 Mpc (10000 kpc)".to_string(),
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
        let hermes_grid =
            HermesGrid::new(config.simulation.n_cells(), config.simulation.box_length);
        let morphis_grid =
            MorphisGrid::<3>::new(config.simulation.n_cells(), config.simulation.box_length);

        let ell_over_m = config.field.length_scale;
        let mass_alpha = config.field.mass;

        let params = FieldParams {
            smoothing_length: ell_over_m * mass_alpha,
            mass_alpha,
        };

        let alpha = init::colliding_halos_field(
            &hermes_grid,
            cosmology,
            &params,
            config.time.scale_factor_initial(),
            seed,
        );

        let field_state = FieldState {
            grid: morphis_grid,
            alpha: Some(alpha),
            beta: None,
            gamma: None,
            params,
        };

        let gravity = PoissonGravity::new(hermes_grid);
        let dynamics = SchrodingerPoissonDynamics::new(gravity);

        Ok((Content::Fields(field_state), Box::new(dynamics)))
    }
}
