//! Fuzzy dark matter scene: even-subalgebra field in a periodic box.
//!
//! A single dark matter field α evolved by the Schrodinger-Poisson
//! system. The lightest field-theoretic scene — validates the Content::Fields
//! path through the full pipeline.

pub mod init;

use morphis::grid::Grid as MorphisGrid;

use crate::config::Configuration;
use crate::core::content::{Content, FieldParams, FieldState};
use crate::core::dynamics::Dynamics;
use crate::core::schrodinger_dynamics::SchrodingerPoissonDynamics;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid as HermesGrid;
use crate::scenes::Scene;

const SCENE_DEFAULTS: &str = include_str!("defaults.toml");

/// Fuzzy dark matter: wavefunction + self-gravity.
pub struct FuzzyDM;

impl Scene for FuzzyDM {
    fn name(&self) -> &str {
        "fuzzy-dm"
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

        let alpha = init::random_density_field(
            &hermes_grid,
            cosmology,
            &params,
            config.time.scale_factor_initial(),
            config.initialization.perturbation_amplitude,
            config.initialization.band_pass,
            seed,
        );

        let field_state = FieldState {
            grid: morphis_grid,
            alpha: Some(alpha),
            beta: None,
            gamma: None,
            params,
        };

        let dynamics = SchrodingerPoissonDynamics::new();

        Ok((Content::Fields(field_state), Box::new(dynamics)))
    }
}
