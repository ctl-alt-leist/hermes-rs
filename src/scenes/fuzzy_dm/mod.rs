//! Fuzzy dark matter scene: wavefunction in a periodic box.
//!
//! A single even-subalgebra wavefunction evolved by the Schrodinger-Poisson
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
        _cosmology: &Cosmology,
        _seed: u64,
    ) -> Result<(Content, Box<dyn Dynamics>), HermesError> {
        let hermes_grid = HermesGrid::new(config.simulation.n_cells, config.simulation.box_length);
        let morphis_grid =
            MorphisGrid::<3>::new(config.simulation.n_cells, config.simulation.box_length);

        // Parameters tuned for gravity-dispersion balance at the box
        // fundamental mode. ell/m ~ 480 gives dispersive pressure
        // comparable to gravitational collapse, with CFL ~ 0.65
        // (stable for 300 steps over 20 Gyr).
        let params = FieldParams {
            smoothing_length: 480.0,
            mass_alpha: 1.0,
        };

        let psi = init::gaussian_wavepacket(&hermes_grid, &params);

        let field_state = FieldState {
            grid: morphis_grid,
            psi: Some(psi),
            beta: None,
            gamma: None,
            params,
        };

        let dynamics = SchrodingerPoissonDynamics::new();

        Ok((Content::Fields(field_state), Box::new(dynamics)))
    }
}
