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

        // Calibrated so de Broglie wavelength at sigma_v ~ 300 km/s
        // matches the cell size: ell/m = sigma_v * dx / (2 pi).
        // sigma_v = 307 kpc/Gyr, dx = box_length / n_cells.
        // Calibrated so de Broglie wavelength at the characteristic velocity
        // matches the cell size: ell/m = sigma_v * dx / (2 pi).
        let dx = config.simulation.box_length / config.simulation.n_cells as f64;
        let sigma_v = 307.0; // 300 km/s in kpc/Gyr
        let ell_over_m = sigma_v * dx / (2.0 * std::f64::consts::PI);
        let mass_alpha = 1e10;

        let params = FieldParams {
            smoothing_length: ell_over_m * mass_alpha,
            mass_alpha,
        };

        let psi = init::nonlinear_modes(&hermes_grid, &params, _cosmology);

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
