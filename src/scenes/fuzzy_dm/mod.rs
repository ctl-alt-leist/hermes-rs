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
        cosmology: &Cosmology,
        seed: u64,
    ) -> Result<(Content, Box<dyn Dynamics>), HermesError> {
        let hermes_grid =
            HermesGrid::new(config.simulation.n_cells(), config.simulation.box_length);
        let morphis_grid =
            MorphisGrid::<3>::new(config.simulation.n_cells(), config.simulation.box_length);

        // Calibrated so the quantum Jeans length is ~ L_box / 4.
        // Modes larger than lambda_J collapse under gravity; smaller
        // modes are stabilized by quantum pressure. This gives a few
        // Jeans lengths in the box for visible filament/core formation.
        //
        // lambda_J = 2 pi / (16 pi G rho / (ell/m)^2)^(1/4)
        // For lambda_J ~ 2500 kpc: ell/m ~ 480 kpc^2/Gyr.
        let ell_over_m = 2000.0;
        let mass_alpha = 1e10;

        let params = FieldParams {
            smoothing_length: ell_over_m * mass_alpha,
            mass_alpha,
        };

        let psi = init::random_density_field(
            &hermes_grid,
            cosmology,
            &params,
            config.time.scale_factor_initial,
            seed,
        );

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
