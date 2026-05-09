//! Initial condition generators.
//!
//! Each submodule provides initialization functions for a specific
//! physical scenario. The initialization method is selected by the
//! `method` field in `[simulation.initialization]`:
//!
//!   - `"zeldovich"` — Zel'dovich displacement from CDM power spectrum
//!   - `"nfw-group"` — colliding NFW dark matter halos
//!
//! The top-level `initialize_from_config()` dispatches to the correct
//! init function based on the EngineConfig.

pub mod nfw;
pub mod nfw_field;
pub mod zeldovich;
pub mod zeldovich_field;

// Re-export the most-used items.
pub use zeldovich::{power_spectrum, transfer_function, zeldovich_init};

use crate::config::EngineConfig;
use crate::core::content::{Content, FieldParams, FieldState};
use crate::core::dynamics::Dynamics;
use crate::core::mixed_dynamics::MixedDynamics;
use crate::core::pm_dynamics::ParticleMeshDynamics;
use crate::core::schrodinger_dynamics::SchrodingerPoissonDynamics;
use crate::engine::coupling::poisson::PoissonGravity;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;

/// Initialize content and dynamics from an EngineConfig.
///
/// Dispatches on `config.simulation.initialization.method` and the
/// ontology (particles vs fields) to select the right initialization
/// function and dynamics module.
pub fn initialize_from_config(
    config: &EngineConfig,
    cosmology: &Cosmology,
) -> Result<(Content, Box<dyn Dynamics>), HermesError> {
    let method = config.simulation.initialization.method.as_str();
    let has_particles = config.ontology.has_particles();
    let has_fields = config.ontology.has_fields();
    let seed = config.simulation.initialization.seed;

    let n_cells = config.simulation.grid.n_cells;
    let box_length = config.simulation.grid.box_length;
    let grid = Grid::new(n_cells, box_length);

    let scale_factor_initial = config
        .simulation
        .time
        .scale_factor_range
        .map(|r| r[0])
        .unwrap_or(1.0);

    match (method, has_particles, has_fields) {
        ("zeldovich", true, false) => {
            // Particle Zel'dovich from CDM power spectrum.
            let n_per_side = config
                .ontology
                .particles
                .values()
                .next()
                .map(|p| p.n)
                .unwrap_or(32);

            let particles = zeldovich::zeldovich_init(
                n_per_side,
                &grid,
                cosmology,
                scale_factor_initial,
                seed,
            )?;
            let dynamics = ParticleMeshDynamics::new(grid);

            Ok((Content::Particles(particles), Box::new(dynamics)))
        }

        ("zeldovich", false, true) => {
            // Field Zel'dovich wavefunction(s).
            let (ell, mass) = field_params_from_config(config)?;
            let morphis_grid = morphis::grid::Grid::<3>::new(n_cells, box_length);

            let params = FieldParams {
                smoothing_length: ell,
                mass_alpha: mass,
            };

            let spectrum = config.simulation.initialization.spectrum.as_str();
            let perturbation_amplitude = config.simulation.initialization.perturbation_amplitude;

            let n_fields = config.ontology.fields.len();
            let density_fraction = 1.0 / n_fields as f64;

            let init_field = |frac: f64, field_seed: u64| -> Result<_, HermesError> {
                if spectrum == "random" {
                    let band_pass = config.simulation.initialization.band_pass;
                    Ok(zeldovich_field::random_density_field(
                        &grid,
                        cosmology,
                        &params,
                        scale_factor_initial,
                        perturbation_amplitude,
                        band_pass,
                        field_seed,
                        frac,
                    ))
                } else {
                    zeldovich_field::zeldovich_wavefunction(
                        &grid,
                        cosmology,
                        &params,
                        scale_factor_initial,
                        perturbation_amplitude,
                        field_seed,
                        frac,
                    )
                }
            };

            // Initialize alpha (first field species).
            let alpha = init_field(density_fraction, seed)?;

            // Initialize beta (second field species, if present) with a different seed.
            let beta = if n_fields >= 2 {
                Some(init_field(density_fraction, seed + 1)?)
            } else {
                None
            };

            let field_state = FieldState {
                grid: morphis_grid,
                alpha: Some(alpha),
                beta,
                gamma: None,
                params,
            };

            let gravity = PoissonGravity::new(grid);
            let dynamics = SchrodingerPoissonDynamics::new(gravity);

            Ok((Content::Fields(field_state), Box::new(dynamics)))
        }

        ("nfw-group", true, false) => {
            // NFW halo particles.
            let n_per_side = config
                .ontology
                .particles
                .values()
                .next()
                .map(|p| p.n)
                .unwrap_or(32);

            let halos = halo_configs_from_config(config);
            let particles = nfw::colliding_halos_init(
                n_per_side,
                &grid,
                cosmology,
                scale_factor_initial,
                seed,
                &halos,
                1.0,
            )?;
            let dynamics = ParticleMeshDynamics::new(grid);

            Ok((Content::Particles(particles), Box::new(dynamics)))
        }

        ("nfw-group", false, true) => {
            // NFW halo field.
            let (ell, mass) = field_params_from_config(config)?;
            let morphis_grid = morphis::grid::Grid::<3>::new(n_cells, box_length);

            let params = FieldParams {
                smoothing_length: ell,
                mass_alpha: mass,
            };

            let halos = halo_configs_from_config(config);
            let alpha = nfw_field::colliding_halos_field(
                &grid,
                cosmology,
                &params,
                scale_factor_initial,
                seed,
                &halos,
                1.0,
            );

            let field_state = FieldState {
                grid: morphis_grid,
                alpha: Some(alpha),
                beta: None,
                gamma: None,
                params,
            };

            let gravity = PoissonGravity::new(grid);
            let dynamics = SchrodingerPoissonDynamics::new(gravity);

            Ok((Content::Fields(field_state), Box::new(dynamics)))
        }

        ("nfw-group", true, true) => {
            // Mixed: particles and fields each carry half the cosmological density.
            let n_per_side = config
                .ontology
                .particles
                .values()
                .next()
                .map(|p| p.n)
                .unwrap_or(32);

            let (ell, mass) = field_params_from_config(config)?;
            let morphis_grid = morphis::grid::Grid::<3>::new(n_cells, box_length);

            let params = FieldParams {
                smoothing_length: ell,
                mass_alpha: mass,
            };

            let halos = halo_configs_from_config(config);

            let particles = nfw::colliding_halos_init(
                n_per_side,
                &grid,
                cosmology,
                scale_factor_initial,
                seed,
                &halos,
                0.5,
            )?;

            let alpha = nfw_field::colliding_halos_field(
                &grid,
                cosmology,
                &params,
                scale_factor_initial,
                seed,
                &halos,
                0.5,
            );

            let field_state = FieldState {
                grid: morphis_grid,
                alpha: Some(alpha),
                beta: None,
                gamma: None,
                params,
            };

            let gravity = PoissonGravity::new(grid);
            let dynamics = MixedDynamics::new(gravity);

            Ok((
                Content::Mixed {
                    particles,
                    fields: field_state,
                },
                Box::new(dynamics),
            ))
        }

        _ => Err(HermesError::Config(format!(
            "unsupported initialization: method={method}, particles={has_particles}, fields={has_fields}"
        ))),
    }
}

/// Extract field parameters (smoothing_length, mass) from EngineConfig.
fn field_params_from_config(config: &EngineConfig) -> Result<(f64, f64), HermesError> {
    let field_spec = config
        .ontology
        .fields
        .values()
        .next()
        .ok_or_else(|| HermesError::Config("no field species defined".to_string()))?;

    let length_scale = field_spec
        .length_scale
        .ok_or_else(|| HermesError::Config("field species requires length_scale".to_string()))?;
    let mass = field_spec
        .mass
        .ok_or_else(|| HermesError::Config("field species requires mass".to_string()))?;

    let ell = length_scale * mass;

    Ok((ell, mass))
}

/// Convert TOML halo specs to the init code's HaloConfig format.
fn halo_configs_from_config(config: &EngineConfig) -> Vec<nfw::HaloConfig> {
    let specs = &config.simulation.initialization.halos;
    if specs.is_empty() {
        return nfw::default_halo_configs();
    }

    specs
        .iter()
        .map(|s| nfw::HaloConfig {
            mass_fraction: s.mass_fraction,
            concentration: s.concentration,
        })
        .collect()
}
