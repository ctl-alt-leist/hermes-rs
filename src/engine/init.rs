/// Engine initialization and stepping: config → Engine → run.
///
/// Reads an EngineConfig, initializes a SimulationState from the
/// init library, constructs sectors and solver from the ontology,
/// and assembles an Engine. The stepping loop runs the engine
/// through a scale-factor schedule and sends snapshots through
/// the pipeline.
use std::collections::BTreeMap;
use std::sync::Arc;

use crate::config::EngineConfig;
use crate::engine::Engine;
use crate::engine::sector::Sector;
use crate::engine::sector::gross_pitaevskii::GrossPitaevskiiSector;
use crate::engine::sector::particle::ParticleSector;
use crate::engine::sector::schrodinger::SchrodingerSector;
use crate::engine::solver::GravitySolver;
use crate::engine::state::{FieldEntry, SimulationState};
use crate::error::HermesError;
use crate::io::snapshot::Snapshot;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::initial::{nfw, nfw_field, zeldovich, zeldovich_field};
use crate::physics::integrator::scale_factor_schedule;
use crate::run::pipeline::SnapshotSender;

// ============================================================================
// Engine construction from config
// ============================================================================

/// Build an Engine from an EngineConfig.
///
/// Initializes the simulation state (fields and particles) from the
/// init library, constructs sectors from the ontology, and wires up
/// the gravity solver if couplings require it.
pub fn init_engine(config: &EngineConfig) -> Result<(Engine, Cosmology), HermesError> {
    let cosmology = Cosmology::from_engine_config(config)?;

    let n_cells = config.simulation.grid.n_cells;
    let box_length = config.simulation.grid.box_length;
    let grid = Grid::new(n_cells, box_length);
    let morphis_grid = morphis::grid::Grid::<3>::new(n_cells, box_length);

    let scale_factor_initial = config
        .simulation
        .time
        .scale_factor_range
        .map(|r| r[0])
        .unwrap_or(1.0);

    let seed = config.simulation.initialization.seed;
    let method = config.simulation.initialization.method.as_str();

    // Initialize species into SimulationState.
    let (fields, particles) = init_species(
        config,
        &grid,
        &cosmology,
        scale_factor_initial,
        seed,
        method,
    )?;

    let state = SimulationState {
        particles,
        fields,
        grid: grid.clone(),
        morphis_grid,
        time: scale_factor_initial,
        step: 0,
    };

    // Construct sectors from ontology.
    let sectors = build_sectors(config);

    // Construct solver from couplings.
    let solver = if config.ontology.has_gravity() {
        if state.has_particles() {
            Some(GravitySolver::with_particles(morphis_grid, grid))
        } else {
            Some(GravitySolver::new(morphis_grid))
        }
    } else {
        None
    };

    let engine = Engine::new(state, sectors, solver, Some(cosmology.clone()));

    Ok((engine, cosmology))
}

// ============================================================================
// Species initialization
// ============================================================================

/// Initialize all field and particle species from the config.
///
/// Each species can override the global initialization method and
/// parameters via its own `[initialization]` sub-table. The global
/// `[simulation.initialization]` provides defaults.
#[allow(clippy::type_complexity)]
fn init_species(
    config: &EngineConfig,
    grid: &Grid,
    cosmology: &Cosmology,
    scale_factor_initial: f64,
    seed: u64,
    global_method: &str,
) -> Result<
    (
        BTreeMap<String, FieldEntry>,
        BTreeMap<String, crate::physics::particles::Particles>,
    ),
    HermesError,
> {
    let global_init = &config.simulation.initialization;
    let n_total_species = config.ontology.fields.len() + config.ontology.particles.len();
    let default_density_fraction = if n_total_species > 0 {
        1.0 / n_total_species as f64
    } else {
        1.0
    };

    let mut fields = BTreeMap::new();
    let mut particles = BTreeMap::new();

    // Initialize field species.
    for (k, (name, spec)) in config.ontology.fields.iter().enumerate() {
        let mass = spec
            .mass
            .ok_or_else(|| HermesError::Config(format!("field '{name}' requires mass")))?;
        let length_scale = spec.length_scale.unwrap_or(1.0);
        let ell = length_scale * mass;
        let field_seed = seed + k as u64;

        let species_init = spec.initialization.as_ref();

        // Resolve per-species overrides with global fallbacks.
        let method = species_init
            .and_then(|s| s.method.as_deref())
            .unwrap_or(global_method);
        let density_fraction = species_init
            .and_then(|s| s.density_fraction)
            .unwrap_or(default_density_fraction);
        let spectrum = species_init
            .and_then(|s| s.spectrum.as_deref())
            .unwrap_or(global_init.spectrum.as_str());
        let amplitude = species_init
            .and_then(|s| s.perturbation_amplitude)
            .unwrap_or(global_init.perturbation_amplitude);
        let band_pass = species_init
            .and_then(|s| s.band_pass)
            .unwrap_or(global_init.band_pass);

        let params = crate::core::content::FieldParams {
            smoothing_length: ell,
            mass_alpha: mass,
        };

        let data = match method {
            "zeldovich" => {
                if spectrum == "random" {
                    zeldovich_field::random_density_field(
                        grid,
                        cosmology,
                        &params,
                        scale_factor_initial,
                        amplitude,
                        band_pass,
                        field_seed,
                        density_fraction,
                    )
                } else {
                    zeldovich_field::zeldovich_wavefunction(
                        grid,
                        cosmology,
                        &params,
                        scale_factor_initial,
                        amplitude,
                        field_seed,
                        density_fraction,
                    )?
                }
            }

            "nfw-group" => {
                let halos = halo_configs_from_config(config);
                nfw_field::colliding_halos_field(
                    grid,
                    cosmology,
                    &params,
                    scale_factor_initial,
                    field_seed,
                    &halos,
                    density_fraction,
                )
            }

            "uniform" => {
                let morphis_grid = morphis::grid::Grid::<3>::new(grid.n_cells, grid.box_length);
                let g = morphis::metric::euclidean::<3>();
                let rho_mean = cosmology.density_matter() * density_fraction;
                let uniform_amplitude = (rho_mean / mass).sqrt();
                morphis::even_field::EvenField::from_fn(&morphis_grid, g, |_| {
                    (uniform_amplitude, 0.0)
                })
            }

            "gaussian-packet" => {
                let center = species_init
                    .and_then(|s| s.center)
                    .or(global_init.center)
                    .unwrap_or([0.5, 0.5, 0.5]);
                let width = species_init
                    .and_then(|s| s.width)
                    .or(global_init.width)
                    .unwrap_or(0.05);
                let momentum = species_init
                    .and_then(|s| s.momentum)
                    .or(global_init.momentum)
                    .unwrap_or([0.0, 0.0, 0.0]);

                let morphis_grid = morphis::grid::Grid::<3>::new(grid.n_cells, grid.box_length);
                let g = morphis::metric::euclidean::<3>();
                let rho_mean = cosmology.density_matter() * density_fraction;
                let nu = ell / mass;
                let box_length = grid.box_length;

                morphis::even_field::EvenField::from_fn(&morphis_grid, g, |pos| {
                    let dx = pos[0] - center[0] * box_length;
                    let dy = pos[1] - center[1] * box_length;
                    let dz = pos[2] - center[2] * box_length;
                    let r2 = dx * dx + dy * dy + dz * dz;
                    let sigma = width * box_length;

                    let rho = rho_mean * (-r2 / (2.0 * sigma * sigma)).exp()
                        / (2.0 * std::f64::consts::PI * sigma * sigma).powf(1.5);
                    let amplitude = (rho / mass).sqrt().max(1e-30);

                    let phase =
                        (momentum[0] * pos[0] + momentum[1] * pos[1] + momentum[2] * pos[2]) / nu;

                    (amplitude * phase.cos(), amplitude * phase.sin())
                })
            }

            _ => {
                return Err(HermesError::Config(format!(
                    "unsupported field initialization method '{method}' for species '{name}'"
                )));
            }
        };

        fields.insert(
            name.clone(),
            FieldEntry {
                data,
                smoothing_length: ell,
                mass,
                self_interaction: spec.self_interaction,
            },
        );
    }

    // Initialize particle species.
    let n_fields = config.ontology.fields.len();
    for (k, (name, spec)) in config.ontology.particles.iter().enumerate() {
        let particle_seed = seed + (n_fields + k) as u64;

        let species_init = spec.initialization.as_ref();
        let method = species_init
            .and_then(|s| s.method.as_deref())
            .unwrap_or(global_method);
        let density_fraction = species_init
            .and_then(|s| s.density_fraction)
            .unwrap_or(default_density_fraction);

        let p = match method {
            "zeldovich" => zeldovich::zeldovich_init(
                spec.n,
                grid,
                cosmology,
                scale_factor_initial,
                particle_seed,
            )?,

            "nfw-group" => {
                let halos = halo_configs_from_config(config);
                nfw::colliding_halos_init(
                    spec.n,
                    grid,
                    cosmology,
                    scale_factor_initial,
                    particle_seed,
                    &halos,
                    density_fraction,
                )?
            }

            _ => {
                return Err(HermesError::Config(format!(
                    "unsupported particle initialization method '{method}' for species '{name}'"
                )));
            }
        };

        particles.insert(name.clone(), p);
    }

    Ok((fields, particles))
}

// ============================================================================
// Sector construction
// ============================================================================

/// Build sectors from the ontology.
fn build_sectors(config: &EngineConfig) -> Vec<Box<dyn Sector>> {
    let mut sectors: Vec<Box<dyn Sector>> = Vec::new();

    for (name, spec) in &config.ontology.fields {
        let sector: Box<dyn Sector> = if spec.self_interaction.is_some() {
            Box::new(GrossPitaevskiiSector::new(name.clone()))
        } else {
            Box::new(SchrodingerSector::new(name.clone()))
        };
        sectors.push(sector);
    }

    for name in config.ontology.particles.keys() {
        sectors.push(Box::new(ParticleSector::new(name.clone())));
    }

    sectors
}

// ============================================================================
// Stepping loop
// ============================================================================

/// Run the engine through its full schedule, sending snapshots to the pipeline.
pub fn run_engine(
    engine: &mut Engine,
    config: &EngineConfig,
    cosmology: &Cosmology,
    sender: &SnapshotSender,
    on_step: impl Fn(usize, f64),
) -> Result<(), HermesError> {
    let time_config = &config.simulation.time;
    let write_interval = config.output.snapshots.interval.max(1);

    let (schedule, dt_for_finalize) = if let Some(range) = time_config.scale_factor_range {
        let schedule = scale_factor_schedule(
            range[0],
            range[1],
            time_config.n_steps,
            &time_config.stepping,
        );
        (schedule, None)
    } else if let Some(range) = time_config.time_range {
        let dt = (range[1] - range[0]) / time_config.n_steps as f64;
        let schedule: Vec<f64> = (0..=time_config.n_steps)
            .map(|k| range[0] + k as f64 * dt)
            .collect();
        (schedule, Some(dt))
    } else {
        return Err(HermesError::Config(
            "time config requires scale_factor_range or time_range".to_string(),
        ));
    };

    let n_steps = time_config.n_steps;

    // Send initial snapshot.
    let snapshot = Arc::new(Snapshot::capture_from_state(&engine.state, 0, schedule[0]));
    sender.send(snapshot);

    // Stepping loop.
    for k in 0..n_steps {
        let a_prev = schedule[k];
        let a_next = schedule[k + 1];

        // Compute dt from scale factor step.
        let dt = if let Some(fixed_dt) = dt_for_finalize {
            // Static spacetime: uniform dt.
            fixed_dt
        } else {
            // FLRW: dt from scale factor and Hubble.
            let a_mid = (a_prev + a_next) / 2.0;
            (a_next - a_prev) / (a_mid * cosmology.hubble_parameter(a_mid))
        };

        let scale_factor = (a_prev + a_next) / 2.0;
        engine.step(scale_factor, dt)?;

        engine.state.time = a_next;

        on_step(engine.state.step, a_next);

        // Send snapshot at write interval or final step.
        let is_final = k + 1 == n_steps;
        if engine.state.step.is_multiple_of(write_interval) || is_final {
            if is_final {
                engine.finalize(scale_factor, dt)?;
            }
            let snapshot = Arc::new(Snapshot::capture_from_state(
                &engine.state,
                engine.state.step,
                a_next,
            ));
            sender.send(snapshot);
        }
    }

    sender.done();

    Ok(())
}

// ============================================================================
// Helpers
// ============================================================================

/// Convert TOML halo specs to init code's HaloConfig format.
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
