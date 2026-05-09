//! Simulation runner — thin orchestrator over the pipeline.
//!
//! Routes CLI flags to the appropriate pipeline configuration:
//! headless, live, playback, or record. The simulation always runs
//! on a spawned thread; the main thread owns the event loop.

use std::sync::mpsc;
use std::thread;
use std::time::Instant;

use crate::config::Configuration;
use crate::error::HermesError;
use crate::run::cli::Cli;
use crate::run::pipeline::{
    self, PipelineMessage, SnapshotSender, spawn_disk_writer, spawn_router,
};
use crate::scenes::scene_by_name;

use std::path::PathBuf;

/// Resolve a scene argument to a TOML file path.
///
/// Tries, in order:
///   1. The argument as-is (e.g. "scenes/cosmic-web-pm.toml")
///   2. The argument with .toml appended (e.g. "scenes/cosmic-web-pm" → "scenes/cosmic-web-pm.toml")
///
/// Returns None if no file is found (falls back to legacy Scene trait).
fn resolve_scene_toml(name: &str) -> Option<PathBuf> {
    let path = PathBuf::from(name);
    if path.exists() && path.extension().is_some() {
        return Some(path);
    }

    let with_ext = PathBuf::from(format!("{name}.toml"));
    if with_ext.exists() {
        return Some(with_ext);
    }

    None
}

/// Run based on CLI arguments.
pub fn run(cli: &Cli) -> Result<(), HermesError> {
    if let Some(ref dir) = cli.playback {
        return run_playback(dir, cli);
    }

    if let Some(ref dir) = cli.resume {
        return run_resume(dir, cli);
    }

    // Config-driven path: scenes/<name>.toml → Engine.
    if let Some(scene_path) = resolve_scene_toml(&cli.scene) {
        let mut engine_config = crate::config::load_scene_config(&scene_path)?;

        // Apply CLI overrides.
        if let Some(steps) = cli.steps {
            engine_config.simulation.time.n_steps = steps;
        }
        if let Some(grid) = cli.grid {
            engine_config.simulation.grid.n_cells = grid;
        }
        engine_config.simulation.initialization.seed = cli.seed;

        if !cli.quiet {
            let cosmology =
                crate::physics::cosmology::Cosmology::from_engine_config(&engine_config)?;
            let config = Configuration::from_engine_config(&engine_config, &cosmology);
            print_header(&config, cli);
        }

        #[cfg(feature = "vis")]
        if cli.live {
            return run_engine_live(engine_config, cli);
        }

        return run_engine_headless(engine_config, cli);
    }

    // Legacy fallback: Scene trait.
    let scene = scene_by_name(&cli.scene)?;
    let config = load_config(cli, scene.default_overrides())?;

    if !cli.quiet {
        print_header(&config, cli);
    }

    #[cfg(feature = "vis")]
    if cli.live {
        return run_live(config, cli);
    }

    #[cfg(not(feature = "vis"))]
    if cli.live {
        return Err(HermesError::Config(
            "live viewer requires --features vis".to_string(),
        ));
    }

    run_headless(config, cli)
}

// ============================================================================
/// Build a Simulation from either the config-driven path (scenes/<name>.toml)
/// or the legacy Scene trait path.
fn build_simulation(
    scene_name: &str,
    config: Configuration,
    seed: u64,
) -> Result<crate::core::simulation::Simulation, HermesError> {
    if let Some(scene_path) = resolve_scene_toml(scene_name) {
        let mut engine_config = crate::config::load_scene_config(&scene_path)?;

        // Transfer any overrides from the legacy config into the engine config.
        engine_config.simulation.time.n_steps = config.time.n_steps;
        engine_config.simulation.grid.n_cells = config.simulation.n_cells();
        engine_config.simulation.grid.box_length = config.simulation.box_length;
        engine_config.simulation.initialization.seed = seed;

        crate::core::simulation::Simulation::from_engine_config(&engine_config)
    } else {
        let scene = scene_by_name(scene_name)?;
        crate::core::simulation::Simulation::from_scene(&*scene, config, seed)
    }
}

// ============================================================================
// Engine-based runners (TOML config → Engine → pipeline)
// ============================================================================

fn run_engine_headless(
    engine_config: crate::config::EngineConfig,
    cli: &Cli,
) -> Result<(), HermesError> {
    let save_dir = cli.save_directory();
    let quiet = cli.quiet;
    let scene_name = cli.scene.clone();
    let total_steps = engine_config.simulation.time.n_steps;
    let write_interval = engine_config.output.snapshots.interval.max(1);

    // Copy scene TOML into output directory for --resume.
    if let Some(ref dir) = save_dir
        && let Some(scene_path) = resolve_scene_toml(&scene_name)
    {
        let dest = PathBuf::from(dir).join("scene.toml");
        std::fs::create_dir_all(dir).ok();
        std::fs::copy(&scene_path, &dest).ok();
    }

    // Channels.
    let (sim_tx, router_rx) = mpsc::sync_channel::<PipelineMessage>(512);
    let sender = SnapshotSender::new(sim_tx);

    use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

    let multi = MultiProgress::new();

    let mut consumer_senders: Vec<pipeline::ConsumerConfig> = Vec::new();
    let mut handles: Vec<thread::JoinHandle<()>> = Vec::new();

    let sim_bar = if !quiet {
        let bar = multi.add(ProgressBar::new(total_steps as u64));
        bar.set_style(
            ProgressStyle::with_template(
                "{spinner:.cyan} [{elapsed_precise}] [{bar:40.cyan/dark.grey}] step {pos}/{len} z={msg} ({eta} remaining)",
            )
            .unwrap()
            .progress_chars("=> "),
        );
        Some(bar)
    } else {
        None
    };

    if let Some(ref dir) = save_dir {
        if !quiet {
            println!("Saving snapshots to {dir}/");
        }
        let (disk_tx, disk_rx) = mpsc::sync_channel::<PipelineMessage>(512);
        consumer_senders.push(pipeline::ConsumerConfig {
            tx: disk_tx,
            droppable: false,
        });

        let disk_progress = if !quiet {
            let total_snapshots = total_steps / write_interval + 1;
            let bar = multi.add(ProgressBar::new(total_snapshots as u64));
            bar.set_style(
                ProgressStyle::with_template(
                    "{spinner:.green}            [{bar:40.green/dark.grey}] {pos}/{len} snapshots saved",
                )
                .unwrap()
                .progress_chars("=> "),
            );
            Some(bar)
        } else {
            None
        };

        let dir_owned = dir.clone();
        handles.push(spawn_disk_writer(disk_rx, dir_owned, disk_progress));
    }

    let router_handle = spawn_router(router_rx, consumer_senders);

    let run_start = Instant::now();

    let sim_handle = thread::Builder::new()
        .name("simulation".to_string())
        .spawn(move || -> Result<(usize, f64), HermesError> {
            let (mut engine, cosmology) = crate::engine::init::init_engine(&engine_config)?;

            if !quiet {
                let progress_bar = sim_bar.unwrap();

                crate::engine::init::run_engine(
                    &mut engine,
                    &engine_config,
                    &cosmology,
                    &sender,
                    |step, scale_factor| {
                        let redshift = 1.0 / scale_factor - 1.0;
                        progress_bar.set_position(step as u64);
                        progress_bar.set_message(format!("{redshift:.1}"));
                    },
                )?;

                progress_bar.finish_and_clear();
            } else {
                crate::engine::init::run_engine(
                    &mut engine,
                    &engine_config,
                    &cosmology,
                    &sender,
                    |_, _| {},
                )?;
            }

            Ok((engine.state.step, engine.state.time))
        })
        .expect("failed to spawn simulation thread");

    let (n_steps_done, final_scale_factor) =
        sim_handle.join().expect("simulation thread panicked")?;

    let run_time = run_start.elapsed();

    let _ = router_handle.join();
    for h in handles {
        let _ = h.join();
    }

    if !quiet {
        println!(
            "Completed {} steps in {:.2}s ({:.1} ms/step)",
            n_steps_done,
            run_time.as_secs_f64(),
            run_time.as_secs_f64() * 1000.0 / n_steps_done.max(1) as f64,
        );
        println!("Final state: z = {:.1}", 1.0 / final_scale_factor - 1.0);
    }

    Ok(())
}

#[cfg(feature = "vis")]
fn run_engine_live(
    engine_config: crate::config::EngineConfig,
    cli: &Cli,
) -> Result<(), HermesError> {
    let save_dir = cli.save_directory();
    let quiet = cli.quiet;
    let scene_name = cli.scene.clone();
    let total_steps = engine_config.simulation.time.n_steps;
    let write_interval = engine_config.output.snapshots.interval.max(1);

    // Copy scene TOML into output directory.
    if let Some(ref dir) = save_dir
        && let Some(scene_path) = resolve_scene_toml(&scene_name)
    {
        let dest = PathBuf::from(dir).join("scene.toml");
        std::fs::create_dir_all(dir).ok();
        std::fs::copy(&scene_path, &dest).ok();
    }

    let cosmology_for_vis =
        crate::physics::cosmology::Cosmology::from_engine_config(&engine_config)?;
    let vis_config =
        Configuration::from_engine_config(&engine_config, &cosmology_for_vis).visualization;
    let box_length = engine_config.simulation.grid.box_length;

    let (sim_tx, router_rx) = mpsc::sync_channel::<PipelineMessage>(512);
    let sender = SnapshotSender::new(sim_tx);

    let mut consumer_senders: Vec<pipeline::ConsumerConfig> = Vec::new();
    let mut handles: Vec<thread::JoinHandle<()>> = Vec::new();

    // Disk writer (if saving).
    if let Some(ref dir) = save_dir {
        if !quiet {
            println!("Saving snapshots to {dir}/");
        }
        let (disk_tx, disk_rx) = mpsc::sync_channel::<PipelineMessage>(512);
        consumer_senders.push(pipeline::ConsumerConfig {
            tx: disk_tx,
            droppable: false,
        });
        let dir_owned = dir.clone();
        handles.push(spawn_disk_writer(disk_rx, dir_owned, None));
    }

    // Precompute thread → viewer.
    let (precompute_tx, precompute_rx) = mpsc::sync_channel::<PipelineMessage>(512);
    consumer_senders.push(pipeline::ConsumerConfig {
        tx: precompute_tx,
        droppable: true,
    });

    let (frame_tx, frame_rx) = mpsc::sync_channel::<pipeline::ViewerMessage>(4);
    handles.push(pipeline::spawn_precompute(
        precompute_rx,
        frame_tx,
        box_length,
        vis_config.clone(),
    ));

    let router_handle = spawn_router(router_rx, consumer_senders);

    // Simulation on spawned thread.
    let sim_handle = thread::Builder::new()
        .name("simulation".to_string())
        .spawn(move || -> Result<(), HermesError> {
            let (mut engine, cosmology) = crate::engine::init::init_engine(&engine_config)?;

            crate::engine::init::run_engine(
                &mut engine,
                &engine_config,
                &cosmology,
                &sender,
                |_, _| {},
            )
        })
        .expect("failed to spawn simulation thread");

    // Main thread: viewer event loop.
    if !quiet {
        println!("Starting live viewer + simulation...");
        println!("(close the viewer window to exit)");
        println!();
    }
    pipeline::run_viewer_main_thread(frame_rx, vis_config);

    let _ = sim_handle.join().expect("simulation thread panicked");
    let _ = router_handle.join();
    for h in handles {
        let _ = h.join();
    }

    Ok(())
}

// ============================================================================
// Legacy headless: simulation on spawned thread, main blocks on join
// ============================================================================

fn run_headless(config: Configuration, cli: &Cli) -> Result<(), HermesError> {
    let save_dir = cli.save_directory();
    let seed = cli.seed;
    let quiet = cli.quiet;
    let scene_name = cli.scene.clone();

    // Write scene.toml into the output directory so --resume can find it.
    if let Some(ref dir) = save_dir
        && let Some(scene_path) = resolve_scene_toml(&scene_name)
    {
        let dest = PathBuf::from(dir).join("scene.toml");
        std::fs::create_dir_all(dir).ok();
        std::fs::copy(&scene_path, &dest).ok();
    }
    let total_steps = config.time.n_steps;
    // Simulation → Router channel.
    let (sim_tx, router_rx) = mpsc::sync_channel::<PipelineMessage>(512);
    let sender = SnapshotSender::new(sim_tx);

    // Multi-progress for parallel bars.
    use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

    let multi = MultiProgress::new();

    // Build consumer list.
    let mut consumer_senders: Vec<pipeline::ConsumerConfig> = Vec::new();
    let mut handles: Vec<thread::JoinHandle<()>> = Vec::new();

    // Simulation progress bar first (top line).
    let sim_bar = if !quiet {
        let bar = multi.add(ProgressBar::new(total_steps as u64));
        bar.set_style(
            ProgressStyle::with_template(
                "{spinner:.cyan} [{elapsed_precise}] [{bar:40.cyan/dark.grey}] step {pos}/{len} z={msg} ({eta} remaining)",
            )
            .unwrap()
            .progress_chars("=> "),
        );
        Some(bar)
    } else {
        None
    };

    if let Some(ref dir) = save_dir {
        if !quiet {
            println!("Saving snapshots to {dir}/");
        }
        let (disk_tx, disk_rx) = mpsc::sync_channel::<PipelineMessage>(512);
        consumer_senders.push(pipeline::ConsumerConfig {
            tx: disk_tx,
            droppable: false,
        });

        // Disk progress bar second (bottom line).
        let disk_progress = if !quiet {
            let write_interval = config.output.write_interval.max(1);
            let total_snapshots = total_steps / write_interval + 1;
            let bar = multi.add(ProgressBar::new(total_snapshots as u64));
            bar.set_style(
                ProgressStyle::with_template(
                    "{spinner:.green}            [{bar:40.green/dark.grey}] {pos}/{len} snapshots saved",
                )
                .unwrap()
                .progress_chars("=> "),
            );
            Some(bar)
        } else {
            None
        };

        let dir_owned = dir.clone();
        handles.push(spawn_disk_writer(disk_rx, dir_owned, disk_progress));
    }

    let router_handle = spawn_router(router_rx, consumer_senders);

    // Spawn simulation.
    let run_start = Instant::now();

    let sim_handle = thread::Builder::new()
        .name("simulation".to_string())
        .spawn(
            move || -> Result<crate::core::simulation::Simulation, HermesError> {
                let mut sim = build_simulation(&scene_name, config, seed)?;

                if !quiet {
                    let progress_bar = sim_bar.unwrap();

                    sim.run_into_pipeline(&sender, |step, scale_factor| {
                        let redshift = 1.0 / scale_factor - 1.0;
                        progress_bar.set_position(step as u64);
                        progress_bar.set_message(format!("{redshift:.1}"));
                    })?;

                    progress_bar.finish_and_clear();
                } else {
                    sim.run_into_pipeline(&sender, |_, _| {})?;
                }

                Ok(sim)
            },
        )
        .expect("failed to spawn simulation thread");

    // Main thread: block waiting for simulation.
    let sim = sim_handle.join().expect("simulation thread panicked")?;

    let run_time = run_start.elapsed();

    let _ = router_handle.join();
    for h in handles {
        let _ = h.join();
    }

    if !quiet {
        println!(
            "Completed {} steps in {:.2}s ({:.1} ms/step)",
            sim.step,
            run_time.as_secs_f64(),
            run_time.as_secs_f64() * 1000.0 / sim.step as f64,
        );

        if let Some(diag) = sim.latest_diagnostics() {
            println!();
            println!("Final state (z = {:.1}):", 1.0 / sim.scale_factor - 1.0);
            println!("  Mass:     {:.4e} M_☉", diag.mass_total);
            println!("  |P|:      {:.4e}", diag.momentum_magnitude());
            println!("  E_kin:    {:.4e}", diag.energy_kinetic);
            println!("  E_pot:    {:.4e}", diag.energy_potential);
            println!("  |L|:      {:.4e}", diag.angular_momentum_magnitude());
        }
    }

    Ok(())
}

// ============================================================================
// Live: simulation on spawned thread, viewer on main
// ============================================================================

#[cfg(feature = "vis")]
fn run_live(config: Configuration, cli: &Cli) -> Result<(), HermesError> {
    let save_dir = cli.save_directory();
    let seed = cli.seed;
    let quiet = cli.quiet;
    let scene_name = cli.scene.clone();
    let box_length = config.simulation.box_length;
    let vis_config = config.visualization.clone();

    if !quiet {
        println!("Starting live viewer + simulation...");
        println!("(close the viewer window to exit)");
        println!();
    }

    // Simulation → Router.
    let (sim_tx, router_rx) = mpsc::sync_channel::<PipelineMessage>(512);
    let sender = SnapshotSender::new(sim_tx);

    let mut consumer_senders: Vec<pipeline::ConsumerConfig> = Vec::new();
    let mut handles: Vec<thread::JoinHandle<()>> = Vec::new();

    // Disk writer (optional).
    if let Some(ref dir) = save_dir {
        let (disk_tx, disk_rx) = mpsc::sync_channel::<PipelineMessage>(512);
        consumer_senders.push(pipeline::ConsumerConfig {
            tx: disk_tx,
            droppable: false,
        });
        handles.push(spawn_disk_writer(disk_rx, dir.clone(), None));
    }

    // Precompute → Viewer channel.
    let (precompute_tx, precompute_rx) = mpsc::sync_channel::<PipelineMessage>(4);
    consumer_senders.push(pipeline::ConsumerConfig {
        tx: precompute_tx,
        droppable: true,
    });

    let (frame_tx, frame_rx) = mpsc::sync_channel::<pipeline::ViewerMessage>(4);
    handles.push(pipeline::spawn_precompute(
        precompute_rx,
        frame_tx,
        box_length,
        vis_config.clone(),
    ));

    // Router.
    let router_handle = spawn_router(router_rx, consumer_senders);

    // Simulation thread.
    let sim_handle = thread::Builder::new()
        .name("simulation".to_string())
        .spawn(move || -> Result<(), HermesError> {
            let mut sim = build_simulation(&scene_name, config, seed)?;
            sim.run_into_pipeline(&sender, |_, _| {})?;

            Ok(())
        })
        .expect("failed to spawn simulation thread");

    // Main thread: viewer event loop.
    pipeline::run_viewer_main_thread(frame_rx, vis_config);

    // Clean up.
    let _ = sim_handle.join();
    let _ = router_handle.join();
    for h in handles {
        let _ = h.join();
    }

    Ok(())
}

// ============================================================================
// Resume: continue from the last snapshot
// ============================================================================

fn run_resume(dir: &str, cli: &Cli) -> Result<(), HermesError> {
    use crate::io::snapshot::load_snapshot;

    // Find and load the last snapshot.
    let paths = pipeline::find_snapshot_paths(dir);
    let last_path = paths
        .last()
        .ok_or_else(|| HermesError::Config(format!("no snapshots found in {dir}/")))?;

    let snapshot = load_snapshot(last_path)?;
    let particles = snapshot
        .to_particles()
        .ok_or_else(|| HermesError::Config("resume requires particle content".to_string()))?;

    let a_start = snapshot.scale_factor;

    if !cli.quiet {
        println!(
            "Resuming from {} (step {}, a = {:.4}, z = {:.1})",
            last_path.display(),
            snapshot.step,
            a_start,
            1.0 / a_start - 1.0,
        );
    }

    // Load config from <dir>/scene.toml (written when the scene was saved).
    // Falls back to the --scene flag or legacy path if not found.
    let scene_toml_in_dir = PathBuf::from(dir).join("scene.toml");
    let mut config = if scene_toml_in_dir.exists() {
        let engine_config = crate::config::load_scene_config(&scene_toml_in_dir)?;
        let cosmology = crate::physics::cosmology::Cosmology::from_engine_config(&engine_config)?;
        Configuration::from_engine_config(&engine_config, &cosmology)
    } else if let Some(toml_path) = resolve_scene_toml(&cli.scene) {
        let engine_config = crate::config::load_scene_config(&toml_path)?;
        let cosmology = crate::physics::cosmology::Cosmology::from_engine_config(&engine_config)?;
        Configuration::from_engine_config(&engine_config, &cosmology)
    } else {
        let scene = scene_by_name(&cli.scene)?;
        load_config(cli, scene.default_overrides())?
    };

    // Override the initial scale factor to the snapshot's value.
    config.time.scale_factor_range[0] = a_start;

    // --steps is required for resume: how many additional steps to run.
    if let Some(steps) = cli.steps {
        config.time.n_steps = steps;
    }

    if !cli.quiet {
        println!(
            "Evolving a = {:.4} → {:.4} in {} steps ({})",
            config.time.scale_factor_initial(),
            config.time.scale_factor_final(),
            config.time.n_steps,
            config.time.scale_factor_stepping,
        );
        println!();
    }

    let start_step = snapshot.step;

    let mut sim = crate::core::simulation::Simulation::resume_particles(
        particles,
        a_start,
        start_step,
        config.clone(),
    )?;

    // Save to the same directory we resumed from (unless --save overrides).
    let save_dir = match &cli.save {
        Some(Some(d)) => Some(d.clone()),
        Some(None) => Some(dir.to_string()),
        None => None,
    };
    let quiet = cli.quiet;
    let total_steps = config.time.n_steps;

    let (sim_tx, router_rx) = mpsc::sync_channel::<PipelineMessage>(512);
    let sender = SnapshotSender::new(sim_tx);

    use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
    let multi = MultiProgress::new();

    let mut consumer_senders: Vec<pipeline::ConsumerConfig> = Vec::new();
    let mut handles: Vec<thread::JoinHandle<()>> = Vec::new();

    let sim_bar = if !quiet {
        let bar = multi.add(ProgressBar::new(total_steps as u64));
        bar.set_style(
            ProgressStyle::with_template(
                "{spinner:.cyan} [{elapsed_precise}] [{bar:40.cyan/dark.grey}] step {pos}/{len} a={msg} ({eta} remaining)",
            )
            .unwrap()
            .progress_chars("=> "),
        );
        Some(bar)
    } else {
        None
    };

    if let Some(ref dir) = save_dir {
        if !quiet {
            println!("Saving snapshots to {dir}/");
        }
        let (disk_tx, disk_rx) = mpsc::sync_channel::<PipelineMessage>(512);
        consumer_senders.push(pipeline::ConsumerConfig {
            tx: disk_tx,
            droppable: false,
        });

        let disk_progress = if !quiet {
            let write_interval = config.output.write_interval.max(1);
            let total_snapshots = total_steps / write_interval + 1;
            let bar = multi.add(ProgressBar::new(total_snapshots as u64));
            bar.set_style(
                ProgressStyle::with_template(
                    "{spinner:.green}            [{bar:40.green/dark.grey}] {pos}/{len} snapshots saved",
                )
                .unwrap()
                .progress_chars("=> "),
            );
            Some(bar)
        } else {
            None
        };

        let dir_owned = dir.clone();
        handles.push(spawn_disk_writer(disk_rx, dir_owned, disk_progress));
    }

    let router_handle = spawn_router(router_rx, consumer_senders);
    let run_start = Instant::now();

    let sim_handle = thread::Builder::new()
        .name("simulation".to_string())
        .spawn(
            move || -> Result<crate::core::simulation::Simulation, HermesError> {
                if !quiet {
                    let progress_bar = sim_bar.unwrap();
                    sim.run_into_pipeline(&sender, |step, scale_factor| {
                        progress_bar.set_position(step as u64);
                        progress_bar.set_message(format!("{scale_factor:.4}"));
                    })?;
                    progress_bar.finish_and_clear();
                } else {
                    sim.run_into_pipeline(&sender, |_, _| {})?;
                }

                Ok(sim)
            },
        )
        .expect("failed to spawn simulation thread");

    let sim = sim_handle.join().expect("simulation thread panicked")?;
    let run_time = run_start.elapsed();

    let _ = router_handle.join();
    for h in handles {
        let _ = h.join();
    }

    if !quiet {
        println!(
            "Completed {} steps in {:.2}s ({:.1} ms/step)",
            sim.step,
            run_time.as_secs_f64(),
            run_time.as_secs_f64() * 1000.0 / sim.step as f64,
        );
    }

    Ok(())
}

// ============================================================================
// Playback / Record
// ============================================================================

fn run_playback(dir: &str, cli: &Cli) -> Result<(), HermesError> {
    if let Some(ref output_path) = cli.record {
        return record_to_gif(dir, output_path, cli);
    }

    #[cfg(not(feature = "vis"))]
    {
        let _ = (dir, cli);
        Err(HermesError::Config(
            "playback requires --features vis".to_string(),
        ))
    }

    #[cfg(feature = "vis")]
    {
        if !cli.quiet {
            println!("Playback from {dir}/");
        }

        // Infer scene TOML from data directory: <dir>.toml is the sibling config.
        // Falls back to global defaults if no TOML is found.
        let scene_toml = PathBuf::from(format!("{dir}.toml"));
        let vis_config = if scene_toml.exists() {
            let engine_config = crate::config::load_scene_config(&scene_toml)?;
            let cosmology =
                crate::physics::cosmology::Cosmology::from_engine_config(&engine_config)?;
            Configuration::from_engine_config(&engine_config, &cosmology).visualization
        } else if let Some(toml_path) = resolve_scene_toml(&cli.scene) {
            let engine_config = crate::config::load_scene_config(&toml_path)?;
            let cosmology =
                crate::physics::cosmology::Cosmology::from_engine_config(&engine_config)?;
            Configuration::from_engine_config(&engine_config, &cosmology).visualization
        } else {
            crate::config::VisualizationConfig::default()
        };

        pipeline::run_playback_viewer(dir, cli.fps, &vis_config)
    }
}

fn record_to_gif(dir: &str, output_path: &str, cli: &Cli) -> Result<(), HermesError> {
    use crate::colormap::colormap_hot;
    use crate::io::snapshot::load_snapshot;

    let snapshot_paths = pipeline::find_snapshot_paths(dir);
    let total = snapshot_paths.len();
    if total == 0 {
        return Err(HermesError::Config(format!("no snapshots found in {dir}/")));
    }

    if !cli.quiet {
        println!("Recording {total} frames to {output_path}...");
    }

    // Infer vis config from sibling TOML of the data directory.
    let scene_toml = PathBuf::from(format!("{dir}.toml"));
    let vis = if scene_toml.exists() {
        let engine_config = crate::config::load_scene_config(&scene_toml)?;
        let cosmology = crate::physics::cosmology::Cosmology::from_engine_config(&engine_config)?;
        Configuration::from_engine_config(&engine_config, &cosmology).visualization
    } else {
        crate::config::VisualizationConfig::default()
    };

    let width = vis.gif_resolution;
    let height = vis.gif_resolution;
    let point_radius = vis.gif_point_radius;

    let first = load_snapshot(&snapshot_paths[0])?;
    let box_length = first
        .positions()
        .unwrap_or(&[])
        .iter()
        .flat_map(|pos: &morphis::vector::Vector<3>| {
            (1..=3).map(move |d| pos.component(&[d]).abs())
        })
        .fold(0.0_f64, f64::max)
        * 1.1;
    let scale = 1.0 / box_length;

    if let Some(parent) = std::path::Path::new(output_path).parent() {
        std::fs::create_dir_all(parent).ok();
    }

    let output_file = std::fs::File::create(output_path)
        .map_err(|e| HermesError::Config(format!("failed to create {output_path}: {e}")))?;
    let mut encoder = image::codecs::gif::GifEncoder::new(output_file);
    let frame_delay = image::Delay::from_saturating_duration(std::time::Duration::from_millis(
        1000 / cli.fps.max(1),
    ));

    let progress = if !cli.quiet {
        use indicatif::{ProgressBar, ProgressStyle};
        let progress_bar = ProgressBar::new(total as u64);
        progress_bar.set_style(
            ProgressStyle::with_template(
                "{spinner:.cyan} [{bar:40.cyan/dark.grey}] {pos}/{len} frames",
            )
            .unwrap()
            .progress_chars("=> "),
        );
        Some(progress_bar)
    } else {
        None
    };

    for path in &snapshot_paths {
        let snapshot = load_snapshot(path)?;

        let speeds: Vec<f64> = snapshot
            .momenta()
            .unwrap()
            .iter()
            .map(|mom| mom.norm())
            .collect();
        let speed_max = speeds.iter().copied().fold(1e-30_f64, f64::max);
        let speed_min = speeds.iter().copied().fold(f64::MAX, f64::min);
        let speed_range = (speed_max - speed_min).max(1e-30);

        let mut pixels = vec![0_u8; (width * height * 4) as usize];
        for pixel in pixels.chunks_exact_mut(4) {
            pixel[3] = 255;
        }

        for (pos, &speed) in snapshot.positions().unwrap().iter().zip(speeds.iter()) {
            let x_norm = pos.component(&[1]) * scale;
            let y_norm = pos.component(&[2]) * scale;

            let pixel_x = (x_norm * width as f64) as i32;
            let pixel_y = (y_norm * height as f64) as i32;

            let normalized = ((speed - speed_min) / speed_range).clamp(0.0, 1.0);
            let color = colormap_hot(normalized);
            let r = (color[0] * 255.0) as u8;
            let g = (color[1] * 255.0) as u8;
            let b = (color[2] * 255.0) as u8;

            for dy in -point_radius..=point_radius {
                for dx in -point_radius..=point_radius {
                    let px = pixel_x + dx;
                    let py = pixel_y + dy;
                    if px >= 0 && px < width as i32 && py >= 0 && py < height as i32 {
                        let offset = ((py as u32 * width + px as u32) * 4) as usize;
                        pixels[offset] = pixels[offset].saturating_add(r);
                        pixels[offset + 1] = pixels[offset + 1].saturating_add(g);
                        pixels[offset + 2] = pixels[offset + 2].saturating_add(b);
                    }
                }
            }
        }

        let frame = image::Frame::from_parts(
            image::RgbaImage::from_raw(width, height, pixels).unwrap(),
            0,
            0,
            frame_delay,
        );

        encoder
            .encode_frame(frame)
            .map_err(|e| HermesError::Config(format!("GIF encode failed: {e}")))?;

        if let Some(ref progress_bar) = progress {
            progress_bar.inc(1);
        }
    }

    if let Some(progress_bar) = progress {
        progress_bar.finish_and_clear();
    }

    if !cli.quiet {
        println!("Saved {output_path} ({total} frames)");
    }

    Ok(())
}

// ============================================================================
// Config loading
// ============================================================================

fn load_config(
    cli: &Cli,
    scene_defaults: Option<toml::Value>,
) -> Result<Configuration, HermesError> {
    let file_override = if let Some(ref path) = cli.config {
        let content = std::fs::read_to_string(path)?;
        let value: toml::Value = toml::from_str(&content)
            .map_err(|e| HermesError::Config(format!("failed to parse {path}: {e}")))?;
        Some(value)
    } else {
        None
    };

    let mut overrides = toml::map::Map::new();

    if let Some(steps) = cli.steps {
        let mut time = toml::map::Map::new();
        time.insert("n_steps".to_string(), toml::Value::Integer(steps as i64));
        overrides.insert("time".to_string(), toml::Value::Table(time));
    }
    {
        let mut sim = toml::map::Map::new();
        if let Some(particles) = cli.particles {
            sim.insert(
                "n_particles".to_string(),
                toml::Value::Integer(particles as i64),
            );
        }
        if let Some(grid) = cli.grid {
            sim.insert("n_grid".to_string(), toml::Value::Integer(grid as i64));
        }
        if !sim.is_empty() {
            overrides.insert("simulation".to_string(), toml::Value::Table(sim));
        }
    }

    let programmatic = if overrides.is_empty() {
        None
    } else {
        Some(toml::Value::Table(overrides))
    };

    // Four-tier merge: global defaults → scene defaults → user file → CLI overrides.
    // build_configuration does: defaults → config_file → overrides.
    // We insert scene defaults as the config_file tier, and merge the actual
    // user file into overrides if both are present.
    match (scene_defaults, file_override) {
        (Some(scene), Some(file)) => {
            // Scene defaults as first override, then user file, then CLI.
            let mut combined = scene;
            crate::config::deep_merge_public(&mut combined, &file);
            if let Some(ref prog) = programmatic {
                crate::config::deep_merge_public(&mut combined, prog);
            }
            crate::config::build_configuration(Some(&combined), None)
        }
        (Some(scene), None) => {
            crate::config::build_configuration(Some(&scene), programmatic.as_ref())
        }
        (None, file) => crate::config::build_configuration(file.as_ref(), programmatic.as_ref()),
    }
}

fn print_header(config: &Configuration, cli: &Cli) {
    println!("Hermes — {}", cli.scene);
    println!("{}", "=".repeat(40));
    println!("Grid:       {}³ cells", config.simulation.n_cells());
    println!(
        "Particles:  {}³ = {}",
        config.simulation.n_particles,
        config.simulation.n_particles.pow(3)
    );
    println!(
        "Box:        {:.0} kpc ({:.0} Mpc)",
        config.simulation.box_length,
        config.simulation.box_length / 1000.0
    );
    println!(
        "Redshift:   z = {:.0} → z = {:.1}",
        1.0 / config.time.scale_factor_initial() - 1.0,
        1.0 / config.time.scale_factor_final() - 1.0,
    );
    println!("Steps:      {}", config.time.n_steps);
    println!("Seed:       {}", cli.seed);
    println!();
}
