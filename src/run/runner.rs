//! Simulation runner — orchestrates init, stepping, observers, and post-run.

use std::time::Instant;

use crate::config::{Configuration, build_configuration};
use crate::error::HermesError;
use crate::io::observer::{FileObserver, Observer};
use crate::run::cli::Cli;
use crate::scenes::scene_by_name;

/// Run a simulation based on CLI arguments.
pub fn run(cli: &Cli) -> Result<(), HermesError> {
    // Handle playback mode separately.
    if let Some(ref dir) = cli.playback {
        return run_playback(dir, cli);
    }

    // Build configuration.
    let config = load_config(cli)?;

    if !cli.quiet {
        print_header(&config, cli);
    }

    // Look up the scene.
    let scene = scene_by_name(&cli.scene)?;

    // Initialize.
    let start = Instant::now();
    let mut sim = scene.initialize(&config, cli.seed)?;

    if !cli.quiet {
        println!("Initialized in {:.2}s", start.elapsed().as_secs_f64());
        println!();
    }

    // Determine save directory.
    let save_dir = cli.save_directory();

    // Live mode: viewer on main thread, simulation on background thread.
    #[cfg(feature = "vis")]
    if cli.live {
        if !cli.quiet {
            println!("Starting live viewer + simulation...");
            println!("(close the viewer window to exit)");
            println!();
        }

        return crate::vis::run_with_live_viewer(config, cli.seed, save_dir.as_deref(), 4);
    }

    #[cfg(not(feature = "vis"))]
    if cli.live {
        return Err(HermesError::Config(
            "live viewer requires --features vis".to_string(),
        ));
    }

    // Headless mode: run with observers.
    let mut observers: Vec<Box<dyn Observer>> = Vec::new();

    if let Some(ref dir) = save_dir {
        if !cli.quiet {
            println!("Saving snapshots to {dir}/");
        }
        observers.push(Box::new(FileObserver::new(dir)));
    }

    let total_steps = config.time.n_steps;
    let run_start = Instant::now();

    if !cli.quiet {
        use indicatif::{ProgressBar, ProgressStyle};

        let progress_bar = ProgressBar::new(total_steps as u64);
        progress_bar.set_style(
            ProgressStyle::with_template(
                "{spinner:.cyan} [{elapsed_precise}] [{bar:40.cyan/dark.grey}] step {pos}/{len} z={msg} ({eta} remaining)"
            )
            .unwrap()
            .progress_chars("=> "),
        );

        let n_steps = sim.run_with_progress(&mut observers, |step, scale_factor| {
            let redshift = 1.0 / scale_factor - 1.0;
            progress_bar.set_position(step as u64);
            progress_bar.set_message(format!("{redshift:.1}"));
        })?;

        let run_time = run_start.elapsed();
        progress_bar.finish_and_clear();
        println!(
            "Completed {n_steps} steps in {:.2}s ({:.1} ms/step)",
            run_time.as_secs_f64(),
            run_time.as_secs_f64() * 1000.0 / n_steps as f64,
        );
    } else {
        let n_steps = sim.run(&mut observers)?;
        let run_time = run_start.elapsed();
        let _ = (n_steps, run_time);
    }

    if !cli.quiet
        && let Some(diag) = sim.latest_diagnostics()
    {
        println!();
        println!("Final state (z = {:.1}):", 1.0 / sim.scale_factor - 1.0);
        println!("  Mass:     {:.4e} M_☉", diag.mass_total);
        println!("  |P|:      {:.4e}", diag.momentum_magnitude());
        println!("  E_kin:    {:.4e}", diag.energy_kinetic);
        println!("  E_pot:    {:.4e}", diag.energy_potential);
        println!("  |L|:      {:.4e}", diag.angular_momentum_magnitude());
    }

    Ok(())
}

/// Load and merge configuration from CLI arguments.
fn load_config(cli: &Cli) -> Result<Configuration, HermesError> {
    let file_override = if let Some(ref path) = cli.config_file {
        let content = std::fs::read_to_string(path)?;
        let value: toml::Value = toml::from_str(&content)
            .map_err(|e| HermesError::Config(format!("failed to parse {path}: {e}")))?;
        Some(value)
    } else {
        None
    };

    // Build programmatic overrides from CLI flags.
    let mut overrides = toml::map::Map::new();

    if cli.steps.is_some() || cli.particles.is_some() {
        if let Some(steps) = cli.steps {
            let mut time = toml::map::Map::new();
            time.insert("n_steps".to_string(), toml::Value::Integer(steps as i64));
            overrides.insert("time".to_string(), toml::Value::Table(time));
        }
        if let Some(particles) = cli.particles {
            let mut sim = toml::map::Map::new();
            sim.insert(
                "n_particles".to_string(),
                toml::Value::Integer(particles as i64),
            );
            // Also set n_cells to match particles for simplicity.
            sim.insert(
                "n_cells".to_string(),
                toml::Value::Integer(particles as i64),
            );
            overrides.insert("simulation".to_string(), toml::Value::Table(sim));
        }
    }

    let programmatic = if overrides.is_empty() {
        None
    } else {
        Some(toml::Value::Table(overrides))
    };

    build_configuration(file_override.as_ref(), programmatic.as_ref())
}

fn print_header(config: &Configuration, cli: &Cli) {
    println!("Hermes — {}", cli.scene);
    println!("{}", "=".repeat(40));
    println!("Grid:       {}³ cells", config.simulation.n_cells);
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
        1.0 / config.time.scale_factor_initial - 1.0,
        1.0 / config.time.scale_factor_final - 1.0,
    );
    println!("Steps:      {}", config.time.n_steps);
    println!("Seed:       {}", cli.seed);
    println!();
}

/// Play back saved snapshots.
fn run_playback(dir: &str, cli: &Cli) -> Result<(), HermesError> {
    #[cfg(not(feature = "vis"))]
    {
        let _ = (dir, cli);
        Err(HermesError::Config(
            "playback requires --features vis".to_string(),
        ))
    }

    #[cfg(feature = "vis")]
    {
        use kiss3d::light::Light;
        use kiss3d::nalgebra::Point3;
        use kiss3d::window::Window;

        use crate::io::snapshot::load_snapshot;
        use crate::vis::colormap::colormap_hot;

        // Load all snapshots.
        let mut snapshots = Vec::new();
        let mut step = 0_usize;
        loop {
            let path = format!("{dir}/snapshot_{step:05}.bin");
            match load_snapshot(std::path::Path::new(&path)) {
                Ok(snapshot) => {
                    snapshots.push(snapshot);
                    step += 1;
                }
                Err(_) => break,
            }
        }

        if snapshots.is_empty() {
            return Err(HermesError::Config(format!("no snapshots found in {dir}/")));
        }

        if !cli.quiet {
            println!("Loaded {} snapshots from {dir}/", snapshots.len());
            println!(
                "Scale factor: {:.4} → {:.4}",
                snapshots.first().unwrap().scale_factor,
                snapshots.last().unwrap().scale_factor,
            );
            println!("Particles: {}", snapshots[0].particle_count());
            println!("Precomputing display data...");
        }

        // Estimate box length from particle positions.
        let box_length = snapshots[0]
            .positions
            .iter()
            .flat_map(|pos| (0..3).map(move |d| pos.component(&[d]).abs()))
            .fold(0.0_f64, f64::max)
            * 1.1;
        let scale = 1.0 / box_length as f32;

        // Precompute frames.
        struct FrameData {
            positions: Vec<[f32; 3]>,
            colors: Vec<[f32; 3]>,
        }

        let frames: Vec<FrameData> = snapshots
            .iter()
            .map(|snapshot| {
                let speeds: Vec<f64> = snapshot.momenta.iter().map(|mom| mom.norm()).collect();
                let speed_max = speeds.iter().copied().fold(1e-30_f64, f64::max);
                let speed_min = speeds.iter().copied().fold(f64::MAX, f64::min);
                let speed_range = (speed_max - speed_min).max(1e-30);

                let mut positions = Vec::with_capacity(snapshot.particle_count());
                let mut colors = Vec::with_capacity(snapshot.particle_count());

                for (pos, &speed) in snapshot.positions.iter().zip(speeds.iter()) {
                    positions.push([
                        pos.component(&[0]) as f32 * scale - 0.5,
                        pos.component(&[1]) as f32 * scale - 0.5,
                        pos.component(&[2]) as f32 * scale - 0.5,
                    ]);
                    let normalized = ((speed - speed_min) / speed_range).clamp(0.0, 1.0);
                    colors.push(colormap_hot(normalized));
                }

                FrameData { positions, colors }
            })
            .collect();

        if !cli.quiet {
            println!("Playing at ~15 fps (close window to exit)");
        }

        let mut window = Window::new_with_size("hermes — playback", 1024, 768);
        window.set_background_color(0.0, 0.0, 0.0);
        window.set_light(Light::StickToCamera);
        window.set_point_size(3.0);

        let n_frames = frames.len();
        let mut frame_index = 0_usize;
        let frame_duration = std::time::Duration::from_millis(67);
        let mut last_frame_time = std::time::Instant::now();

        while window.render() {
            if last_frame_time.elapsed() >= frame_duration {
                frame_index = (frame_index + 1) % n_frames;
                last_frame_time = std::time::Instant::now();
            }

            let frame = &frames[frame_index];
            for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
                let point = Point3::new(pos[0], pos[1], pos[2]);
                let color_point = Point3::new(color[0], color[1], color[2]);
                window.draw_point(&point, &color_point);
            }
        }

        Ok(())
    }
}
