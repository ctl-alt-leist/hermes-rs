//! Cosmological N-body simulation with live visualization.
//!
//! Runs a 32³ particle-mesh simulation from z ≈ 49 to z = 0 with:
//!   - Live 3D viewer updating as the simulation runs
//!   - Snapshots saved to ./data/ for post-hoc analysis
//!   - Post-run plots: density slice, power spectrum, conservation
//!
//! Run with:
//!   cargo run --example cosmic_web --features vis --release

use std::path::Path;
use std::time::Instant;

use hermes_rs::config::build_configuration;
use hermes_rs::io::observer::{FileObserver, Observer};
use hermes_rs::physics::simulation::Simulation;
use hermes_rs::vis;
use hermes_rs::vis::LiveObserver;

fn main() {
    let config = build_configuration(None, Some(&small_universe())).expect("config");

    println!("Hermes — Cosmological PM Simulation");
    println!("====================================");
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
    println!();

    let box_length = config.simulation.box_length;

    // Initialize.
    let start = Instant::now();
    let mut sim = Simulation::from_config(config, 42).expect("initialization failed");

    let init_time = start.elapsed();
    println!("Initialized in {:.2}s", init_time.as_secs_f64());

    let initial_diag = sim.diagnostics_history[0].clone();
    println!("  Mass:     {:.4e} M_☉", initial_diag.mass_total);
    println!("  |P|:      {:.4e}", initial_diag.momentum_magnitude());
    println!("  E_kin:    {:.4e}", initial_diag.energy_kinetic);
    println!("  E_pot:    {:.4e}", initial_diag.energy_potential);
    println!();

    // Set up observers: live 3D viewer + file snapshots.
    let live_observer = LiveObserver::new(box_length, 4);
    let file_observer = FileObserver::new("data");
    let mut observers: Vec<Box<dyn Observer>> =
        vec![Box::new(live_observer), Box::new(file_observer)];

    // Run with live visualization.
    println!("Running simulation with live viewer...");
    let run_start = Instant::now();
    let n_steps = sim.run(&mut observers).expect("simulation failed");
    let run_time = run_start.elapsed();

    println!(
        "Completed {n_steps} steps in {:.2}s ({:.1} ms/step)",
        run_time.as_secs_f64(),
        run_time.as_secs_f64() * 1000.0 / n_steps as f64,
    );

    let final_diag = sim.latest_diagnostics().unwrap();
    println!();
    println!("Final state (z = {:.1}):", 1.0 / sim.scale_factor - 1.0);
    println!(
        "  Mass:     {:.4e} M_☉ (Δ = {:.2e})",
        final_diag.mass_total,
        (final_diag.mass_total - initial_diag.mass_total) / initial_diag.mass_total,
    );
    println!("  |P|:      {:.4e}", final_diag.momentum_magnitude());
    println!("  E_kin:    {:.4e}", final_diag.energy_kinetic);
    println!("  E_pot:    {:.4e}", final_diag.energy_potential);
    println!(
        "  |L|:      {:.4e}",
        final_diag.angular_momentum_magnitude()
    );

    // Post-run plots.
    println!();
    println!("Generating plots...");

    let output = Path::new("output");
    std::fs::create_dir_all(output).ok();

    vis::render_density_slice(
        &sim.particles,
        &sim.grid,
        sim.grid.box_length / 2.0,
        sim.grid.box_length / 4.0,
        &output.join("density_slice.png"),
        512,
    )
    .expect("density slice failed");
    println!("  output/density_slice.png");

    vis::plot_power_spectrum(
        &sim.particles,
        &sim.grid,
        &output.join("power_spectrum.png"),
    )
    .expect("power spectrum failed");
    println!("  output/power_spectrum.png");

    vis::plot_conservation(&sim.diagnostics_history, &output.join("conservation.png"))
        .expect("conservation plot failed");
    println!("  output/conservation.png");
}

/// Small universe config: 32³ grid/particles, 50 Mpc box, 100 steps.
fn small_universe() -> toml::Value {
    toml::from_str(
        r#"
        [cosmology]
        hubble       = 0.674
        omega_m      = 0.315
        omega_b      = 0.0493
        omega_r      = 9.15e-5
        omega_k      = 0.0
        omega_lambda = 0.6849085
        sigma_8      = 0.811
        spectral_index = 0.965

        [simulation]
        n_cells      = 32
        n_particles  = 32
        box_length   = 50000.0

        [time]
        scale_factor_initial = 0.02
        scale_factor_final   = 1.0
        n_steps              = 100
        stepping             = "log_a"

        [output]
        directory         = "output"
        snapshot_interval = 10
        "#,
    )
    .unwrap()
}
