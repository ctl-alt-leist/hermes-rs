//! Cosmological N-body simulation with live 3D visualization.
//!
//! Runs a 32³ particle-mesh simulation from z ≈ 49 to z = 0 with:
//!   - Live 3D viewer updating in real time as the simulation runs
//!   - Snapshots saved to ./data/ for post-hoc analysis
//!
//! Run with:
//!   cargo run --example cosmic_web --features vis --release

use hermes_rs::config::build_configuration;
use hermes_rs::vis;

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
    println!("Starting live viewer + simulation...");
    println!("(close the viewer window to exit)");
    println!();

    // Viewer on main thread, simulation on background thread.
    // Snapshots saved to ./data/ automatically.
    vis::run_with_live_viewer(config, 42, Some("data"), 4).expect("simulation failed");
}

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
        n_steps              = 300
        stepping             = "log_a"

        [output]
        directory         = "output"
        snapshot_interval = 30
        "#,
    )
    .unwrap()
}
