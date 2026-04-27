//! Run a simulation and save all snapshots to ./data/ for later playback.
//!
//! No visualization — just compute and save.
//!
//! Run with:
//!   cargo run --example run_and_save --release

use std::time::Instant;

use hermes_rs::config::build_configuration;
use hermes_rs::io::observer::{FileObserver, Observer};
use hermes_rs::physics::simulation::Simulation;

fn main() {
    let config = build_configuration(None, Some(&config())).expect("config");

    println!("Hermes — Save simulation snapshots");
    println!("===================================");
    println!("Grid:       {}³", config.simulation.n_cells);
    println!("Particles:  {}³", config.simulation.n_particles);
    println!("Steps:      {}", config.time.n_steps);
    println!();

    let start = Instant::now();
    let mut sim = Simulation::from_config(config, 42).expect("init failed");
    println!("Initialized in {:.2}s", start.elapsed().as_secs_f64());

    let file_observer = FileObserver::new("data");
    let mut observers: Vec<Box<dyn Observer>> = vec![Box::new(file_observer)];

    let run_start = Instant::now();
    let n_steps = sim.run(&mut observers).expect("run failed");
    println!(
        "Completed {n_steps} steps in {:.2}s",
        run_start.elapsed().as_secs_f64()
    );
}

fn config() -> toml::Value {
    toml::from_str(
        r#"
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
        snapshot_interval = 30
        "#,
    )
    .unwrap()
}
