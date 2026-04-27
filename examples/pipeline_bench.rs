//! Compare old observer-based path vs new pipeline path.
//!
//! Run with:
//!   cargo run --example pipeline_bench --release

use std::sync::Arc;
use std::sync::mpsc;
use std::time::Instant;

use hermes_rs::config::build_configuration;
use hermes_rs::io::observer::{FileObserver, NullObserver, Observer};
use hermes_rs::io::snapshot::Snapshot;
use hermes_rs::physics::simulation::Simulation;
use hermes_rs::run::pipeline::{PipelineMessage, SnapshotSender, spawn_disk_writer, spawn_router};

fn config(particles: usize, steps: usize) -> hermes_rs::config::Configuration {
    let overrides: toml::Value = toml::from_str(&format!(
        r#"
        [simulation]
        n_cells = {particles}
        n_particles = {particles}
        [time]
        n_steps = {steps}
        scale_factor_initial = 0.02
        scale_factor_final = 1.0
        [output]
        snapshot_interval = 10
        "#
    ))
    .unwrap();

    build_configuration(None, Some(&overrides)).unwrap()
}

fn main() {
    println!("Pipeline vs Observer Benchmark");
    println!("==============================\n");

    for &n in &[16, 32] {
        let steps = 50;
        println!("--- {}³ particles, {steps} steps ---\n", n);

        // 1. Observer path: no observers (baseline simulation speed).
        let cfg = config(n, steps);
        let mut sim = Simulation::from_config(cfg, 42).unwrap();
        let start = Instant::now();
        sim.run(&mut []).unwrap();
        let baseline = start.elapsed();
        println!(
            "  Baseline (no observers):     {:>7.1} ms  ({:.1} ms/step)",
            baseline.as_secs_f64() * 1000.0,
            baseline.as_secs_f64() * 1000.0 / steps as f64
        );

        // 2. Observer path: NullObserver (measures snapshot capture overhead).
        let cfg = config(n, steps);
        let mut sim = Simulation::from_config(cfg, 42).unwrap();
        let mut observers: Vec<Box<dyn Observer>> = vec![Box::new(NullObserver)];
        let start = Instant::now();
        sim.run(&mut observers).unwrap();
        let null_observer = start.elapsed();
        let capture_overhead = null_observer.as_secs_f64() - baseline.as_secs_f64();
        println!(
            "  Observer (NullObserver):      {:>7.1} ms  (capture overhead: {:.1} ms)",
            null_observer.as_secs_f64() * 1000.0,
            capture_overhead * 1000.0,
        );

        // 3. Observer path: FileObserver (measures capture + sync write).
        let cfg = config(n, steps);
        let mut sim = Simulation::from_config(cfg, 42).unwrap();
        let dir = format!("data/bench-observer-{n}");
        let mut observers: Vec<Box<dyn Observer>> = vec![Box::new(FileObserver::new(&dir))];
        let start = Instant::now();
        sim.run(&mut observers).unwrap();
        let file_observer = start.elapsed();
        println!(
            "  Observer (FileObserver):      {:>7.1} ms  (sync write overhead: {:.1} ms)",
            file_observer.as_secs_f64() * 1000.0,
            (file_observer.as_secs_f64() - baseline.as_secs_f64()) * 1000.0,
        );

        // 4. Pipeline path: no consumers (measures Arc<Snapshot> overhead).
        let cfg = config(n, steps);
        let mut sim = Simulation::from_config(cfg, 42).unwrap();
        let (sim_tx, router_rx) = mpsc::sync_channel::<PipelineMessage>(8);
        let sender = SnapshotSender::new(sim_tx);
        let router_handle = spawn_router(router_rx, vec![]);
        let start = Instant::now();
        sim.run_into_pipeline(&sender, |_, _| {}).unwrap();
        let pipeline_empty = start.elapsed();
        let _ = router_handle.join();
        println!(
            "  Pipeline (no consumers):     {:>7.1} ms  (Arc overhead: {:.1} ms)",
            pipeline_empty.as_secs_f64() * 1000.0,
            (pipeline_empty.as_secs_f64() - baseline.as_secs_f64()) * 1000.0,
        );

        // 5. Pipeline path: disk writer (measures async write).
        let cfg = config(n, steps);
        let mut sim = Simulation::from_config(cfg, 42).unwrap();
        let (sim_tx, router_rx) = mpsc::sync_channel::<PipelineMessage>(8);
        let sender = SnapshotSender::new(sim_tx);
        let dir = format!("data/bench-pipeline-{n}");
        let (disk_tx, disk_rx) = mpsc::sync_channel::<PipelineMessage>(16);
        let router_handle = spawn_router(router_rx, vec![disk_tx]);
        let disk_handle = spawn_disk_writer(disk_rx, dir);
        let start = Instant::now();
        sim.run_into_pipeline(&sender, |_, _| {}).unwrap();
        let pipeline_disk = start.elapsed();
        let _ = router_handle.join();
        let _ = disk_handle.join();
        println!(
            "  Pipeline (disk writer):      {:>7.1} ms  (sim not blocked by writes)",
            pipeline_disk.as_secs_f64() * 1000.0,
        );

        // Summary.
        let speedup = file_observer.as_secs_f64() / pipeline_disk.as_secs_f64();
        println!();
        println!(
            "  Speedup (pipeline+disk vs observer+disk): {:.2}x",
            speedup
        );
        println!();

        // Clean up bench data.
        let _ = std::fs::remove_dir_all(format!("data/bench-observer-{n}"));
        let _ = std::fs::remove_dir_all(format!("data/bench-pipeline-{n}"));
    }
}
