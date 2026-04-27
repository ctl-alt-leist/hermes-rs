//! Profile where time is spent in the live viewer pipeline.
//!
//! Run with:
//!   cargo run --example profile_viewer --features vis --release

use std::time::Instant;

use hermes_rs::io::snapshot::Snapshot;
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::initial::zeldovich_init;
use hermes_rs::physics::poisson::PoissonSolver;

fn main() {
    let grid = Grid::new(32, 50_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);
    let particles = zeldovich_init(32, &grid, &cosmology, 0.02, 42).unwrap();

    // 1. Time snapshot capture (morphis Vec<Vector<3>> construction).
    let start = Instant::now();
    let mut snapshots = Vec::new();
    for _ in 0..10 {
        let snapshot = Snapshot::capture(&particles, &grid, &cosmology, &mut solver, 0, 0.02);
        snapshots.push(snapshot);
    }
    let snapshot_time = start.elapsed();
    println!(
        "Snapshot capture: {:.1} ms each ({} particles, {} morphis vectors)",
        snapshot_time.as_secs_f64() * 100.0,
        particles.count(),
        particles.count() * 2
    );

    // 2. Time snapshot clone (what the channel observer does).
    let snapshot = &snapshots[0];
    let start = Instant::now();
    for _ in 0..100 {
        let _cloned = snapshot.clone();
    }
    let clone_time = start.elapsed();
    println!(
        "Snapshot clone:   {:.3} ms each",
        clone_time.as_secs_f64() * 10.0
    );

    // 3. Time the density estimation (spatial hash binning in update_display_data).
    let n = snapshot.particle_count();
    let scale = 1.0 / 50_000.0_f32;

    let start = Instant::now();
    for _ in 0..100 {
        let mut density_estimate = vec![1.0_f64; n];
        let n_bins = 16_usize;
        let bin_size = 1.0 / n_bins as f64;
        let mut bin_counts = vec![0_usize; n_bins * n_bins * n_bins];

        for pos in snapshot.positions.iter() {
            let x = (pos.component(&[0]) * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let y = (pos.component(&[1]) * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let z = (pos.component(&[2]) * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let bx = (x / bin_size) as usize;
            let by = (y / bin_size) as usize;
            let bz = (z / bin_size) as usize;
            bin_counts[bx * n_bins * n_bins + by * n_bins + bz] += 1;
        }

        let mean_count = n as f64 / (n_bins * n_bins * n_bins) as f64;

        for (p, pos) in snapshot.positions.iter().enumerate() {
            let x = (pos.component(&[0]) * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let y = (pos.component(&[1]) * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let z = (pos.component(&[2]) * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let bx = (x / bin_size) as usize;
            let by = (y / bin_size) as usize;
            let bz = (z / bin_size) as usize;
            let count = bin_counts[bx * n_bins * n_bins + by * n_bins + bz] as f64;
            density_estimate[p] = count / mean_count;
        }
        std::hint::black_box(&density_estimate);
    }
    let density_time = start.elapsed();
    println!(
        "Density estimate: {:.3} ms each",
        density_time.as_secs_f64() * 10.0
    );

    // 4. Time morphis component extraction (the hot path in update_display_data).
    let start = Instant::now();
    for _ in 0..100 {
        let mut positions = vec![[0.0_f32; 3]; n];
        for (p, pos) in snapshot.positions.iter().enumerate() {
            positions[p] = [
                pos.component(&[0]) as f32,
                pos.component(&[1]) as f32,
                pos.component(&[2]) as f32,
            ];
        }
        std::hint::black_box(&positions);
    }
    let extract_time = start.elapsed();
    println!(
        "Component extract: {:.3} ms each ({n} × 3 morphis .component() calls)",
        extract_time.as_secs_f64() * 10.0
    );

    // 5. For comparison: time flat array extraction.
    let start = Instant::now();
    for _ in 0..100 {
        let mut positions = vec![[0.0_f32; 3]; n];
        for p in 0..n {
            positions[p] = [
                particles.position[[0, p]] as f32,
                particles.position[[1, p]] as f32,
                particles.position[[2, p]] as f32,
            ];
        }
        std::hint::black_box(&positions);
    }
    let flat_time = start.elapsed();
    println!(
        "Flat array extract: {:.3} ms each ({n} × 3 direct array accesses)",
        flat_time.as_secs_f64() * 10.0
    );

    println!();
    println!(
        "Morphis overhead: {:.1}x vs flat arrays",
        extract_time.as_secs_f64() / flat_time.as_secs_f64()
    );

    // 6. Time the full Snapshot::capture including diagnostics (Poisson solve).
    let start = Instant::now();
    let _snapshot = Snapshot::capture(&particles, &grid, &cosmology, &mut solver, 0, 0.02);
    let full_capture_time = start.elapsed();
    println!();
    println!(
        "Full Snapshot::capture (incl. diagnostics + Poisson): {:.1} ms",
        full_capture_time.as_secs_f64() * 1000.0
    );
    println!("  → This runs every snapshot_interval steps during simulation");
}
