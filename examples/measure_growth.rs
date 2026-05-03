//! Measure delta_rms and delta_max from saved snapshots.
//!
//! Usage:
//!   cargo run --example measure_growth --release -- data/cosmic-web-field-growing

use hermes_rs::io::snapshot::{SnapshotContent, load_snapshot};
use hermes_rs::run::pipeline::find_snapshot_paths;

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Usage: measure_growth <snapshot_dir>");
        std::process::exit(1);
    });

    let paths = find_snapshot_paths(&dir);
    if paths.is_empty() {
        eprintln!("No snapshots found in {dir}/");
        return;
    }

    println!(
        "{:>6}  {:>6}  {:>12}  {:>12}  {:>12}  {:>12}",
        "step", "a", "delta_rms", "delta_max", "delta_min", "rho_mean"
    );

    // Sample every ~10% of snapshots to keep output compact.
    let stride = (paths.len() / 20).max(1);

    for (n, path) in paths.iter().enumerate() {
        if n % stride != 0 && n != paths.len() - 1 {
            continue;
        }

        let snapshot = match load_snapshot(path) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("  skip {}: {e}", path.display());
                continue;
            }
        };

        if let SnapshotContent::Fields { ref density, .. } = snapshot.content {
            let rho_mean = density.iter().sum::<f64>() / density.len() as f64;

            let deltas: Vec<f64> = density.iter().map(|&rho| rho / rho_mean - 1.0).collect();
            let delta_rms =
                (deltas.iter().map(|d| d * d).sum::<f64>() / deltas.len() as f64).sqrt();
            let delta_max = deltas.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let delta_min = deltas.iter().copied().fold(f64::INFINITY, f64::min);

            println!(
                "{:6}  {:6.4}  {:12.6e}  {:12.6e}  {:12.6e}  {:12.6e}",
                snapshot.step, snapshot.scale_factor, delta_rms, delta_max, delta_min, rho_mean
            );
        }
    }
}
