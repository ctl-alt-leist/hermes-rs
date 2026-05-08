//! Inspect particle positions and momenta from saved snapshots.
//!
//! Usage: cargo run --example inspect_particles --release -- data/cosmic-web-pm

use hermes_rs::io::snapshot::load_snapshot;
use hermes_rs::run::pipeline::find_snapshot_paths;

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Usage: inspect_particles <snapshot_dir>");
        std::process::exit(1);
    });

    let paths = find_snapshot_paths(&dir);
    let stride = (paths.len() / 10).max(1);

    println!(
        "{:>6}  {:>6}  {:>12}  {:>12}  {:>12}  {:>12}",
        "step", "a", "max_disp_kpc", "disp_%_space", "mean_|p|", "max_|p|"
    );

    for (n, path) in paths.iter().enumerate() {
        if n % stride != 0 && n != paths.len() - 1 {
            continue;
        }

        let snap = load_snapshot(path).unwrap();
        if let Some(species) = snap.particles.first() {
            let positions = &species.positions;
            let momenta = &species.momenta;
            let n_p = (positions.len() as f64).cbrt().round() as usize;
            let spacing = 100000.0 / n_p as f64; // assumes 100 Mpc box

            let mut max_disp = 0.0_f64;
            for (p, pos) in positions.iter().enumerate() {
                let m0 = p / (n_p * n_p);
                let m1 = (p / n_p) % n_p;
                let m2 = p % n_p;
                let x0 = (m0 as f64 + 0.5) * spacing;
                let y0 = (m1 as f64 + 0.5) * spacing;
                let z0 = (m2 as f64 + 0.5) * spacing;

                // Periodic distance
                let wrap = |d: f64, l: f64| -> f64 {
                    let d = d.abs();
                    if d > l / 2.0 { l - d } else { d }
                };
                let dx = wrap(pos.component(&[1]) - x0, 100000.0);
                let dy = wrap(pos.component(&[2]) - y0, 100000.0);
                let dz = wrap(pos.component(&[3]) - z0, 100000.0);
                let disp = (dx * dx + dy * dy + dz * dz).sqrt();
                if disp > max_disp {
                    max_disp = disp;
                }
            }

            let max_mom = momenta.iter().map(|m| m.norm()).fold(0.0_f64, f64::max);
            let mean_mom = momenta.iter().map(|m| m.norm()).sum::<f64>() / momenta.len() as f64;

            println!(
                "{:6}  {:6.4}  {:12.1}  {:12.2}  {:12.4e}  {:12.4e}",
                snap.step,
                snap.scale_factor,
                max_disp,
                max_disp / spacing * 100.0,
                mean_mom,
                max_mom
            );
        }
    }
}
