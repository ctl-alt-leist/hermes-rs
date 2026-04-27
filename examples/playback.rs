//! Play back saved snapshots from ./data/ in a 3D viewer.
//!
//! Loads all snapshot files, then plays them at a controlled framerate.
//! The viewer is not bound to simulation speed — it renders at 60fps
//! and advances through snapshots at a configurable pace.
//!
//! Run with:
//!   cargo run --example playback --features vis --release

use std::time::Instant;

use kiss3d::light::Light;
use kiss3d::nalgebra::Point3;
use kiss3d::window::Window;

use hermes_rs::io::snapshot::{Snapshot, load_snapshot};
use hermes_rs::vis::colormap::colormap_hot;

fn main() {
    // Load all snapshots from ./data/.
    let mut snapshots = Vec::new();
    let mut step = 0_usize;

    loop {
        let path = format!("data/snapshot_{step:05}.bin");
        match load_snapshot(std::path::Path::new(&path)) {
            Ok(snapshot) => {
                snapshots.push(snapshot);
                step += 1;
            }
            Err(_) => break,
        }
    }

    if snapshots.is_empty() {
        eprintln!("No snapshots found in ./data/");
        eprintln!("Run `cargo run --example run_and_save --release` first.");
        std::process::exit(1);
    }

    println!("Loaded {} snapshots", snapshots.len());
    println!(
        "Scale factor: {:.4} → {:.4}",
        snapshots.first().unwrap().scale_factor,
        snapshots.last().unwrap().scale_factor,
    );
    println!("Particles: {}", snapshots[0].particle_count());

    // Precompute display data for all frames.
    let box_length = estimate_box_length(&snapshots[0]);
    let scale = 1.0 / box_length as f32;

    println!("Precomputing display data...");
    let frames: Vec<FrameData> = snapshots
        .iter()
        .map(|snapshot| precompute_frame(snapshot, scale))
        .collect();
    println!("Ready. Playing at ~15 fps (close window to exit)");
    println!();

    // Playback loop.
    let mut window = Window::new_with_size("hermes — playback", 1024, 768);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(3.0);

    let n_frames = frames.len();
    let mut frame_index = 0_usize;
    let frame_duration = std::time::Duration::from_millis(67); // ~15 fps
    let mut last_frame_time = Instant::now();
    let looping = true;

    while window.render() {
        // Advance frame at controlled rate.
        if last_frame_time.elapsed() >= frame_duration {
            frame_index += 1;
            if frame_index >= n_frames {
                if looping {
                    frame_index = 0;
                } else {
                    frame_index = n_frames - 1;
                }
            }
            last_frame_time = Instant::now();
        }

        let frame = &frames[frame_index];

        // Draw all particles.
        for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
            let point = Point3::new(pos[0], pos[1], pos[2]);
            let color_point = Point3::new(color[0], color[1], color[2]);
            window.draw_point(&point, &color_point);
        }
    }
}

struct FrameData {
    positions: Vec<[f32; 3]>,
    colors: Vec<[f32; 3]>,
}

fn precompute_frame(snapshot: &Snapshot, scale: f32) -> FrameData {
    let n = snapshot.particle_count();

    // Color by velocity magnitude (momentum norm as proxy).
    // Fast particles in dense regions glow bright; slow ones in voids stay dark.
    let speeds: Vec<f64> = snapshot.momenta.iter().map(|mom| mom.norm()).collect();

    let speed_max = speeds.iter().copied().fold(1e-30_f64, f64::max);
    let speed_min = speeds.iter().copied().fold(f64::MAX, f64::min);
    let speed_range = (speed_max - speed_min).max(1e-30);

    let mut positions = Vec::with_capacity(n);
    let mut colors = Vec::with_capacity(n);

    for p in 0..n {
        positions.push([
            snapshot.positions[p].component(&[0]) as f32 * scale - 0.5,
            snapshot.positions[p].component(&[1]) as f32 * scale - 0.5,
            snapshot.positions[p].component(&[2]) as f32 * scale - 0.5,
        ]);

        let normalized = ((speeds[p] - speed_min) / speed_range).clamp(0.0, 1.0);
        colors.push(colormap_hot(normalized));
    }

    FrameData { positions, colors }
}

fn estimate_box_length(snapshot: &Snapshot) -> f64 {
    snapshot
        .positions
        .iter()
        .flat_map(|pos| (0..3).map(move |d| pos.component(&[d]).abs()))
        .fold(0.0_f64, f64::max)
        * 1.1
}
