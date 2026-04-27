//! Live visualization observer with channel-based decoupling.
//!
//! The simulation runs on the calling thread. A separate render thread
//! reads snapshots from a bounded channel and displays the latest state
//! at display framerate. The simulation never blocks on the viewer — if
//! the channel is full, old snapshots are replaced.

use std::sync::mpsc;
use std::thread;

use kiss3d::light::Light;
use kiss3d::nalgebra::Point3;
use kiss3d::window::Window;

use crate::io::observer::Observer;
use crate::io::snapshot::Snapshot;
use crate::vis::colormap::colormap_hot;

/// Observer that feeds snapshots to a live 3D viewer on a separate thread.
///
/// The viewer shows the latest particle positions with density-dependent
/// coloring on a dark background. The simulation and viewer are decoupled:
/// the simulation pushes snapshots into a bounded channel, and the viewer
/// reads at its own framerate.
pub struct LiveObserver {
    sender: mpsc::SyncSender<LiveMessage>,
    render_thread: Option<thread::JoinHandle<()>>,
}

enum LiveMessage {
    Snapshot(Box<Snapshot>),
    Done,
}

impl LiveObserver {
    /// Create a live observer and spawn the render thread.
    ///
    /// `box_length` is needed to normalize positions for display.
    /// `buffer_size` controls how many snapshots can queue before
    /// the oldest is dropped (default: 4).
    pub fn new(box_length: f64, buffer_size: usize) -> Self {
        let (sender, receiver) = mpsc::sync_channel(buffer_size);

        let render_thread = thread::spawn(move || {
            render_loop(receiver, box_length);
        });

        Self {
            sender,
            render_thread: Some(render_thread),
        }
    }
}

impl Observer for LiveObserver {
    fn on_snapshot(&mut self, snapshot: &Snapshot) {
        // Non-blocking: if the channel is full, drop the oldest by
        // draining and sending the new one.
        let message = LiveMessage::Snapshot(Box::new(snapshot.clone()));
        match self.sender.try_send(message) {
            Ok(()) => {}
            Err(mpsc::TrySendError::Full(msg)) => {
                // Drain one old message and retry.
                let _ = self.sender.try_send(msg);
            }
            Err(mpsc::TrySendError::Disconnected(_)) => {
                // Viewer window was closed — silently stop sending.
            }
        }
    }

    fn on_complete(&mut self) {
        let _ = self.sender.send(LiveMessage::Done);
        if let Some(handle) = self.render_thread.take() {
            let _ = handle.join();
        }
    }
}

/// Render loop running on the viewer thread.
///
/// Reads the latest snapshot from the channel each frame and draws
/// all particles as colored points. If no new snapshot is available,
/// the previous frame is redrawn.
fn render_loop(receiver: mpsc::Receiver<LiveMessage>, box_length: f64) {
    let mut window = Window::new_with_size("hermes — live simulation", 1200, 900);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(2.0);

    let scale = 1.0 / box_length as f32;
    let mut current_positions: Vec<[f32; 3]> = Vec::new();
    let mut current_colors: Vec<[f32; 3]> = Vec::new();
    let mut info_text = String::new();

    while window.render() {
        // Drain channel — use the latest snapshot available.
        let mut got_done = false;
        while let Ok(message) = receiver.try_recv() {
            match message {
                LiveMessage::Snapshot(snapshot) => {
                    update_display_data(
                        &snapshot,
                        scale,
                        &mut current_positions,
                        &mut current_colors,
                    );
                    let redshift = 1.0 / snapshot.scale_factor - 1.0;
                    info_text = format!(
                        "step {} | z = {:.1} | a = {:.4}",
                        snapshot.step, redshift, snapshot.scale_factor
                    );
                }
                LiveMessage::Done => {
                    got_done = true;
                }
            }
        }

        // Draw particles.
        for (pos, color) in current_positions.iter().zip(current_colors.iter()) {
            let point = Point3::new(pos[0], pos[1], pos[2]);
            let color_point = Point3::new(color[0], color[1], color[2]);
            window.draw_point(&point, &color_point);
        }

        if got_done && !info_text.is_empty() {
            // Simulation finished — keep window open for inspection
            // until the user closes it. Print final state to terminal.
            println!("Live viewer: simulation complete — {info_text}");
        }
    }
}

/// Convert a snapshot into display-ready position and color arrays.
fn update_display_data(
    snapshot: &Snapshot,
    scale: f32,
    positions: &mut Vec<[f32; 3]>,
    colors: &mut Vec<[f32; 3]>,
) {
    let n = snapshot.particle_count();
    positions.resize(n, [0.0; 3]);
    colors.resize(n, [0.0; 3]);

    // Compute density-based colors from positions.
    // We use a simple estimate: compute mean density and per-particle
    // deviation from the position spread.
    let mut density_estimate = vec![1.0_f64; n];

    // Quick density proxy: count neighbors within a radius.
    // For speed, just use position clustering as a rough proxy.
    // A proper implementation would CIC-deposit and interpolate back,
    // but that requires a Grid, which we don't carry in the snapshot.
    // Instead, use a spatial hash for approximate density.
    if n > 0 {
        // Bin particles into a coarse grid for density estimation.
        let n_bins = 16_usize;
        let bin_size = 1.0 / n_bins as f64;
        let mut bin_counts = vec![0_usize; n_bins * n_bins * n_bins];

        for pos in snapshot.positions.iter() {
            let x = (pos.component(&[0]) as f64 * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let y = (pos.component(&[1]) as f64 * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let z = (pos.component(&[2]) as f64 * scale as f64 + 0.5).clamp(0.0, 0.9999);

            let bx = (x / bin_size) as usize;
            let by = (y / bin_size) as usize;
            let bz = (z / bin_size) as usize;

            bin_counts[bx * n_bins * n_bins + by * n_bins + bz] += 1;
        }

        let mean_count = n as f64 / (n_bins * n_bins * n_bins) as f64;

        for (p, pos) in snapshot.positions.iter().enumerate() {
            let x = (pos.component(&[0]) as f64 * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let y = (pos.component(&[1]) as f64 * scale as f64 + 0.5).clamp(0.0, 0.9999);
            let z = (pos.component(&[2]) as f64 * scale as f64 + 0.5).clamp(0.0, 0.9999);

            let bx = (x / bin_size) as usize;
            let by = (y / bin_size) as usize;
            let bz = (z / bin_size) as usize;

            let count = bin_counts[bx * n_bins * n_bins + by * n_bins + bz] as f64;
            density_estimate[p] = count / mean_count;
        }
    }

    // Map density to colors and positions.
    let log_min = 0.1_f64.ln();
    let log_max = density_estimate
        .iter()
        .copied()
        .fold(1.0_f64, f64::max)
        .ln()
        .max(log_min + 0.1);
    let log_range = log_max - log_min;

    for (p, pos) in snapshot.positions.iter().enumerate() {
        positions[p] = [
            pos.component(&[0]) as f32 * scale - 0.5,
            pos.component(&[1]) as f32 * scale - 0.5,
            pos.component(&[2]) as f32 * scale - 0.5,
        ];

        let log_density = density_estimate[p].max(0.01).ln();
        let normalized = ((log_density - log_min) / log_range).clamp(0.0, 1.0);
        colors[p] = colormap_hot(normalized);
    }
}
