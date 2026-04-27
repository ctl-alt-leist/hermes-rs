//! Live visualization with the viewer on the main thread.
//!
//! On macOS, the window event loop must run on the main thread. The
//! simulation runs on a spawned background thread and sends precomputed
//! display frames through a bounded channel.
//!
//! The simulation thread converts each snapshot to a `DisplayFrame`
//! (flat `[f32; 3]` arrays for positions and colors) before sending.
//! The render thread just draws — no morphis extraction, no density
//! computation, no per-frame allocation.

use std::sync::mpsc;
use std::thread;

use kiss3d::light::Light;
use kiss3d::nalgebra::Point3;
use kiss3d::window::Window;

use crate::config::Configuration;
use crate::io::observer::{FileObserver, Observer};
use crate::io::snapshot::Snapshot;
use crate::physics::simulation::Simulation;
use crate::vis::colormap::colormap_hot;

/// Precomputed display data — everything the renderer needs, no morphis.
struct DisplayFrame {
    positions: Vec<[f32; 3]>,
    colors: Vec<[f32; 3]>,
}

enum LiveMessage {
    Frame(Box<DisplayFrame>),
    Done,
}

/// Run a simulation with a live 3D viewer.
///
/// The viewer runs on the calling (main) thread. The simulation runs on
/// a spawned background thread, precomputing display frames and sending
/// them through a bounded channel. Close the window to exit.
pub fn run_with_live_viewer(
    config: Configuration,
    seed: u64,
    snapshot_dir: Option<&str>,
    buffer_size: usize,
) -> Result<(), crate::error::HermesError> {
    let box_length = config.simulation.box_length;
    let snapshot_dir_owned = snapshot_dir.map(|s| s.to_string());
    let (sender, receiver) = mpsc::sync_channel::<LiveMessage>(buffer_size);

    let sim_handle = thread::spawn(move || -> Result<Simulation, crate::error::HermesError> {
        let mut sim = Simulation::from_config(config, seed)?;

        let channel_observer = ChannelObserver {
            sender: sender.clone(),
            box_length,
        };
        let mut observers: Vec<Box<dyn Observer>> = vec![Box::new(channel_observer)];

        if let Some(dir) = snapshot_dir_owned {
            observers.push(Box::new(FileObserver::new(dir)));
        }

        sim.run(&mut observers)?;

        let _ = sender.send(LiveMessage::Done);

        Ok(sim)
    });

    render_loop(receiver);

    match sim_handle.join() {
        Ok(Ok(sim)) => {
            println!(
                "Simulation complete: {} steps, a = {:.4}",
                sim.step, sim.scale_factor
            );
            Ok(())
        }
        Ok(Err(e)) => Err(e),
        Err(_) => Err(crate::error::HermesError::Config(
            "simulation thread panicked".to_string(),
        )),
    }
}

// ============================================================================
// Channel observer — precomputes display data on the simulation thread
// ============================================================================

struct ChannelObserver {
    sender: mpsc::SyncSender<LiveMessage>,
    box_length: f64,
}

impl Observer for ChannelObserver {
    fn on_snapshot(&mut self, snapshot: &Snapshot) {
        let frame = precompute_frame(snapshot, self.box_length);
        let message = LiveMessage::Frame(Box::new(frame));
        match self.sender.try_send(message) {
            Ok(()) => {}
            Err(mpsc::TrySendError::Full(_)) => {}
            Err(mpsc::TrySendError::Disconnected(_)) => {}
        }
    }

    fn on_complete(&mut self) {}
}

/// Convert a snapshot to display-ready arrays on the simulation thread.
fn precompute_frame(snapshot: &Snapshot, box_length: f64) -> DisplayFrame {
    let n = snapshot.particle_count();
    let scale = 1.0 / box_length as f32;

    // Color by velocity magnitude (momentum norm).
    let speeds: Vec<f64> = snapshot.momenta.iter().map(|mom| mom.norm()).collect();

    let speed_max = speeds.iter().copied().fold(1e-30_f64, f64::max);
    let speed_min = speeds.iter().copied().fold(f64::MAX, f64::min);
    let speed_range = (speed_max - speed_min).max(1e-30);

    let mut positions = Vec::with_capacity(n);
    let mut colors = Vec::with_capacity(n);

    for (pos, &speed) in snapshot.positions.iter().zip(speeds.iter()) {
        positions.push([
            pos.component(&[0]) as f32 * scale - 0.5,
            pos.component(&[1]) as f32 * scale - 0.5,
            pos.component(&[2]) as f32 * scale - 0.5,
        ]);

        let normalized = ((speed - speed_min) / speed_range).clamp(0.0, 1.0);
        colors.push(colormap_hot(normalized));
    }

    DisplayFrame { positions, colors }
}

// ============================================================================
// Render loop (main thread) — just draws precomputed data
// ============================================================================

fn render_loop(receiver: mpsc::Receiver<LiveMessage>) {
    let mut window = Window::new_with_size("hermes — live simulation", 1024, 768);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(3.0);

    let mut current_positions: Vec<[f32; 3]> = Vec::new();
    let mut current_colors: Vec<[f32; 3]> = Vec::new();

    while window.render() {
        // Drain channel — use the latest frame.
        while let Ok(message) = receiver.try_recv() {
            match message {
                LiveMessage::Frame(frame) => {
                    current_positions = frame.positions;
                    current_colors = frame.colors;
                }
                LiveMessage::Done => {}
            }
        }

        for (pos, color) in current_positions.iter().zip(current_colors.iter()) {
            let point = Point3::new(pos[0], pos[1], pos[2]);
            let color_point = Point3::new(color[0], color[1], color[2]);
            window.draw_point(&point, &color_point);
        }
    }
}
