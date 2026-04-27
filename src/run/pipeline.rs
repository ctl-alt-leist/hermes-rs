//! Pipeline threading architecture.
//!
//! The simulation, disk I/O, and visualization run on independent threads
//! connected by bounded channels. The simulation produces `Arc<Snapshot>`
//! into a channel; a router fans out to consumers (disk writer, viewer
//! precompute); the main thread owns the event loop.
//!
//! ```text
//! Sim Thread ──ch(8)──> Router ──ch(16)──> Disk Writer
//!                          └────ch(4)──> Precompute ──ch(4)──> Main (Viewer)
//! ```

use std::sync::Arc;
use std::sync::mpsc::{Receiver, SyncSender};
use std::thread;

use rayon::prelude::*;

use crate::colormap::colormap_hot;
use crate::io::snapshot::{Snapshot, save_snapshot};

// ============================================================================
// Message types
// ============================================================================

/// Message flowing through the simulation → router pipeline.
pub enum PipelineMessage {
    /// A simulation snapshot, Arc-wrapped for zero-copy fan-out.
    Snapshot(Arc<Snapshot>),
    /// Simulation has finished.
    Done,
}

/// Message flowing from precompute → viewer (main thread).
pub enum ViewerMessage {
    /// A display-ready frame.
    Frame(Box<DisplayFrame>),
    /// Simulation has finished. Viewer stays open for inspection.
    Done,
}

/// Display-ready frame — flat f32 arrays, no morphis types.
pub struct DisplayFrame {
    pub positions: Vec<[f32; 3]>,
    pub colors: Vec<[f32; 3]>,
    pub step: usize,
    pub scale_factor: f64,
}

// ============================================================================
// Simulation sender
// ============================================================================

/// Channel-based output for the simulation thread.
///
/// Replaces the Observer trait in the pipeline architecture. The
/// simulation sends `Arc<Snapshot>` into the channel; the router
/// handles fan-out to consumers.
pub struct SnapshotSender {
    tx: SyncSender<PipelineMessage>,
}

impl SnapshotSender {
    pub fn new(tx: SyncSender<PipelineMessage>) -> Self {
        Self { tx }
    }

    /// Send a snapshot into the pipeline. Non-blocking: drops if full.
    pub fn send(&self, snapshot: Arc<Snapshot>) {
        let _ = self.tx.try_send(PipelineMessage::Snapshot(snapshot));
    }

    /// Signal that the simulation is complete.
    pub fn done(&self) {
        let _ = self.tx.send(PipelineMessage::Done);
    }
}

// ============================================================================
// Consumer
// ============================================================================

/// A downstream consumer connected to the router via a channel.
struct Consumer {
    tx: SyncSender<PipelineMessage>,
}

impl Consumer {
    fn try_send(&self, msg: PipelineMessage) {
        let _ = self.tx.try_send(msg);
    }

    fn send_blocking(&self, msg: PipelineMessage) {
        let _ = self.tx.send(msg);
    }
}

// ============================================================================
// Router
// ============================================================================

/// Spawn the router thread: receives from simulation, fans out to consumers.
///
/// Snapshots are Arc::clone'd to each consumer (reference count only).
/// Uses try_send for snapshots (drop if consumer is slow) and blocking
/// send for Done (ensure every consumer shuts down).
pub fn spawn_router(
    rx: Receiver<PipelineMessage>,
    consumers: Vec<SyncSender<PipelineMessage>>,
) -> thread::JoinHandle<()> {
    let consumers: Vec<Consumer> = consumers.into_iter().map(|tx| Consumer { tx }).collect();

    thread::Builder::new()
        .name("pipeline-router".to_string())
        .spawn(move || {
            for msg in rx {
                match msg {
                    PipelineMessage::Snapshot(snapshot) => {
                        for consumer in &consumers {
                            consumer.try_send(PipelineMessage::Snapshot(Arc::clone(&snapshot)));
                        }
                    }
                    PipelineMessage::Done => {
                        for consumer in &consumers {
                            consumer.send_blocking(PipelineMessage::Done);
                        }
                        break;
                    }
                }
            }
        })
        .expect("failed to spawn router thread")
}

// ============================================================================
// Disk writer
// ============================================================================

/// Spawn a disk writer thread. Receives snapshots and saves to bincode files.
pub fn spawn_disk_writer(
    rx: Receiver<PipelineMessage>,
    directory: String,
) -> thread::JoinHandle<()> {
    thread::Builder::new()
        .name("disk-writer".to_string())
        .spawn(move || {
            let mut n_saved = 0_usize;
            for msg in rx {
                match msg {
                    PipelineMessage::Snapshot(snapshot) => {
                        let filename = format!("snapshot-{:05}.bin", snapshot.step);
                        let path = std::path::PathBuf::from(&directory).join(filename);
                        match save_snapshot(&snapshot, &path) {
                            Ok(()) => n_saved += 1,
                            Err(e) => eprintln!("warning: failed to save snapshot: {e}"),
                        }
                    }
                    PipelineMessage::Done => break,
                }
            }
            if n_saved > 0 {
                println!("DiskWriter: saved {n_saved} snapshots to {directory}");
            }
        })
        .expect("failed to spawn disk writer thread")
}

// ============================================================================
// Precompute (snapshot → display frame)
// ============================================================================

/// Spawn a precompute thread: converts Arc<Snapshot> to DisplayFrame
/// using rayon-parallel operations, then sends to the viewer.
#[cfg(feature = "vis")]
pub fn spawn_precompute(
    rx: Receiver<PipelineMessage>,
    frame_tx: SyncSender<ViewerMessage>,
    box_length: f64,
) -> thread::JoinHandle<()> {
    thread::Builder::new()
        .name("precompute".to_string())
        .spawn(move || {
            for msg in rx {
                match msg {
                    PipelineMessage::Snapshot(snapshot) => {
                        let frame = precompute_frame_rayon(&snapshot, box_length);
                        let _ = frame_tx.try_send(ViewerMessage::Frame(Box::new(frame)));
                    }
                    PipelineMessage::Done => {
                        let _ = frame_tx.send(ViewerMessage::Done);
                        break;
                    }
                }
            }
        })
        .expect("failed to spawn precompute thread")
}

/// Convert a snapshot to a display-ready frame using rayon parallelism.
pub fn precompute_frame_rayon(snapshot: &Snapshot, box_length: f64) -> DisplayFrame {
    let scale = 1.0 / box_length as f32;

    let speeds: Vec<f64> = snapshot.momenta.par_iter().map(|mom| mom.norm()).collect();

    let speed_max = speeds.par_iter().copied().reduce(|| 1e-30_f64, f64::max);
    let speed_min = speeds.par_iter().copied().reduce(|| f64::MAX, f64::min);
    let speed_range = (speed_max - speed_min).max(1e-30);

    let (positions, colors): (Vec<_>, Vec<_>) = snapshot
        .positions
        .par_iter()
        .zip(speeds.par_iter())
        .map(|(pos, &speed)| {
            let p = [
                pos.component(&[0]) as f32 * scale - 0.5,
                pos.component(&[1]) as f32 * scale - 0.5,
                pos.component(&[2]) as f32 * scale - 0.5,
            ];
            let normalized = ((speed - speed_min) / speed_range).clamp(0.0, 1.0);
            let c = colormap_hot(normalized);
            (p, c)
        })
        .unzip();

    DisplayFrame {
        positions,
        colors,
        step: snapshot.step,
        scale_factor: snapshot.scale_factor,
    }
}

// ============================================================================
// Viewer (main thread)
// ============================================================================

/// Run the viewer event loop on the main thread.
///
/// Receives DisplayFrames from the precompute thread via channel.
/// Draws the latest frame each render tick. Stays open after simulation
/// completes. Close the window to exit.
#[cfg(feature = "vis")]
pub fn run_viewer_main_thread(frame_rx: Receiver<ViewerMessage>) {
    use kiss3d::light::Light;
    use kiss3d::nalgebra::Point3;
    use kiss3d::window::Window;

    let mut window = Window::new_with_size("hermes — live simulation", 1024, 768);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(3.0);

    let mut current_frame: Option<Box<DisplayFrame>> = None;

    while window.render() {
        while let Ok(msg) = frame_rx.try_recv() {
            match msg {
                ViewerMessage::Frame(frame) => current_frame = Some(frame),
                ViewerMessage::Done => {}
            }
        }

        if let Some(ref frame) = current_frame {
            for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
                let point = Point3::new(pos[0], pos[1], pos[2]);
                let color_point = Point3::new(color[0], color[1], color[2]);
                window.draw_point(&point, &color_point);
            }
        }
    }
}

/// Run the playback viewer on the main thread with a loader thread.
///
/// Loads snapshots from disk one at a time on a background thread,
/// precomputes DisplayFrames, and sends them to the viewer. Memory
/// usage is bounded by the channel capacity.
#[cfg(feature = "vis")]
pub fn run_playback_viewer(dir: &str, fps: u64) -> Result<(), crate::error::HermesError> {
    use std::sync::mpsc as playback_mpsc;

    use kiss3d::light::Light;
    use kiss3d::nalgebra::Point3;
    use kiss3d::window::Window;

    use crate::io::snapshot::load_snapshot;

    let total = count_snapshots(dir);
    if total == 0 {
        return Err(crate::error::HermesError::Config(format!(
            "no snapshots found in {dir}/"
        )));
    }

    println!("Loading {total} snapshots from {dir}/...");

    // Estimate box length from first snapshot.
    let first = load_snapshot(std::path::Path::new(&format!("{dir}/snapshot-00000.bin")))?;
    let box_length = first
        .positions
        .iter()
        .flat_map(|pos: &morphis::vector::Vector<3>| (0..3).map(move |d| pos.component(&[d]).abs()))
        .fold(0.0_f64, f64::max)
        * 1.1;

    // Loader thread: load + precompute frames, send via channel.
    let (frame_tx, frame_rx) = playback_mpsc::sync_channel::<ViewerMessage>(32);
    let dir_owned = dir.to_string();

    let loader_handle = thread::Builder::new()
        .name("playback-loader".to_string())
        .spawn(move || {
            // Precompute all frames, sending each as ready.
            // Snapshots are discarded after precompute — only DisplayFrames
            // stay in the channel buffer.
            for step in 0..total {
                let path = format!("{dir_owned}/snapshot-{step:05}.bin");
                match crate::io::snapshot::load_snapshot(std::path::Path::new(&path)) {
                    Ok(snapshot) => {
                        let frame = precompute_frame_rayon(&snapshot, box_length);
                        // Blocking send — playback should not drop frames.
                        if frame_tx
                            .send(ViewerMessage::Frame(Box::new(frame)))
                            .is_err()
                        {
                            break; // Viewer closed.
                        }
                    }
                    Err(e) => eprintln!("warning: failed to load {path}: {e}"),
                }
            }
            let _ = frame_tx.send(ViewerMessage::Done);
        })
        .expect("failed to spawn playback loader");

    // Collect all frames for looping playback.
    // The loader sends them as they're ready; we buffer them for replay.
    let mut frames: Vec<Box<DisplayFrame>> = Vec::with_capacity(total);
    let mut loading_done = false;

    let mut window = Window::new_with_size("hermes — playback", 1024, 768);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(3.0);

    let frame_duration = std::time::Duration::from_millis(1000 / fps.max(1));
    let mut last_frame_time = std::time::Instant::now();
    let mut frame_index = 0_usize;

    while window.render() {
        // Drain any newly loaded frames.
        while let Ok(msg) = frame_rx.try_recv() {
            match msg {
                ViewerMessage::Frame(frame) => frames.push(frame),
                ViewerMessage::Done => loading_done = true,
            }
        }

        // Advance frame at controlled rate.
        if !frames.is_empty() && last_frame_time.elapsed() >= frame_duration {
            frame_index += 1;
            if frame_index >= frames.len() {
                if loading_done {
                    frame_index = 0; // Loop.
                } else {
                    frame_index = frames.len() - 1; // Hold last while loading.
                }
            }
            last_frame_time = std::time::Instant::now();
        }

        // Draw current frame.
        if let Some(frame) = frames.get(frame_index) {
            for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
                let point = Point3::new(pos[0], pos[1], pos[2]);
                let color_point = Point3::new(color[0], color[1], color[2]);
                window.draw_point(&point, &color_point);
            }
        }
    }

    let _ = loader_handle.join();

    Ok(())
}

// ============================================================================
// Utilities
// ============================================================================

/// Count snapshot files in a directory.
pub fn count_snapshots(dir: &str) -> usize {
    let mut count = 0;
    loop {
        let path = format!("{dir}/snapshot-{count:05}.bin");
        if !std::path::Path::new(&path).exists() {
            break;
        }
        count += 1;
    }

    count
}
