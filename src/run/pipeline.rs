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
    /// If true, frames are dropped when the channel is full (viewer).
    /// If false, the router blocks until the consumer catches up (disk writer).
    droppable: bool,
}

impl Consumer {
    fn send(&self, msg: PipelineMessage) {
        if self.droppable {
            let _ = self.tx.try_send(msg);
        } else {
            let _ = self.tx.send(msg);
        }
    }

    fn send_blocking(&self, msg: PipelineMessage) {
        let _ = self.tx.send(msg);
    }
}

// ============================================================================
// Router
// ============================================================================

/// Configuration for a consumer channel.
pub struct ConsumerConfig {
    pub tx: SyncSender<PipelineMessage>,
    /// If true, frames are dropped when channel is full (viewer).
    /// If false, the router blocks until the consumer catches up (disk).
    pub droppable: bool,
}

/// Spawn the router thread: receives from simulation, fans out to consumers.
///
/// Snapshots are Arc::clone'd to each consumer (reference count only).
/// Non-droppable consumers (disk writer) receive every frame; droppable
/// consumers (viewer) have frames silently dropped when their channel is full.
pub fn spawn_router(
    rx: Receiver<PipelineMessage>,
    consumers: Vec<ConsumerConfig>,
) -> thread::JoinHandle<()> {
    let consumers: Vec<Consumer> = consumers
        .into_iter()
        .map(|c| Consumer {
            tx: c.tx,
            droppable: c.droppable,
        })
        .collect();

    thread::Builder::new()
        .name("pipeline-router".to_string())
        .spawn(move || {
            for msg in rx {
                match msg {
                    PipelineMessage::Snapshot(snapshot) => {
                        for consumer in &consumers {
                            consumer.send(PipelineMessage::Snapshot(Arc::clone(&snapshot)));
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

/// Convert a snapshot to a display-ready frame.
///
/// Dispatches on content type: particles use velocity colormap,
/// fields sample density at grid points with brightness colormap.
pub fn precompute_frame_rayon(snapshot: &Snapshot, box_length: f64) -> DisplayFrame {
    use crate::io::snapshot::SnapshotContent;

    let scale = 1.0 / box_length as f32;

    match &snapshot.content {
        SnapshotContent::Particles {
            positions, momenta, ..
        } => {
            let speeds: Vec<f64> = momenta.iter().map(|mom| mom.norm()).collect();
            let speed_max = speeds.iter().copied().fold(1e-30_f64, f64::max);
            let speed_min = speeds.iter().copied().fold(f64::MAX, f64::min);
            let speed_range = (speed_max - speed_min).max(1e-30);

            let mut out_positions = Vec::with_capacity(positions.len());
            let mut out_colors = Vec::with_capacity(positions.len());

            for (pos, &speed) in positions.iter().zip(speeds.iter()) {
                out_positions.push([
                    pos.component(&[0]) as f32 * scale - 0.5,
                    pos.component(&[1]) as f32 * scale - 0.5,
                    pos.component(&[2]) as f32 * scale - 0.5,
                ]);
                let normalized = ((speed - speed_min) / speed_range).clamp(0.0, 1.0);
                out_colors.push(colormap_hot(normalized));
            }

            DisplayFrame {
                positions: out_positions,
                colors: out_colors,
                step: snapshot.step,
                scale_factor: snapshot.scale_factor,
            }
        }
        SnapshotContent::Fields { density, n_cells } => {
            let n = *n_cells;
            if n == 0 || density.is_empty() {
                return DisplayFrame {
                    positions: Vec::new(),
                    colors: Vec::new(),
                    step: snapshot.step,
                    scale_factor: snapshot.scale_factor,
                };
            }

            // Density-weighted random sampling: draw points with probability
            // proportional to density. Dense regions get many points,
            // empty regions get none. Positions randomized within each cell
            // to avoid grid artifacts.
            use rand::Rng;
            use rand::SeedableRng;
            use rand_chacha::ChaCha20Rng;

            let cell_size = 1.0 / n as f64;
            let total_density: f64 = density.iter().sum();
            let n_sample_points = 30_000; // target point count

            let density_max = density.iter().copied().fold(0.0_f64, f64::max);
            let density_min = density
                .iter()
                .copied()
                .filter(|&d| d > 0.0)
                .fold(f64::MAX, f64::min);
            let log_min = density_min.max(1e-30).ln();
            let log_max = density_max.max(1e-30).ln();
            let log_range = (log_max - log_min).max(1e-10);

            let mut out_positions = Vec::with_capacity(n_sample_points);
            let mut out_colors = Vec::with_capacity(n_sample_points);

            // Use step as seed so each frame has different sampling noise.
            let mut rng = ChaCha20Rng::seed_from_u64(snapshot.step as u64);

            for _ in 0..n_sample_points {
                // Pick a random cell weighted by density.
                let target = rng.random_range(0.0..total_density);
                let mut cumulative = 0.0;
                let mut chosen = 0_usize;

                for (flat_idx, &rho) in density.iter().enumerate() {
                    cumulative += rho;
                    if cumulative >= target {
                        chosen = flat_idx;
                        break;
                    }
                }

                let m0 = chosen / (n * n);
                let m1 = (chosen / n) % n;
                let m2 = chosen % n;

                // Random position within the cell.
                let jitter_x: f64 = rng.random_range(0.0..1.0);
                let jitter_y: f64 = rng.random_range(0.0..1.0);
                let jitter_z: f64 = rng.random_range(0.0..1.0);

                let x = ((m0 as f64 + jitter_x) * cell_size - 0.5) as f32;
                let y = ((m1 as f64 + jitter_y) * cell_size - 0.5) as f32;
                let z = ((m2 as f64 + jitter_z) * cell_size - 0.5) as f32;

                out_positions.push([x, y, z]);

                let rho = density[chosen];
                let log_rho = rho.max(1e-30).ln();
                let normalized = ((log_rho - log_min) / log_range).clamp(0.0, 1.0);
                let brightness = normalized * normalized;
                out_colors.push(colormap_hot(brightness));
            }

            DisplayFrame {
                positions: out_positions,
                colors: out_colors,
                step: snapshot.step,
                scale_factor: snapshot.scale_factor,
            }
        }
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
    use kiss3d::camera::ArcBall;
    use kiss3d::light::Light;
    use kiss3d::nalgebra::Point3;
    use kiss3d::window::Window;

    let mut window = Window::new_with_size("hermes — live simulation", 1024, 768);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(5.0);

    // Camera at a 3/4 angle, zoomed out to see the full box.
    let eye = Point3::new(0.8, 0.6, 1.0);
    let at = Point3::new(0.0, 0.0, 0.0);
    let mut camera = ArcBall::new(eye, at);

    let mut current_frame: Option<Box<DisplayFrame>> = None;

    while window.render_with_camera(&mut camera) {
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

    let snapshot_paths = find_snapshot_paths(dir);
    let total = snapshot_paths.len();
    if total == 0 {
        return Err(crate::error::HermesError::Config(format!(
            "no snapshots found in {dir}/"
        )));
    }

    println!("Loading {total} snapshots from {dir}/...");

    // Estimate box length from first snapshot.
    let first = load_snapshot(&snapshot_paths[0])?;
    let box_length = estimate_box_length(&first);

    // Loader thread: load + precompute frames, send via channel.
    let (frame_tx, frame_rx) = playback_mpsc::sync_channel::<ViewerMessage>(32);

    let loader_handle = thread::Builder::new()
        .name("playback-loader".to_string())
        .spawn(move || {
            for path in &snapshot_paths {
                match crate::io::snapshot::load_snapshot(path) {
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
                    Err(e) => eprintln!("warning: failed to load {}: {e}", path.display()),
                }
            }
            let _ = frame_tx.send(ViewerMessage::Done);
        })
        .expect("failed to spawn playback loader");

    // Collect all frames from the loader thread before starting playback.
    // This ensures smooth rendering — no I/O during the render loop.
    let mut frames: Vec<Box<DisplayFrame>> = Vec::with_capacity(total);

    {
        use indicatif::{ProgressBar, ProgressStyle};

        let progress = ProgressBar::new(total as u64);
        progress.set_style(
            ProgressStyle::with_template(
                "{spinner:.cyan} [{bar:40.cyan/dark.grey}] {pos}/{len} frames loaded",
            )
            .unwrap()
            .progress_chars("=> "),
        );

        for msg in frame_rx {
            match msg {
                ViewerMessage::Frame(frame) => {
                    frames.push(frame);
                    progress.set_position(frames.len() as u64);
                }
                ViewerMessage::Done => break,
            }
        }

        progress.finish_and_clear();
    }

    let _ = loader_handle.join();

    println!(
        "Playing {} frames at ~{fps} fps (close window to exit)",
        frames.len()
    );

    // Render loop — pure playback from precomputed data, no I/O.
    let mut window = Window::new_with_size("hermes — playback", 1024, 768);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(5.0);

    // Camera at a 3/4 angle, zoomed out to see the full box.
    let eye = Point3::new(0.8, 0.6, 1.0);
    let at = Point3::new(0.0, 0.0, 0.0);
    let mut camera = kiss3d::camera::ArcBall::new(eye, at);

    let n_frames = frames.len();
    let frame_duration = std::time::Duration::from_millis(1000 / fps.max(1));
    let mut last_frame_time = std::time::Instant::now();
    let mut frame_index = 0_usize;

    while window.render_with_camera(&mut camera) {
        if last_frame_time.elapsed() >= frame_duration {
            frame_index = (frame_index + 1) % n_frames;
            last_frame_time = std::time::Instant::now();
        }

        let frame = &frames[frame_index];
        for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
            let point = Point3::new(pos[0], pos[1], pos[2]);
            let color_point = Point3::new(color[0], color[1], color[2]);
            window.draw_point(&point, &color_point);
        }
    }

    Ok(())
}

// ============================================================================
// Utilities
// ============================================================================

/// Find all snapshot files in a directory, sorted by name.
///
/// Scans for all `snapshot-*.bin` files rather than assuming sequential
/// numbering (the pipeline may drop frames under load).
pub fn find_snapshot_paths(dir: &str) -> Vec<std::path::PathBuf> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };

    let mut paths: Vec<std::path::PathBuf> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("snapshot-") && n.ends_with(".bin"))
        })
        .collect();

    paths.sort();

    paths
}

/// Count snapshot files in a directory.
pub fn count_snapshots(dir: &str) -> usize {
    find_snapshot_paths(dir).len()
}

#[cfg(feature = "vis")]
/// Estimate the box length from a snapshot.
///
/// For particles: maximum coordinate extent.
/// For fields: uses n_cells as a proxy (assumes unit box, scale = 1).
fn estimate_box_length(snapshot: &crate::io::snapshot::Snapshot) -> f64 {
    use crate::io::snapshot::SnapshotContent;

    match &snapshot.content {
        SnapshotContent::Particles { positions, .. } => {
            positions
                .iter()
                .flat_map(|pos: &morphis::vector::Vector<3>| {
                    (0..3).map(move |d| pos.component(&[d]).abs())
                })
                .fold(0.0_f64, f64::max)
                * 1.1
        }
        SnapshotContent::Fields { n_cells, .. } => {
            // Field snapshots don't carry the box length explicitly.
            // Use 1.0 as the normalized box — the precompute maps to [-0.5, 0.5].
            *n_cells as f64
        }
    }
}
