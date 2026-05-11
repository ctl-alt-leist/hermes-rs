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

use crate::colormap::colormap_by_name;
use crate::config::VisualizationConfig;
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

/// How the viewer should render a frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RenderMode {
    /// Discrete points (particles).
    Points,
    /// Soft volumetric blobs with additive blending (fields).
    Volumetric,
}

/// Display-ready frame — flat f32 arrays, no morphis types.
pub struct DisplayFrame {
    pub positions: Vec<[f32; 3]>,
    pub colors: Vec<[f32; 3]>,
    pub render_mode: RenderMode,
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
///
/// If a progress bar is provided, it is incremented after each save
/// and finished when the writer receives Done.
pub fn spawn_disk_writer(
    rx: Receiver<PipelineMessage>,
    directory: String,
    progress: Option<indicatif::ProgressBar>,
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
                            Ok(()) => {
                                n_saved += 1;
                                if let Some(ref bar) = progress {
                                    bar.set_position(n_saved as u64);
                                }
                            }
                            Err(e) => eprintln!("warning: failed to save snapshot: {e}"),
                        }
                    }
                    PipelineMessage::Done => break,
                }
            }
            if let Some(ref bar) = progress {
                bar.finish_and_clear();
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
    vis_config: VisualizationConfig,
) -> thread::JoinHandle<()> {
    thread::Builder::new()
        .name("precompute".to_string())
        .spawn(move || {
            for msg in rx {
                match msg {
                    PipelineMessage::Snapshot(snapshot) => {
                        let frame = precompute_frame_rayon(&snapshot, box_length, &vis_config);
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
/// Renders all particle and field species into a single DisplayFrame.
/// Each species uses its own colormap from the visualization config.
pub fn precompute_frame_rayon(
    snapshot: &Snapshot,
    box_length: f64,
    vis: &VisualizationConfig,
) -> DisplayFrame {
    let scale = 1.0 / box_length as f32;

    let mut out_positions = Vec::new();
    let mut out_colors = Vec::new();

    // Render mode: volumetric if any fields present, points if particles only.
    let render_mode = if snapshot.has_fields() {
        RenderMode::Volumetric
    } else {
        RenderMode::Points
    };

    // ── Particle species ────────────────────────────────
    for species in &snapshot.particles {
        let species_config = vis.species.get(&species.name);
        let colormap_name = species_config.map(|c| c.colormap.as_str()).unwrap_or("hot");

        let speeds: Vec<f64> = species.momenta.iter().map(|mom| mom.norm()).collect();
        let speed_max = speeds.iter().copied().fold(1e-30_f64, f64::max);
        let speed_min = speeds.iter().copied().fold(f64::MAX, f64::min);
        let speed_range = (speed_max - speed_min).max(1e-30);

        for (pos, &speed) in species.positions.iter().zip(speeds.iter()) {
            out_positions.push([
                pos.component(&[1]) as f32 * scale - 0.5,
                pos.component(&[2]) as f32 * scale - 0.5,
                pos.component(&[3]) as f32 * scale - 0.5,
            ]);
            let normalized = ((speed - speed_min) / speed_range).clamp(0.0, 1.0);
            out_colors.push(colormap_by_name(colormap_name, normalized));
        }
    }

    // ── Field species ───────────────────────────────────
    let n = snapshot.n_cells;
    if n > 0 && !snapshot.fields.is_empty() {
        use rand::Rng;
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;

        let cell_size = 1.0 / n as f64;
        let jitter_amount = vis.jitter * cell_size;
        let total = n * n * n;

        // Precompute grid positions with jitter (shared across all species).
        let mut rng = ChaCha20Rng::seed_from_u64(0);
        let mut positions_grid: Vec<[f32; 3]> = Vec::with_capacity(total);
        for flat_idx in 0..total {
            let m0 = flat_idx / (n * n);
            let m1 = (flat_idx / n) % n;
            let m2 = flat_idx % n;

            let dx: f64 = rng.random_range(-jitter_amount..jitter_amount);
            let dy: f64 = rng.random_range(-jitter_amount..jitter_amount);
            let dz: f64 = rng.random_range(-jitter_amount..jitter_amount);

            let x = ((m0 as f64 + 0.5) * cell_size + dx - 0.5) as f32;
            let y = ((m1 as f64 + 0.5) * cell_size + dy - 0.5) as f32;
            let z = ((m2 as f64 + 0.5) * cell_size + dz - 0.5) as f32;

            positions_grid.push([x, y, z]);
        }

        for field_snapshot in &snapshot.fields {
            let species_config = vis.species.get(&field_snapshot.name);
            let colormap_name = species_config.map(|c| c.colormap.as_str()).unwrap_or("hot");
            let colormap_range = species_config
                .and_then(|c| c.colormap_range)
                .unwrap_or(vis.colormap_range);

            let density = &field_snapshot.density;
            let density_mean = density.iter().sum::<f64>() / density.len().max(1) as f64;
            let log_min = (colormap_range[0] * density_mean).max(1e-30).ln();
            let log_max = (colormap_range[1] * density_mean).ln();
            let log_range = (log_max - log_min).max(1e-10);

            for (flat_idx, &rho) in density.iter().enumerate() {
                out_positions.push(positions_grid[flat_idx]);

                let log_rho = rho.max(1e-30).ln();
                let normalized = ((log_rho - log_min) / log_range).clamp(0.0, 1.0);
                out_colors.push(colormap_by_name(colormap_name, normalized));
            }
        }
    }

    DisplayFrame {
        positions: out_positions,
        colors: out_colors,
        render_mode,
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
///
/// Uses kiss3d's State trait to inject the volumetric renderer for
/// field content alongside the standard point renderer for particles.
#[cfg(feature = "vis")]
pub fn run_viewer_main_thread(frame_rx: Receiver<ViewerMessage>, vis_config: VisualizationConfig) {
    use kiss3d::window::Window;

    let mut window = Window::new_with_size("hermes — live simulation", 1024, 768);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_point_size(vis_config.point_size);

    let mut state = ViewerState::new(frame_rx, &vis_config);

    while window.render_with_state(&mut state) {}
}

// ============================================================================
// Viewer state (kiss3d State trait)
// ============================================================================

/// Viewer state implementing kiss3d's State trait.
///
/// Manages the camera, volumetric renderer, and frame reception.
/// Dispatches on RenderMode: particle frames use draw_point,
/// field frames use the volumetric renderer with additive blending.
#[cfg(feature = "vis")]
struct ViewerState {
    camera: kiss3d::camera::ArcBall,
    volumetric: crate::visuals::volumetric_renderer::VolumetricRenderer,
    frame_rx: Receiver<ViewerMessage>,
    current_frame: Option<Box<DisplayFrame>>,
}

#[cfg(feature = "vis")]
impl ViewerState {
    fn new(frame_rx: Receiver<ViewerMessage>, vis: &VisualizationConfig) -> Self {
        use kiss3d::nalgebra::Point3;

        let e = vis.camera_eye();
        let eye = Point3::new(e[0], e[1], e[2]);
        let at = Point3::new(0.0, 0.0, 0.0);

        Self {
            camera: kiss3d::camera::ArcBall::new(eye, at),
            volumetric: crate::visuals::volumetric_renderer::VolumetricRenderer::new(
                vis.blob_size,
                vis.blob_alpha,
                vis.blob_falloff,
            ),
            frame_rx,
            current_frame: None,
        }
    }

    /// Create from preloaded frames (for playback).
    fn from_frames(
        frames: Vec<DisplayFrame>,
        fps: u64,
        vis: &VisualizationConfig,
    ) -> PlaybackState {
        use kiss3d::nalgebra::Point3;

        let e = vis.camera_eye();
        let eye = Point3::new(e[0], e[1], e[2]);
        let at = Point3::new(0.0, 0.0, 0.0);

        PlaybackState {
            camera: kiss3d::camera::ArcBall::new(eye, at),
            volumetric: crate::visuals::volumetric_renderer::VolumetricRenderer::new(
                vis.blob_size,
                vis.blob_alpha,
                vis.blob_falloff,
            ),
            frames,
            frame_index: 0,
            frame_duration: std::time::Duration::from_millis(1000 / fps.max(1)),
            last_frame_time: std::time::Instant::now(),
            paused: false,
        }
    }
}

#[cfg(feature = "vis")]
impl kiss3d::window::State for ViewerState {
    fn step(&mut self, window: &mut kiss3d::window::Window) {
        use kiss3d::light::Light;
        use kiss3d::nalgebra::Point3;

        window.set_light(Light::StickToCamera);

        // Drain the channel for the latest frame.
        while let Ok(msg) = self.frame_rx.try_recv() {
            match msg {
                ViewerMessage::Frame(frame) => self.current_frame = Some(frame),
                ViewerMessage::Done => {}
            }
        }

        if let Some(ref frame) = self.current_frame {
            match frame.render_mode {
                RenderMode::Points => {
                    for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
                        let point = Point3::new(pos[0], pos[1], pos[2]);
                        let color_point = Point3::new(color[0], color[1], color[2]);
                        window.draw_point(&point, &color_point);
                    }
                }
                RenderMode::Volumetric => {
                    for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
                        let point = Point3::new(pos[0], pos[1], pos[2]);
                        let color_point = Point3::new(color[0], color[1], color[2]);
                        self.volumetric.draw_point(point, color_point);
                    }
                }
            }
        }
    }

    fn cameras_and_effect_and_renderer(
        &mut self,
    ) -> (
        Option<&mut dyn kiss3d::camera::Camera>,
        Option<&mut dyn kiss3d::planar_camera::PlanarCamera>,
        Option<&mut dyn kiss3d::renderer::Renderer>,
        Option<&mut dyn kiss3d::post_processing::PostProcessingEffect>,
    ) {
        (
            Some(&mut self.camera),
            None,
            Some(&mut self.volumetric),
            None,
        )
    }
}

/// Playback state for pre-loaded frames.
///
/// Keyboard controls:
///   Space — play / pause
///   Left  — previous frame
///   Right — next frame
///   Down  — back 10%
///   Up    — forward 10%
///   Home  — jump to first frame
///   End   — jump to last frame
#[cfg(feature = "vis")]
struct PlaybackState {
    camera: kiss3d::camera::ArcBall,
    volumetric: crate::visuals::volumetric_renderer::VolumetricRenderer,
    frames: Vec<DisplayFrame>,
    frame_index: usize,
    frame_duration: std::time::Duration,
    last_frame_time: std::time::Instant,
    paused: bool,
}

#[cfg(feature = "vis")]
impl kiss3d::window::State for PlaybackState {
    fn step(&mut self, window: &mut kiss3d::window::Window) {
        use kiss3d::event::{Action, Key};
        use kiss3d::light::Light;
        use kiss3d::nalgebra::Point3;

        window.set_light(Light::StickToCamera);

        if self.frames.is_empty() {
            return;
        }

        // Keyboard controls.
        let n = self.frames.len();
        let jump = (n / 10).max(1);

        if window.get_key(Key::Space) == Action::Press {
            self.paused = !self.paused;
        }

        if window.get_key(Key::Left) == Action::Press && self.frame_index > 0 {
            self.frame_index -= 1;
            self.paused = true;
        }
        if window.get_key(Key::Right) == Action::Press && self.frame_index < n - 1 {
            self.frame_index += 1;
            self.paused = true;
        }
        if window.get_key(Key::Down) == Action::Press {
            self.frame_index = self.frame_index.saturating_sub(jump);
            self.paused = true;
        }
        if window.get_key(Key::Up) == Action::Press {
            self.frame_index = (self.frame_index + jump).min(n - 1);
            self.paused = true;
        }
        if window.get_key(Key::Home) == Action::Press {
            self.frame_index = 0;
            self.paused = true;
        }
        if window.get_key(Key::End) == Action::Press {
            self.frame_index = n - 1;
            self.paused = true;
        }

        if !self.paused && self.last_frame_time.elapsed() >= self.frame_duration {
            self.frame_index = (self.frame_index + 1) % n;
            self.last_frame_time = std::time::Instant::now();
        }

        let frame = &self.frames[self.frame_index];
        match frame.render_mode {
            RenderMode::Points => {
                for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
                    let point = Point3::new(pos[0], pos[1], pos[2]);
                    let color_point = Point3::new(color[0], color[1], color[2]);
                    window.draw_point(&point, &color_point);
                }
            }
            RenderMode::Volumetric => {
                for (pos, color) in frame.positions.iter().zip(frame.colors.iter()) {
                    let point = Point3::new(pos[0], pos[1], pos[2]);
                    let color_point = Point3::new(color[0], color[1], color[2]);
                    self.volumetric.draw_point(point, color_point);
                }
            }
        }
    }

    fn cameras_and_effect_and_renderer(
        &mut self,
    ) -> (
        Option<&mut dyn kiss3d::camera::Camera>,
        Option<&mut dyn kiss3d::planar_camera::PlanarCamera>,
        Option<&mut dyn kiss3d::renderer::Renderer>,
        Option<&mut dyn kiss3d::post_processing::PostProcessingEffect>,
    ) {
        (
            Some(&mut self.camera),
            None,
            Some(&mut self.volumetric),
            None,
        )
    }
}

/// Run the playback viewer on the main thread with a loader thread.
///
/// Loads snapshots from disk one at a time on a background thread,
/// precomputes DisplayFrames, and sends them to the viewer. Memory
/// usage is bounded by the channel capacity.
#[cfg(feature = "vis")]
pub fn run_playback_viewer(
    dir: &str,
    fps: u64,
    vis_config: &VisualizationConfig,
) -> Result<(), crate::error::HermesError> {
    use std::sync::mpsc as playback_mpsc;

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
    let vis = vis_config.clone();
    let (frame_tx, frame_rx) = playback_mpsc::sync_channel::<ViewerMessage>(32);

    let loader_handle = thread::Builder::new()
        .name("playback-loader".to_string())
        .spawn(move || {
            for path in &snapshot_paths {
                match crate::io::snapshot::load_snapshot(path) {
                    Ok(snapshot) => {
                        let frame = precompute_frame_rayon(&snapshot, box_length, &vis);
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
    let mut frames: Vec<DisplayFrame> = Vec::with_capacity(total);

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
                    frames.push(*frame);
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
    window.set_point_size(vis_config.point_size);

    let mut state = ViewerState::from_frames(frames, fps, vis_config);

    while window.render_with_state(&mut state) {}

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
    if snapshot.has_fields() {
        snapshot.n_cells as f64
    } else {
        snapshot
            .particles
            .iter()
            .flat_map(|s| s.positions.iter())
            .flat_map(|pos: &morphis::vector::Vector<3>| {
                (1..=3).map(move |d| pos.component(&[d]).abs())
            })
            .fold(0.0_f64, f64::max)
            * 1.1
    }
}
