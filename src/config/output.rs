/// Output configuration: what we record, report, and display.
///
/// Parsed from the `[output]` section of the TOML config.
use std::collections::BTreeMap;

use serde::Deserialize;

// ============================================================================
// Output block
// ============================================================================

/// What we record, report, and display.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct OutputBlock {
    #[serde(default)]
    pub snapshots: SnapshotsConfig,
    #[serde(default)]
    pub diagnostics: DiagnosticsConfig,
    #[serde(default)]
    pub logging: LoggingConfig,
    #[serde(default)]
    pub display: DisplayConfig,
}

/// Snapshot output configuration.
#[derive(Debug, Clone, Deserialize)]
pub struct SnapshotsConfig {
    /// Save a snapshot every N steps.
    #[serde(default = "default_snapshot_interval")]
    pub interval: usize,
    /// Output directory. Empty string defaults to data/<scene>/.
    #[serde(default)]
    pub directory: String,
}

fn default_snapshot_interval() -> usize {
    1
}

impl Default for SnapshotsConfig {
    fn default() -> Self {
        Self {
            interval: default_snapshot_interval(),
            directory: String::new(),
        }
    }
}

/// Diagnostics configuration.
#[derive(Debug, Clone, Deserialize)]
pub struct DiagnosticsConfig {
    /// Compute conservation audits every N steps.
    #[serde(default = "default_diagnostics_interval")]
    pub interval: usize,
}

fn default_diagnostics_interval() -> usize {
    10
}

impl Default for DiagnosticsConfig {
    fn default() -> Self {
        Self {
            interval: default_diagnostics_interval(),
        }
    }
}

/// Terminal logging configuration.
#[derive(Debug, Clone, Deserialize)]
pub struct LoggingConfig {
    /// Log level: "error", "warn", "info", "debug", "trace".
    #[serde(default = "default_log_level")]
    pub level: String,
    /// Print step progress to terminal.
    #[serde(default = "default_progress")]
    pub progress: bool,
}

fn default_log_level() -> String {
    "info".to_string()
}
fn default_progress() -> bool {
    true
}

impl Default for LoggingConfig {
    fn default() -> Self {
        Self {
            level: default_log_level(),
            progress: default_progress(),
        }
    }
}

/// Visualization / display configuration.
#[derive(Debug, Clone, Deserialize)]
pub struct DisplayConfig {
    /// Screen-space point size for particle rendering.
    #[serde(default = "default_point_size")]
    pub point_size: f32,
    /// Screen-space blob size for volumetric field rendering.
    #[serde(default = "default_blob_size")]
    pub blob_size: f32,
    /// Per-blob opacity for additive blending.
    #[serde(default = "default_blob_alpha")]
    pub blob_alpha: f32,
    /// Gaussian falloff rate for volumetric blobs.
    #[serde(default = "default_blob_falloff")]
    pub blob_falloff: f32,
    /// Camera distance from origin (box spans [-0.5, 0.5]).
    #[serde(default = "default_camera_distance")]
    pub camera_distance: f32,
    /// Camera direction as [x, y, z] (multiplied by distance).
    #[serde(default = "default_camera_angle")]
    pub camera_angle: [f32; 3],
    /// Colormap range as [floor, ceiling] in units of density / rho_mean.
    #[serde(default = "default_colormap_range")]
    pub colormap_range: [f64; 2],
    /// Grid-point jitter as fraction of cell size.
    #[serde(default = "default_jitter")]
    pub jitter: f64,
    /// Pixel resolution for GIF recording.
    #[serde(default = "default_gif_resolution")]
    pub gif_resolution: u32,
    /// Point radius in pixels for GIF rendering.
    #[serde(default = "default_gif_point_radius")]
    pub gif_point_radius: i32,
    /// Per-species display overrides, keyed by species name.
    ///
    /// Each entry specifies the colormap and optional range for a
    /// field species. Species not listed here use "hot" with the
    /// global colormap_range.
    #[serde(default)]
    pub species: BTreeMap<String, SpeciesDisplayConfig>,
}

/// Per-species visualization settings.
#[derive(Debug, Clone, Deserialize)]
pub struct SpeciesDisplayConfig {
    /// Colormap name: "hot", "cool", "ember", "verdant".
    #[serde(default = "default_colormap")]
    pub colormap: String,
    /// Colormap range as [floor, ceiling] in density / rho_mean.
    /// If absent, falls back to the global colormap_range.
    pub colormap_range: Option<[f64; 2]>,
    /// Per-species blob size override.
    pub blob_size: Option<f32>,
    /// Per-species blob opacity override.
    pub blob_alpha: Option<f32>,
}

fn default_colormap() -> String {
    "hot".to_string()
}

fn default_point_size() -> f32 {
    5.0
}
fn default_blob_size() -> f32 {
    18.0
}
fn default_blob_alpha() -> f32 {
    0.12
}
fn default_blob_falloff() -> f32 {
    10.0
}
fn default_camera_distance() -> f32 {
    1.9
}
fn default_camera_angle() -> [f32; 3] {
    [0.56, 0.42, 0.69]
}
fn default_colormap_range() -> [f64; 2] {
    [0.3, 3.0]
}
fn default_jitter() -> f64 {
    0.3
}
fn default_gif_resolution() -> u32 {
    512
}
fn default_gif_point_radius() -> i32 {
    1
}

impl Default for DisplayConfig {
    fn default() -> Self {
        Self {
            point_size: default_point_size(),
            blob_size: default_blob_size(),
            blob_alpha: default_blob_alpha(),
            blob_falloff: default_blob_falloff(),
            camera_distance: default_camera_distance(),
            camera_angle: default_camera_angle(),
            colormap_range: default_colormap_range(),
            jitter: default_jitter(),
            gif_resolution: default_gif_resolution(),
            gif_point_radius: default_gif_point_radius(),
            species: BTreeMap::new(),
        }
    }
}

impl DisplayConfig {
    /// Camera eye position as [x, y, z].
    pub fn camera_eye(&self) -> [f32; 3] {
        let d = self.camera_distance;
        [
            d * self.camera_angle[0],
            d * self.camera_angle[1],
            d * self.camera_angle[2],
        ]
    }
}
