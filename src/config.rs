/// Configuration loading and management.
///
/// Four-tier hierarchy: embedded defaults → scene defaults → optional config
/// file → CLI overrides. Partial TOML files are deep-merged so that only the
/// fields being overridden need to appear.
use std::path::Path;

use serde::Deserialize;

use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;

// ============================================================================
// Configuration types
// ============================================================================

/// Top-level configuration for a hermes simulation run.
#[derive(Debug, Clone, Deserialize)]
pub struct Configuration {
    pub cosmology: Cosmology,
    pub simulation: SimulationConfig,
    pub time: TimeConfig,
    pub output: OutputConfig,
    #[serde(default)]
    pub visualization: VisualizationConfig,
}

/// Spatial discretization parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct SimulationConfig {
    /// Number of grid cells per side (total cells = n_grid³).
    pub n_grid: usize,
    /// Number of particles per side (total particles = n_particles³).
    pub n_particles: usize,
    /// Comoving box side length in kpc.
    pub box_length: f64,
}

/// Time-stepping parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct TimeConfig {
    /// Initial scale factor (a = 1/(1+z)).
    pub scale_factor_initial: f64,
    /// Final scale factor.
    pub scale_factor_final: f64,
    /// Number of time steps.
    pub n_steps: usize,
    /// Stepping strategy: "log_a" (logarithmic in a) or "linear_a".
    pub stepping: String,
}

/// Output parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct OutputConfig {
    /// Save a snapshot every this many steps.
    pub write_interval: usize,
    /// Compute full diagnostics every this many steps.
    pub diagnostic_interval: usize,
}

/// Visualization parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct VisualizationConfig {
    /// Screen-space point size for particle rendering.
    pub point_size: f32,
    /// Screen-space blob size for volumetric field rendering.
    pub blob_size: f32,
    /// Per-blob opacity for additive blending.
    pub blob_alpha: f32,
    /// Camera distance from origin (box spans [-0.5, 0.5]).
    pub camera_distance: f32,
    /// Density / rho_mean for colormap floor.
    pub colormap_low: f64,
    /// Density / rho_mean for colormap ceiling.
    pub colormap_high: f64,
    /// Grid-point jitter as fraction of cell size.
    pub jitter: f64,
}

impl Default for VisualizationConfig {
    fn default() -> Self {
        Self {
            point_size: 5.0,
            blob_size: 18.0,
            blob_alpha: 0.12,
            camera_distance: 1.9,
            colormap_low: 0.3,
            colormap_high: 3.0,
            jitter: 0.3,
        }
    }
}

// ============================================================================
// Backwards compatibility
// ============================================================================

impl SimulationConfig {
    /// Grid cells per side. Alias for n_grid.
    pub fn n_cells(&self) -> usize {
        self.n_grid
    }
}

// ============================================================================
// Embedded defaults
// ============================================================================

const DEFAULTS_TOML: &str = include_str!("../configs/defaults.toml");

// ============================================================================
// Loading
// ============================================================================

/// Load the embedded default configuration.
pub fn load_defaults() -> Result<Configuration, HermesError> {
    let config: Configuration = toml::from_str(DEFAULTS_TOML)?;
    config.cosmology.validate()?;

    Ok(config)
}

/// Load configuration from a TOML file, merged on top of defaults.
pub fn load_config(path: &Path) -> Result<Configuration, HermesError> {
    let override_str = std::fs::read_to_string(path)?;
    let override_val: toml::Value = toml::from_str(&override_str)?;

    build_configuration(None, Some(&override_val))
}

/// Build a configuration from defaults with optional overrides.
///
/// The `config_file` override is applied first, then `overrides`. Both
/// are deep-merged into the defaults so that only the fields being
/// changed need to appear.
pub fn build_configuration(
    config_file: Option<&toml::Value>,
    overrides: Option<&toml::Value>,
) -> Result<Configuration, HermesError> {
    let mut base: toml::Value = toml::from_str(DEFAULTS_TOML)
        .map_err(|e| HermesError::Config(format!("failed to parse embedded defaults: {e}")))?;

    if let Some(file_val) = config_file {
        deep_merge(&mut base, file_val);
    }
    if let Some(override_val) = overrides {
        deep_merge(&mut base, override_val);
    }

    let config: Configuration = base
        .try_into()
        .map_err(|e| HermesError::Config(format!("failed to deserialize merged config: {e}")))?;

    config.cosmology.validate()?;

    Ok(config)
}

// ============================================================================
// Deep merge
// ============================================================================

/// Recursively merge `source` into `base`, overwriting leaf values.
///
/// Public for use by the runner's four-tier config merge.
pub fn deep_merge_public(base: &mut toml::Value, source: &toml::Value) {
    deep_merge(base, source);
}

/// Recursively merge `source` into `base`, overwriting leaf values.
fn deep_merge(base: &mut toml::Value, source: &toml::Value) {
    if let (toml::Value::Table(base_table), toml::Value::Table(source_table)) = (base, source) {
        for (key, source_val) in source_table {
            if let Some(base_val) = base_table.get_mut(key) {
                if base_val.is_table() && source_val.is_table() {
                    deep_merge(base_val, source_val);
                } else {
                    *base_val = source_val.clone();
                }
            } else {
                base_table.insert(key.clone(), source_val.clone());
            }
        }
    }
}
