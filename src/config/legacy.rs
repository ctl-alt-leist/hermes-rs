/// Legacy configuration types.
///
/// These types parse the old flat TOML format (configs/defaults.toml).
/// They remain in use by the existing simulation driver, scenes, and
/// runner code. New code should use EngineConfig instead.
use std::path::Path;

use serde::Deserialize;

use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;

// ============================================================================
// Configuration types
// ============================================================================

/// Top-level configuration for a hermes simulation run (legacy format).
#[derive(Debug, Clone, Deserialize)]
pub struct Configuration {
    pub cosmology: Cosmology,
    pub simulation: SimulationConfig,
    pub time: TimeConfig,
    pub output: OutputConfig,
    #[serde(default)]
    pub initialization: InitializationConfig,
    #[serde(default)]
    pub field: FieldConfig,
    #[serde(default)]
    pub visualization: VisualizationConfig,
}

/// Spatial discretization parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct SimulationConfig {
    /// Number of grid cells per side (total cells = n_grid^3).
    pub n_grid: usize,
    /// Number of particles per side (total particles = n_particles^3).
    pub n_particles: usize,
    /// Comoving box side length in kpc.
    pub box_length: f64,
}

/// Time-stepping parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct TimeConfig {
    /// Scale factor range [initial, final] where a = 1/(1+z).
    pub scale_factor_range: [f64; 2],
    /// Number of time steps.
    pub n_steps: usize,
    /// Scale factor stepping: "log" (logarithmic) or "linear".
    pub scale_factor_stepping: String,
}

impl TimeConfig {
    /// Initial scale factor.
    pub fn scale_factor_initial(&self) -> f64 {
        self.scale_factor_range[0]
    }

    /// Final scale factor.
    pub fn scale_factor_final(&self) -> f64 {
        self.scale_factor_range[1]
    }
}

/// Output parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct OutputConfig {
    /// Save a snapshot every this many steps.
    pub write_interval: usize,
    /// Compute full diagnostics every this many steps.
    pub diagnostic_interval: usize,
}

/// Initial condition parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct InitializationConfig {
    /// Spectrum source: "power" (CDM P(k)) or "random" (synthetic red spectrum).
    pub spectrum: String,
    /// RMS amplitude of initial density perturbations.
    pub perturbation_amplitude: f64,
    /// Band-pass filter as [k_min / k_fundamental, k_max / k_nyquist].
    pub band_pass: [f64; 2],
}

impl Default for InitializationConfig {
    fn default() -> Self {
        Self {
            spectrum: "power".to_string(),
            perturbation_amplitude: 0.1,
            band_pass: [1.5, 0.5],
        }
    }
}

/// Field theory (Schrodinger-Poisson) parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct FieldConfig {
    /// Smoothing length ratio l/m (kpc^2 / Gyr).
    pub length_scale: f64,
    /// Field mass parameter (M_sun).
    pub mass: f64,
}

impl Default for FieldConfig {
    fn default() -> Self {
        Self {
            length_scale: 2000.0,
            mass: 1e10,
        }
    }
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
    /// Gaussian falloff rate for volumetric blobs.
    pub blob_falloff: f32,
    /// Camera distance from origin (box spans [-0.5, 0.5]).
    pub camera_distance: f32,
    /// Camera direction as [x, y, z] (multiplied by distance).
    pub camera_angle: [f32; 3],
    /// Colormap range as [floor, ceiling] in units of density / rho_mean.
    pub colormap_range: [f64; 2],
    /// Grid-point jitter as fraction of cell size.
    pub jitter: f64,
    /// Pixel resolution for GIF recording.
    pub gif_resolution: u32,
    /// Point radius in pixels for GIF rendering.
    pub gif_point_radius: i32,
}

impl Default for VisualizationConfig {
    fn default() -> Self {
        Self {
            point_size: 5.0,
            blob_size: 18.0,
            blob_alpha: 0.12,
            blob_falloff: 10.0,
            camera_distance: 1.9,
            camera_angle: [0.56, 0.42, 0.69],
            colormap_range: [0.3, 3.0],
            jitter: 0.3,
            gif_resolution: 512,
            gif_point_radius: 1,
        }
    }
}

// ============================================================================
// Conversion from EngineConfig
// ============================================================================

impl Configuration {
    /// Build a legacy Configuration from an EngineConfig.
    ///
    /// This bridges the new config format to the existing simulation
    /// driver, which still reads from the legacy fields.
    pub fn from_engine_config(
        engine: &super::EngineConfig,
        cosmology: &crate::physics::cosmology::Cosmology,
    ) -> Self {
        let grid = &engine.simulation.grid;
        let time = &engine.simulation.time;
        let init = &engine.simulation.initialization;

        // Determine n_particles from the particle species, if any.
        let n_particles = engine
            .ontology
            .particles
            .values()
            .next()
            .map(|p| (p.count as f64).cbrt().round() as usize)
            .unwrap_or(1);

        // Determine field params from the first field species, if any.
        let field_config = engine
            .ontology
            .fields
            .values()
            .next()
            .map(|f| FieldConfig {
                length_scale: f.length_scale.unwrap_or(2000.0),
                mass: f.mass.unwrap_or(1e10),
            })
            .unwrap_or_default();

        let scale_factor_range = time.scale_factor_range.unwrap_or([0.02, 1.0]);

        Self {
            cosmology: cosmology.clone(),
            simulation: SimulationConfig {
                n_grid: grid.n_cells,
                n_particles,
                box_length: grid.box_length,
            },
            time: TimeConfig {
                scale_factor_range,
                n_steps: time.n_steps,
                scale_factor_stepping: time.stepping.clone(),
            },
            output: OutputConfig {
                write_interval: engine.output.snapshots.interval,
                diagnostic_interval: engine.output.diagnostics.interval,
            },
            initialization: InitializationConfig {
                spectrum: init.spectrum.clone(),
                perturbation_amplitude: init.perturbation_amplitude,
                band_pass: init.band_pass,
            },
            field: field_config,
            visualization: VisualizationConfig {
                point_size: engine.output.display.point_size,
                blob_size: engine.output.display.blob_size,
                blob_alpha: engine.output.display.blob_alpha,
                blob_falloff: engine.output.display.blob_falloff,
                camera_distance: engine.output.display.camera_distance,
                camera_angle: engine.output.display.camera_angle,
                colormap_range: engine.output.display.colormap_range,
                jitter: engine.output.display.jitter,
                gif_resolution: engine.output.display.gif_resolution,
                gif_point_radius: engine.output.display.gif_point_radius,
            },
        }
    }
}

// ============================================================================
// Convenience accessors
// ============================================================================

impl SimulationConfig {
    /// Grid cells per side. Alias for n_grid.
    pub fn n_cells(&self) -> usize {
        self.n_grid
    }
}

impl VisualizationConfig {
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

// ============================================================================
// Embedded defaults
// ============================================================================

const DEFAULTS_TOML: &str = include_str!("defaults.toml");

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

/// Recursively merge `source` into `base`, overwriting leaf values.
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
