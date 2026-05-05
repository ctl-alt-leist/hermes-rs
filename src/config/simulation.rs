/// Simulation configuration: how we compute it.
///
/// Parsed from the `[simulation]` section of the TOML config. Covers
/// spatial grid, time stepping, and initial conditions.
use serde::Deserialize;

use crate::error::HermesError;

// ============================================================================
// Simulation block
// ============================================================================

/// Numerical parameters: grid, time stepping, and initialization.
#[derive(Debug, Clone, Deserialize)]
pub struct SimulationBlock {
    pub grid: GridConfig,
    pub time: TimeConfig,
    #[serde(default)]
    pub initialization: InitializationConfig,
}

/// Spatial grid parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct GridConfig {
    /// Grid cells per side (total cells = n_cells^3).
    pub n_cells: usize,
    /// Comoving box side length in kpc.
    pub box_length: f64,
}

/// Time-stepping parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct TimeConfig {
    /// Scale factor range [initial, final]. Used for FLRW spacetimes.
    pub scale_factor_range: Option<[f64; 2]>,
    /// Coordinate time range [initial, final]. Used for static spacetimes.
    pub time_range: Option<[f64; 2]>,
    /// Stepping mode: "log" or "linear".
    #[serde(default = "default_stepping")]
    pub stepping: String,
    /// Number of time steps.
    pub n_steps: usize,
}

fn default_stepping() -> String {
    "log".to_string()
}

impl TimeConfig {
    /// Validate time config against spacetime background.
    pub fn validate(&self, is_expanding: bool) -> Result<(), HermesError> {
        if is_expanding {
            if self.scale_factor_range.is_none() {
                return Err(HermesError::Config(
                    "FLRW spacetime requires scale_factor_range in [simulation.time]".to_string(),
                ));
            }
        } else if self.time_range.is_none() {
            return Err(HermesError::Config(
                "static spacetime requires time_range in [simulation.time]".to_string(),
            ));
        }

        Ok(())
    }
}

/// Initial condition parameters.
#[derive(Debug, Clone, Deserialize)]
pub struct InitializationConfig {
    /// Initialization method: "zeldovich", "nfw-group", "gaussian-packet".
    #[serde(default = "default_method")]
    pub method: String,
    /// Random seed.
    #[serde(default = "default_seed")]
    pub seed: u64,
    /// Power spectrum source: "power", "eisenstein-hu", "random".
    #[serde(default = "default_spectrum")]
    pub spectrum: String,
    /// RMS amplitude of initial density perturbations.
    #[serde(default = "default_perturbation_amplitude")]
    pub perturbation_amplitude: f64,
    /// Band-pass filter as [k_min / k_fundamental, k_max / k_nyquist].
    #[serde(default = "default_band_pass")]
    pub band_pass: [f64; 2],
    /// NFW halo specifications for the nfw-group method.
    #[serde(default)]
    pub halos: Vec<HaloSpec>,
}

/// NFW halo specification for the nfw-group initialization method.
#[derive(Debug, Clone, Deserialize)]
pub struct HaloSpec {
    /// Fraction of total box mass in this halo.
    pub mass_fraction: f64,
    /// NFW concentration parameter.
    pub concentration: f64,
}

fn default_method() -> String {
    "zeldovich".to_string()
}
fn default_seed() -> u64 {
    42
}
fn default_spectrum() -> String {
    "power".to_string()
}
fn default_perturbation_amplitude() -> f64 {
    0.1
}
fn default_band_pass() -> [f64; 2] {
    [1.5, 0.5]
}

impl Default for InitializationConfig {
    fn default() -> Self {
        Self {
            method: default_method(),
            seed: default_seed(),
            spectrum: default_spectrum(),
            perturbation_amplitude: default_perturbation_amplitude(),
            band_pass: default_band_pass(),
            halos: Vec::new(),
        }
    }
}
