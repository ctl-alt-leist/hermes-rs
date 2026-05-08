/// Configuration loading and management.
///
/// Three top-level containers mirror the TOML structure:
///   [ontology]    — what exists and what governs it
///   [simulation]  — how we compute it
///   [output]      — what we record, report, and display
///
/// The legacy `Configuration` type and its loading functions remain
/// available for existing code during the migration period.
pub mod ontology;
pub mod output;
pub mod simulation;

mod legacy;

// Re-export new config types at the config:: level.
pub use ontology::{
    Coupling, FieldGrade, FieldSpecies, LagrangianLegacy, Ontology, ParticleSpecies, Spacetime,
    SpacetimeBackground,
};
pub use output::{DiagnosticsConfig, DisplayConfig, LoggingConfig, OutputBlock, SnapshotsConfig};
pub use simulation::{GridConfig, HaloSpec, InitializationConfig, SimulationBlock, TimeConfig};

// Re-export legacy types so existing `use crate::config::*` still works.
pub use legacy::{
    Configuration, FieldConfig, OutputConfig, SimulationConfig, VisualizationConfig,
    build_configuration, deep_merge_public, load_config, load_defaults,
};

use std::path::Path;

use serde::Deserialize;

use crate::error::HermesError;

// ============================================================================
// Top-level engine config
// ============================================================================

/// Complete engine configuration parsed from the new TOML format.
#[derive(Debug, Clone, Deserialize)]
pub struct EngineConfig {
    pub ontology: Ontology,
    pub simulation: SimulationBlock,
    #[serde(default)]
    pub output: OutputBlock,
}

impl EngineConfig {
    /// Validate the full configuration for internal consistency.
    pub fn validate(&self) -> Result<(), HermesError> {
        self.ontology.validate()?;
        self.simulation
            .time
            .validate(self.ontology.spacetime.is_expanding())?;

        Ok(())
    }
}

// ============================================================================
// Loading
// ============================================================================

const BASE_DEFAULTS_TOML: &str = include_str!("base.toml");

/// Load the base defaults as an EngineConfig.
pub fn load_base_defaults() -> Result<EngineConfig, HermesError> {
    let config: EngineConfig = toml::from_str(BASE_DEFAULTS_TOML)?;
    config.ontology.validate()?;

    Ok(config)
}

/// Load a scene config file, merged on top of base defaults.
pub fn load_scene_config(path: &Path) -> Result<EngineConfig, HermesError> {
    let scene_str = std::fs::read_to_string(path)?;
    let scene_val: toml::Value = toml::from_str(&scene_str)?;

    let mut base: toml::Value = toml::from_str(BASE_DEFAULTS_TOML)
        .map_err(|e| HermesError::Config(format!("failed to parse base defaults: {e}")))?;

    deep_merge(&mut base, &scene_val);

    let config: EngineConfig = base
        .try_into()
        .map_err(|e| HermesError::Config(format!("failed to deserialize merged config: {e}")))?;

    config.ontology.validate()?;
    config
        .simulation
        .time
        .validate(config.ontology.spacetime.is_expanding())?;

    Ok(config)
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
