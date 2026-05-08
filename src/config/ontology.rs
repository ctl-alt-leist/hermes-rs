/// Ontology configuration: what exists and what governs it.
///
/// Parsed from the `[ontology]` section of the TOML config. Defines the
/// background spacetime, the species in the box (particles and fields),
/// and the Lagrangian terms that govern their dynamics and interactions.
use std::collections::BTreeMap;

use serde::Deserialize;

use crate::error::HermesError;

// ============================================================================
// Top-level ontology
// ============================================================================

/// The full ontology of a simulation: spacetime, species, and Lagrangian.
#[derive(Debug, Clone, Deserialize)]
pub struct Ontology {
    pub spacetime: Spacetime,
    #[serde(default)]
    pub particles: BTreeMap<String, ParticleSpecies>,
    #[serde(default)]
    pub fields: BTreeMap<String, FieldSpecies>,
    #[serde(default)]
    pub lagrangian: Lagrangian,
}

// ============================================================================
// Spacetime
// ============================================================================

/// Background spacetime geometry.
#[derive(Debug, Clone, Deserialize)]
pub struct Spacetime {
    /// Background type: "flrw" or "static".
    pub background: SpacetimeBackground,
    /// Hubble constant in km/s/Mpc. Only used for FLRW.
    pub hubble: Option<f64>,
    /// Total matter density parameter. Only used for FLRW.
    pub omega_m: Option<f64>,
    /// Vacuum energy density parameter. Only used for FLRW.
    pub omega_v: Option<f64>,
    /// Baryon density parameter. Only used for FLRW.
    pub omega_b: Option<f64>,
    /// Radiation density parameter. Only used for FLRW.
    pub omega_r: Option<f64>,
    /// Curvature density parameter. Only used for FLRW.
    pub omega_k: Option<f64>,
    /// Amplitude of matter fluctuations at 8 Mpc/h. Only used for FLRW.
    pub sigma_8: Option<f64>,
    /// Primordial power spectrum tilt. Only used for FLRW.
    pub spectral_index: Option<f64>,
}

/// The background geometry type.
#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum SpacetimeBackground {
    Flrw,
    Static,
}

impl Spacetime {
    /// Validate spacetime parameters.
    pub fn validate(&self) -> Result<(), HermesError> {
        match self.background {
            SpacetimeBackground::Flrw => {
                let required = [
                    ("hubble", self.hubble),
                    ("omega_m", self.omega_m),
                    ("omega_v", self.omega_v),
                ];
                for (name, value) in &required {
                    if value.is_none() {
                        return Err(HermesError::Config(format!(
                            "FLRW spacetime requires {name}"
                        )));
                    }
                }
            }
            SpacetimeBackground::Static => {}
        }

        Ok(())
    }

    /// Whether the spacetime is expanding.
    pub fn is_expanding(&self) -> bool {
        self.background == SpacetimeBackground::Flrw
    }
}

// ============================================================================
// Particle species
// ============================================================================

/// A particle species declaration.
#[derive(Debug, Clone, Deserialize)]
pub struct ParticleSpecies {
    /// Total particle count.
    pub count: usize,
    /// Mass per particle in M_sun.
    pub mass: f64,
    /// Electric charge (for electromagnetic coupling). Defaults to 0.
    #[serde(default)]
    pub charge: f64,
}

// ============================================================================
// Field species
// ============================================================================

/// A field species declaration.
#[derive(Debug, Clone, Deserialize)]
pub struct FieldSpecies {
    /// Algebraic grade(s). Single integer for a pure grade (e.g. 0, 2),
    /// array for a multi-grade subspace (e.g. [0, 3] for even subalgebra).
    pub grade: FieldGrade,
    /// Field mass parameter in M_sun. Not required for all field types.
    pub mass: Option<f64>,
    /// Diffusivity l/m in kpc^2 / Gyr. Only for Schrodinger fields.
    pub length_scale: Option<f64>,
    /// Free Lagrangian dynamics: "schrodinger", "euler", "maxwell", "wave".
    pub free: Option<String>,
    /// Propagation speed in km/s. Only for wave fields.
    pub speed: Option<f64>,
    /// Gross-Pitaevskii self-interaction coupling constant.
    /// Units: kpc^3 / Gyr^2 / M_sun. Only for Schrodinger fields.
    pub self_interaction: Option<f64>,
    /// Electric charge (for electromagnetic coupling). Defaults to 0.
    #[serde(default)]
    pub charge: f64,
}

/// Algebraic grade specification: single grade or multi-grade subspace.
#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum FieldGrade {
    /// A single grade (e.g. 0 for scalar, 2 for bivector).
    Single(usize),
    /// Multiple grades forming a subspace (e.g. [0, 3] for even subalgebra).
    Multi(Vec<usize>),
}

impl FieldGrade {
    /// Whether this is the even subalgebra (grades 0 and 3 in 3D).
    pub fn is_even_subalgebra(&self) -> bool {
        match self {
            FieldGrade::Multi(grades) => grades == &[0, 3],
            _ => false,
        }
    }

    /// Whether this is a single grade.
    pub fn is_single(&self) -> bool {
        matches!(self, FieldGrade::Single(_))
    }
}

// ============================================================================
// Lagrangian
// ============================================================================

/// Coupling terms between species.
///
/// Free terms live on each species entry (the `free` key);
/// cross-species interactions live here.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct Lagrangian {
    /// Gravitational coupling. `true` means universal (all massive species).
    #[serde(default)]
    pub gravity: bool,
    /// Electromagnetic coupling. Lists the species names that participate.
    #[serde(default)]
    pub electromagnetic: Vec<String>,
}

// ============================================================================
// Validation
// ============================================================================

impl Ontology {
    /// Validate the full ontology for internal consistency.
    pub fn validate(&self) -> Result<(), HermesError> {
        self.spacetime.validate()?;

        // Validate field species
        for (name, field) in &self.fields {
            if field.free.as_deref() == Some("schrodinger") && field.length_scale.is_none() {
                return Err(HermesError::Config(format!(
                    "field '{name}' has free = \"schrodinger\" but no length_scale"
                )));
            }
            if field.free.as_deref() == Some("wave") && field.speed.is_none() {
                return Err(HermesError::Config(format!(
                    "field '{name}' has free = \"wave\" but no speed"
                )));
            }
        }

        // Validate electromagnetic coupling references
        let all_species: Vec<&str> = self
            .particles
            .keys()
            .chain(self.fields.keys())
            .map(|s| s.as_str())
            .collect();

        for species_name in &self.lagrangian.electromagnetic {
            if !all_species.contains(&species_name.as_str()) {
                return Err(HermesError::Config(format!(
                    "electromagnetic coupling references unknown species '{species_name}'"
                )));
            }
        }

        // Must have at least one species
        if self.particles.is_empty() && self.fields.is_empty() {
            return Err(HermesError::Config(
                "ontology must declare at least one species".to_string(),
            ));
        }

        Ok(())
    }

    /// Whether the simulation has any particle species.
    pub fn has_particles(&self) -> bool {
        !self.particles.is_empty()
    }

    /// Whether the simulation has any field species.
    pub fn has_fields(&self) -> bool {
        !self.fields.is_empty()
    }

    /// Whether gravity is enabled.
    pub fn has_gravity(&self) -> bool {
        self.lagrangian.gravity
    }
}
