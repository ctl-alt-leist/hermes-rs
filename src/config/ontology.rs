/// Ontology configuration: what exists and what governs it.
///
/// Parsed from the `[ontology]` section of the TOML config. Defines the
/// background spacetime, the species in the box (particles and fields),
/// and the coupling terms that govern their interactions.
use std::collections::BTreeMap;

use serde::Deserialize;

use crate::error::HermesError;

// ============================================================================
// Top-level ontology
// ============================================================================

/// The full ontology of a simulation: spacetime, species, and couplings.
#[derive(Debug, Clone, Deserialize)]
pub struct Ontology {
    pub spacetime: Spacetime,
    #[serde(default)]
    pub particles: BTreeMap<String, ParticleSpecies>,
    #[serde(default)]
    pub fields: BTreeMap<String, FieldSpecies>,
    /// Coupling terms between species.
    #[serde(default)]
    pub coupling: Vec<Coupling>,
    /// Legacy lagrangian block (supported for backward compatibility).
    #[serde(default)]
    pub lagrangian: Option<LagrangianLegacy>,
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
    /// Display symbol for diagnostics and visualization labels.
    pub symbol: Option<String>,
    /// Particles per side (total count = n^3).
    pub n: usize,
    /// Mass per particle in M_sun.
    pub mass: f64,
    /// Gravitational softening length in kpc.
    pub softening: Option<f64>,
    /// Deposition kernel: "cic", "tsc", "pcs". Default: "cic".
    #[serde(default = "default_kernel")]
    pub kernel: String,
}

fn default_kernel() -> String {
    "cic".to_string()
}

impl ParticleSpecies {
    /// Total particle count (n^3).
    pub fn total_count(&self) -> usize {
        self.n * self.n * self.n
    }
}

// ============================================================================
// Field species
// ============================================================================

/// A field species declaration.
#[derive(Debug, Clone, Deserialize)]
pub struct FieldSpecies {
    /// Display symbol for diagnostics and visualization labels.
    pub symbol: Option<String>,
    /// Algebraic grade(s). Single integer for a pure grade (e.g. 0, 2),
    /// array for a multi-grade subspace (e.g. [0, 3] for even subalgebra).
    pub grade: FieldGrade,
    /// Field mass parameter in M_sun. Not required for all field types.
    pub mass: Option<f64>,
    /// Diffusivity l/m in kpc^2 / Gyr. Only for Schrodinger fields.
    pub length_scale: Option<f64>,
    /// Free Lagrangian dynamics: "schrodinger", "wave".
    pub free: Option<String>,
    /// Propagation speed in km/s. Only for wave fields.
    pub speed: Option<f64>,
    /// Gross-Pitaevskii self-interaction coupling constant.
    /// Units: kpc^3 / Gyr^2 / M_sun. Only for Schrodinger fields.
    pub self_interaction: Option<f64>,
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
// Couplings
// ============================================================================

/// A coupling term between species.
#[derive(Debug, Clone, Deserialize)]
pub struct Coupling {
    /// Coupling kind: "gravity", "electromagnetic", etc.
    pub kind: String,
    /// Participating species (by name).
    #[serde(default)]
    pub species: Vec<String>,
}

/// Legacy lagrangian block (backward compatibility).
#[derive(Debug, Clone, Default, Deserialize)]
pub struct LagrangianLegacy {
    #[serde(default)]
    pub gravity: bool,
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

        // Validate field species.
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

        // Validate coupling references.
        let all_species: Vec<&str> = self
            .particles
            .keys()
            .chain(self.fields.keys())
            .map(|s| s.as_str())
            .collect();

        for coupling in &self.coupling {
            for species_name in &coupling.species {
                if !all_species.contains(&species_name.as_str()) {
                    return Err(HermesError::Config(format!(
                        "{} coupling references unknown species '{species_name}'",
                        coupling.kind
                    )));
                }
            }
        }

        // Must have at least one species.
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

    /// Whether gravity is enabled (via coupling list or legacy lagrangian).
    pub fn has_gravity(&self) -> bool {
        self.coupling.iter().any(|c| c.kind == "gravity")
            || self.lagrangian.as_ref().is_some_and(|l| l.gravity)
    }
}
