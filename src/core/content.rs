//! Content abstraction: what a simulation carries and evolves.
//!
//! A simulation's state is its content — particles, fields, or both.
//! The content enum enables a single simulation driver to handle all
//! three content kinds without knowing the specifics.

use morphis::even_field::EvenField;
use morphis::field::Field;

use crate::physics::particles::Particles;

/// What a simulation carries.
pub enum Content {
    /// Dark matter particles with PM gravity.
    Particles(Particles),
    /// Grade-aware fields (wavefunctions, bivector fields).
    Fields(FieldState),
    /// Both — coupled through shared gravitational potential.
    Mixed {
        particles: Particles,
        fields: FieldState,
    },
}

impl Content {
    /// Extract particles if present.
    pub fn particles(&self) -> Option<&Particles> {
        match self {
            Content::Particles(p) => Some(p),
            Content::Mixed { particles, .. } => Some(particles),
            Content::Fields(_) => None,
        }
    }

    /// Extract particles mutably if present.
    pub fn particles_mut(&mut self) -> Option<&mut Particles> {
        match self {
            Content::Particles(p) => Some(p),
            Content::Mixed { particles, .. } => Some(particles),
            Content::Fields(_) => None,
        }
    }

    /// Extract field state if present.
    pub fn fields(&self) -> Option<&FieldState> {
        match self {
            Content::Fields(f) => Some(f),
            Content::Mixed { fields, .. } => Some(fields),
            Content::Particles(_) => None,
        }
    }

    /// Extract field state mutably if present.
    pub fn fields_mut(&mut self) -> Option<&mut FieldState> {
        match self {
            Content::Fields(f) => Some(f),
            Content::Mixed { fields, .. } => Some(fields),
            Content::Particles(_) => None,
        }
    }
}

/// Field-theoretic state: wavefunctions and bivector fields.
pub struct FieldState {
    /// Morphis grid for spectral operations.
    pub grid: morphis::grid::Grid<3>,
    /// Dark matter wavefunction (even subalgebra: scalar + pseudoscalar).
    pub psi: Option<EvenField<3>>,
    /// Baryon wavefunction (future).
    pub beta: Option<EvenField<3>>,
    /// Electromagnetic bivector field (future).
    pub gamma: Option<Field<3>>,
    /// Field physics parameters.
    pub params: FieldParams,
}

/// Parameters for field-theoretic dynamics.
pub struct FieldParams {
    /// Effective Planck constant (sets de Broglie wavelength scale).
    pub hbar_eff: f64,
    /// Dark matter field particle mass.
    pub mass_alpha: f64,
}

impl Default for FieldParams {
    fn default() -> Self {
        Self {
            hbar_eff: 1.0,
            mass_alpha: 1.0,
        }
    }
}
