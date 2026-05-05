/// Simulation state: what is in the box at a moment in time.
///
/// A named collection of particle species and field species on a shared
/// grid. The engine evolves this state forward by composing free evolution
/// and coupling steps.
use std::collections::BTreeMap;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;

use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// The simulation state: species on a grid at a moment in time.
pub struct SimulationState {
    /// Named particle species.
    pub particles: BTreeMap<String, Particles>,
    /// Named field species.
    pub fields: BTreeMap<String, FieldEntry>,
    /// Hermes grid (CIC, Poisson solver geometry).
    pub grid: Grid,
    /// Morphis grid (spectral field operations).
    pub morphis_grid: MorphisGrid<3>,
    /// Current scale factor (for FLRW) or coordinate time (for static).
    pub time: f64,
    /// Current step index.
    pub step: usize,
}

/// A single field species in the simulation.
pub struct FieldEntry {
    /// The field data (even subalgebra for now).
    pub data: EvenField<3>,
    /// Smoothing length l (= (l/m) * m).
    pub smoothing_length: f64,
    /// Field mass parameter in M_sun.
    pub mass: f64,
}

impl SimulationState {
    /// Whether the state has any particle species.
    pub fn has_particles(&self) -> bool {
        !self.particles.is_empty()
    }

    /// Whether the state has any field species.
    pub fn has_fields(&self) -> bool {
        !self.fields.is_empty()
    }

    /// Total particle count across all species.
    pub fn total_particle_count(&self) -> usize {
        self.particles.values().map(|p| p.count()).sum()
    }
}
