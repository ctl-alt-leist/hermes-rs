//! Simulation snapshots for storage and visualization.
//!
//! A `Snapshot` captures the complete state of the simulation at one moment:
//! particle positions and momenta as morphis grade-1 vectors, diagnostics,
//! and metadata. Snapshots are the data contract between the simulation,
//! the file system, and the viewer.
//!
//! For disk serialization, morphis vectors are converted to flat arrays at
//! the I/O boundary (morphis types do not implement Serialize). The in-memory
//! representation is morphis-native.

use morphis::vector::Vector;
use serde::{Deserialize, Serialize};

use crate::algebra::{components_from_vector, vector_from_array};
use crate::physics::diagnostics::Diagnostics;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// In-memory snapshot of simulation state at one point in time.
///
/// All physical quantities are morphis objects: positions and momenta
/// are `Vector<3>` (grade-1), angular momentum is `Vector<3>` (grade-2).
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// Time step index.
    pub step: usize,
    /// Scale factor a.
    pub scale_factor: f64,
    /// Particle positions as morphis grade-1 vectors.
    pub positions: Vec<Vector<3>>,
    /// Particle momenta as morphis grade-1 vectors.
    pub momenta: Vec<Vector<3>>,
    /// Mass per particle (M_☉).
    pub mass_particle: f64,
    /// Conservation diagnostics.
    pub diagnostics: DiagnosticsSummary,
}

/// Serializable summary of diagnostics (no morphis types).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiagnosticsSummary {
    pub scale_factor: f64,
    pub mass_total: f64,
    pub momentum_magnitude: f64,
    pub energy_kinetic: f64,
    pub energy_potential: f64,
    pub angular_momentum_magnitude: f64,
}

impl DiagnosticsSummary {
    /// Extract a serializable summary from full diagnostics.
    pub fn from_diagnostics(diagnostics: &Diagnostics) -> Self {
        Self {
            scale_factor: diagnostics.scale_factor,
            mass_total: diagnostics.mass_total,
            momentum_magnitude: diagnostics.momentum_magnitude(),
            energy_kinetic: diagnostics.energy_kinetic,
            energy_potential: diagnostics.energy_potential,
            angular_momentum_magnitude: diagnostics.angular_momentum_magnitude(),
        }
    }
}

/// Serializable form of a snapshot for disk storage.
///
/// Morphis vectors are flattened to `[f64; 3]` arrays at this boundary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnapshotOnDisk {
    pub step: usize,
    pub scale_factor: f64,
    pub positions: Vec<[f64; 3]>,
    pub momenta: Vec<[f64; 3]>,
    pub mass_particle: f64,
    pub diagnostics: DiagnosticsSummary,
}

impl Snapshot {
    /// Capture a full snapshot including diagnostics (Poisson solve).
    ///
    /// This is expensive (~250ms for 32³) because it computes the
    /// gravitational potential for the energy diagnostic. Use
    /// `capture_lightweight` for the viewer path.
    pub fn capture(
        particles: &Particles,
        grid: &Grid,
        cosmology: &crate::physics::cosmology::Cosmology,
        solver: &mut crate::physics::poisson::PoissonSolver,
        step: usize,
        scale_factor: f64,
    ) -> Self {
        let diag = Diagnostics::compute(particles, grid, cosmology, solver, scale_factor);

        Self::capture_with_diagnostics(particles, step, scale_factor, &diag)
    }

    /// Capture a lightweight snapshot — just positions and momenta.
    ///
    /// Skips the Poisson solve entirely. Diagnostics are populated with
    /// only the quantities that are cheap to compute (mass, momentum,
    /// kinetic energy). Suitable for the live viewer path where every
    /// frame matters.
    pub fn capture_lightweight(particles: &Particles, step: usize, scale_factor: f64) -> Self {
        let positions = (0..particles.count())
            .map(|p| particles.position_of(p))
            .collect();

        let momenta = (0..particles.count())
            .map(|p| particles.momentum_of(p))
            .collect();

        let energy_kinetic = particles.kinetic_energy(scale_factor);

        Self {
            step,
            scale_factor,
            positions,
            momenta,
            mass_particle: particles.mass_particle,
            diagnostics: DiagnosticsSummary {
                scale_factor,
                mass_total: particles.total_mass(),
                momentum_magnitude: particles.total_momentum().norm(),
                energy_kinetic,
                energy_potential: 0.0,
                angular_momentum_magnitude: 0.0,
            },
        }
    }

    /// Build a snapshot from pre-computed diagnostics.
    fn capture_with_diagnostics(
        particles: &Particles,
        step: usize,
        scale_factor: f64,
        diagnostics: &Diagnostics,
    ) -> Self {
        let positions = (0..particles.count())
            .map(|p| particles.position_of(p))
            .collect();

        let momenta = (0..particles.count())
            .map(|p| particles.momentum_of(p))
            .collect();

        Self {
            step,
            scale_factor,
            positions,
            momenta,
            mass_particle: particles.mass_particle,
            diagnostics: DiagnosticsSummary::from_diagnostics(diagnostics),
        }
    }

    /// Convert to the serializable disk format.
    pub fn to_disk(&self) -> SnapshotOnDisk {
        SnapshotOnDisk {
            step: self.step,
            scale_factor: self.scale_factor,
            positions: self.positions.iter().map(components_from_vector).collect(),
            momenta: self.momenta.iter().map(components_from_vector).collect(),
            mass_particle: self.mass_particle,
            diagnostics: self.diagnostics.clone(),
        }
    }

    /// Reconstruct from the disk format, restoring morphis vectors.
    pub fn from_disk(disk: SnapshotOnDisk) -> Self {
        Self {
            step: disk.step,
            scale_factor: disk.scale_factor,
            positions: disk
                .positions
                .iter()
                .map(|c| vector_from_array(*c))
                .collect(),
            momenta: disk.momenta.iter().map(|c| vector_from_array(*c)).collect(),
            mass_particle: disk.mass_particle,
            diagnostics: disk.diagnostics,
        }
    }

    /// Number of particles in the snapshot.
    pub fn particle_count(&self) -> usize {
        self.positions.len()
    }
}

/// Save a snapshot to disk as bincode.
///
/// Creates parent directories if they don't exist.
pub fn save_snapshot(
    snapshot: &Snapshot,
    path: &std::path::Path,
) -> Result<(), crate::error::HermesError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let disk = snapshot.to_disk();
    let file = std::fs::File::create(path)?;
    bincode::serialize_into(file, &disk)
        .map_err(|e| crate::error::HermesError::Config(format!("bincode serialize failed: {e}")))?;

    Ok(())
}

/// Load a snapshot from a bincode file on disk.
pub fn load_snapshot(path: &std::path::Path) -> Result<Snapshot, crate::error::HermesError> {
    let file = std::fs::File::open(path)?;
    let disk: SnapshotOnDisk = bincode::deserialize_from(file).map_err(|e| {
        crate::error::HermesError::Config(format!("bincode deserialize failed: {e}"))
    })?;

    Ok(Snapshot::from_disk(disk))
}
