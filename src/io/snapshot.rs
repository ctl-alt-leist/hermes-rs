//! Simulation snapshots for storage and visualization.
//!
//! A `Snapshot` captures the simulation state at one moment. The content
//! can be particles (positions + momenta), fields (density grid), or both.
//! Snapshots are the data contract between the simulation, the file system,
//! and the viewer.

use morphis::vector::Vector;
use serde::{Deserialize, Serialize};

use crate::algebra::{components_from_vector, vector_from_array};
use crate::core::content::Content;
use crate::physics::diagnostics::Diagnostics;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

// ============================================================================
// Snapshot
// ============================================================================

/// In-memory snapshot of simulation state.
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// Time step index.
    pub step: usize,
    /// Scale factor a.
    pub scale_factor: f64,
    /// Content-specific data.
    pub content: SnapshotContent,
    /// Conservation diagnostics.
    pub diagnostics: DiagnosticsSummary,
}

/// Content-specific snapshot data.
#[derive(Debug, Clone)]
pub enum SnapshotContent {
    /// Particle positions and momenta as morphis grade-1 vectors.
    Particles {
        positions: Vec<Vector<3>>,
        momenta: Vec<Vector<3>>,
        mass_particle: f64,
    },
    /// Field density on a grid.
    Fields {
        /// Density values at each grid point (flattened).
        density: Vec<f64>,
        /// Grid cells per side.
        n_cells: usize,
    },
}

impl Snapshot {
    /// Capture a lightweight snapshot from Content.
    ///
    /// For particles: copies positions and momenta.
    /// For fields: extracts density.
    /// Skips expensive diagnostics (Poisson solve).
    pub fn capture_from_content(content: &Content, step: usize, scale_factor: f64) -> Self {
        match content {
            Content::Particles(particles) => {
                Self::capture_lightweight(particles, step, scale_factor)
            }
            Content::Fields(field_state) => {
                let (density, n_cells) = if let Some(ref psi) = field_state.psi {
                    let density_field = psi.density(field_state.params.mass_alpha);
                    let n = field_state.grid.n_cells;
                    let mut density = Vec::with_capacity(n * n * n);
                    for val in density_field.data.iter() {
                        density.push(*val);
                    }
                    (density, n)
                } else {
                    (Vec::new(), 0)
                };

                Self {
                    step,
                    scale_factor,
                    content: SnapshotContent::Fields { density, n_cells },
                    diagnostics: DiagnosticsSummary {
                        scale_factor,
                        mass_total: 0.0,
                        momentum_magnitude: 0.0,
                        energy_kinetic: 0.0,
                        energy_potential: 0.0,
                        angular_momentum_magnitude: 0.0,
                    },
                }
            }
            Content::Mixed { particles, .. } => {
                // For now, snapshot the particle part only.
                Self::capture_lightweight(particles, step, scale_factor)
            }
        }
    }

    /// Capture a lightweight particle snapshot.
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
            content: SnapshotContent::Particles {
                positions,
                momenta,
                mass_particle: particles.mass_particle,
            },
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

    /// Capture a full particle snapshot including diagnostics.
    pub fn capture(
        particles: &Particles,
        grid: &Grid,
        cosmology: &crate::physics::cosmology::Cosmology,
        solver: &mut crate::physics::poisson::PoissonSolver,
        step: usize,
        scale_factor: f64,
    ) -> Self {
        let diag = Diagnostics::compute(particles, grid, cosmology, solver, scale_factor);

        let positions = (0..particles.count())
            .map(|p| particles.position_of(p))
            .collect();

        let momenta = (0..particles.count())
            .map(|p| particles.momentum_of(p))
            .collect();

        Self {
            step,
            scale_factor,
            content: SnapshotContent::Particles {
                positions,
                momenta,
                mass_particle: particles.mass_particle,
            },
            diagnostics: DiagnosticsSummary::from_diagnostics(&diag),
        }
    }

    /// Number of particles (if particle content).
    pub fn particle_count(&self) -> usize {
        match &self.content {
            SnapshotContent::Particles { positions, .. } => positions.len(),
            SnapshotContent::Fields { .. } => 0,
        }
    }

    /// Access particle positions (if particle content).
    pub fn positions(&self) -> Option<&[Vector<3>]> {
        match &self.content {
            SnapshotContent::Particles { positions, .. } => Some(positions),
            _ => None,
        }
    }

    /// Access particle momenta (if particle content).
    pub fn momenta(&self) -> Option<&[Vector<3>]> {
        match &self.content {
            SnapshotContent::Particles { momenta, .. } => Some(momenta),
            _ => None,
        }
    }

    /// Reconstruct a Particles object from a particle snapshot.
    ///
    /// Returns None if the snapshot contains field content.
    pub fn to_particles(&self) -> Option<Particles> {
        match &self.content {
            SnapshotContent::Particles {
                positions,
                momenta,
                mass_particle,
            } => {
                let n = positions.len();
                let mut particles = Particles::zeros(n, *mass_particle);
                for (p, (pos, mom)) in positions.iter().zip(momenta.iter()).enumerate() {
                    particles.set_position(p, pos);
                    particles.set_momentum(p, mom);
                }

                Some(particles)
            }
            _ => None,
        }
    }

    /// Convert to serializable disk format.
    pub fn to_disk(&self) -> SnapshotOnDisk {
        let content = match &self.content {
            SnapshotContent::Particles {
                positions,
                momenta,
                mass_particle,
            } => SnapshotContentOnDisk::Particles {
                positions: positions.iter().map(components_from_vector).collect(),
                momenta: momenta.iter().map(components_from_vector).collect(),
                mass_particle: *mass_particle,
            },
            SnapshotContent::Fields { density, n_cells } => SnapshotContentOnDisk::Fields {
                density: density.clone(),
                n_cells: *n_cells,
            },
        };

        SnapshotOnDisk {
            step: self.step,
            scale_factor: self.scale_factor,
            content,
            diagnostics: self.diagnostics.clone(),
        }
    }

    /// Reconstruct from disk format.
    pub fn from_disk(disk: SnapshotOnDisk) -> Self {
        let content = match disk.content {
            SnapshotContentOnDisk::Particles {
                positions,
                momenta,
                mass_particle,
            } => SnapshotContent::Particles {
                positions: positions.iter().map(|c| vector_from_array(*c)).collect(),
                momenta: momenta.iter().map(|c| vector_from_array(*c)).collect(),
                mass_particle,
            },
            SnapshotContentOnDisk::Fields { density, n_cells } => {
                SnapshotContent::Fields { density, n_cells }
            }
        };

        Self {
            step: disk.step,
            scale_factor: disk.scale_factor,
            content,
            diagnostics: disk.diagnostics,
        }
    }
}

// ============================================================================
// Diagnostics summary
// ============================================================================

/// Serializable summary of diagnostics.
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
    /// Extract from full diagnostics.
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

// ============================================================================
// On-disk format
// ============================================================================

/// Serializable snapshot for disk storage.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnapshotOnDisk {
    pub step: usize,
    pub scale_factor: f64,
    pub content: SnapshotContentOnDisk,
    pub diagnostics: DiagnosticsSummary,
}

/// Serializable content (morphis vectors flattened to [f64; 3]).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SnapshotContentOnDisk {
    Particles {
        positions: Vec<[f64; 3]>,
        momenta: Vec<[f64; 3]>,
        mass_particle: f64,
    },
    Fields {
        density: Vec<f64>,
        n_cells: usize,
    },
}

// ============================================================================
// I/O
// ============================================================================

/// Save a snapshot to disk as bincode.
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
