//! Simulation snapshots for storage and visualization.
//!
//! A `Snapshot` captures the simulation state at one moment. It holds
//! named species of both particles and fields — any combination,
//! including empty. Snapshots are the data contract between the
//! simulation, the file system, and the viewer.

use morphis::vector::Vector;
use serde::{Deserialize, Serialize};

use crate::algebra::{components_from_vector, vector_from_array};
use crate::core::content::Content;
use crate::engine::state::SimulationState;
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
    /// Named particle species.
    pub particles: Vec<ParticleSpeciesSnapshot>,
    /// Named field species.
    pub fields: Vec<FieldSpeciesSnapshot>,
    /// Grid cells per side (for field rendering). Zero if no fields.
    pub n_cells: usize,
    /// Conservation diagnostics.
    pub diagnostics: DiagnosticsSummary,
}

/// Snapshot of a single particle species.
#[derive(Debug, Clone)]
pub struct ParticleSpeciesSnapshot {
    /// Species name.
    pub name: String,
    /// Particle positions as morphis grade-1 vectors.
    pub positions: Vec<Vector<3>>,
    /// Particle momenta as morphis grade-1 vectors.
    pub momenta: Vec<Vector<3>>,
    /// Mass per particle in M_sun.
    pub mass_particle: f64,
}

/// Density snapshot of a single field species.
#[derive(Debug, Clone)]
pub struct FieldSpeciesSnapshot {
    /// Species name.
    pub name: String,
    /// Density values at each grid point (flattened, length n_cells^3).
    pub density: Vec<f64>,
}

impl Snapshot {
    /// Whether this snapshot contains any particle species.
    pub fn has_particles(&self) -> bool {
        !self.particles.is_empty()
    }

    /// Whether this snapshot contains any field species.
    pub fn has_fields(&self) -> bool {
        !self.fields.is_empty()
    }

    /// Total particle count across all species.
    pub fn particle_count(&self) -> usize {
        self.particles.iter().map(|s| s.positions.len()).sum()
    }

    // ========================================================================
    // Capture from the new engine's SimulationState
    // ========================================================================

    /// Capture all named species from the engine's SimulationState.
    pub fn capture_from_state(state: &SimulationState, step: usize, scale_factor: f64) -> Self {
        let mut particle_snapshots = Vec::new();
        for (name, particles) in &state.particles {
            let positions = (0..particles.count())
                .map(|p| particles.position_of(p))
                .collect();
            let momenta = (0..particles.count())
                .map(|p| particles.momentum_of(p))
                .collect();

            particle_snapshots.push(ParticleSpeciesSnapshot {
                name: name.clone(),
                positions,
                momenta,
                mass_particle: particles.mass_particle,
            });
        }

        let mut field_snapshots = Vec::new();
        for (name, field) in &state.fields {
            let density_field = field.data.density(field.mass);
            let density: Vec<f64> = density_field.data.iter().copied().collect();
            field_snapshots.push(FieldSpeciesSnapshot {
                name: name.clone(),
                density,
            });
        }

        let n_cells = state
            .fields
            .values()
            .next()
            .map(|f| f.data.grid.n_cells)
            .unwrap_or(0);

        Self {
            step,
            scale_factor,
            particles: particle_snapshots,
            fields: field_snapshots,
            n_cells,
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

    // ========================================================================
    // Capture from legacy Content (backward compatibility)
    // ========================================================================

    /// Capture a lightweight snapshot from legacy Content.
    pub fn capture_from_content(content: &Content, step: usize, scale_factor: f64) -> Self {
        match content {
            Content::Particles(particles) => {
                Self::capture_particles("dark_matter", particles, step, scale_factor)
            }
            Content::Fields(field_state) => {
                let n = field_state.grid.n_cells;
                let mut field_snapshots = Vec::new();

                if let Some(ref alpha) = field_state.alpha {
                    let density_field = alpha.density(field_state.params.mass_alpha);
                    let density: Vec<f64> = density_field.data.iter().copied().collect();
                    field_snapshots.push(FieldSpeciesSnapshot {
                        name: "dark matter".to_string(),
                        density,
                    });
                }

                if let Some(ref beta) = field_state.beta {
                    let density_field = beta.density(field_state.params.mass_alpha);
                    let density: Vec<f64> = density_field.data.iter().copied().collect();
                    field_snapshots.push(FieldSpeciesSnapshot {
                        name: "baryonic matter".to_string(),
                        density,
                    });
                }

                Self {
                    step,
                    scale_factor,
                    particles: Vec::new(),
                    fields: field_snapshots,
                    n_cells: n,
                    diagnostics: DiagnosticsSummary::zero(scale_factor),
                }
            }
            Content::Mixed { particles, fields } => {
                let positions = (0..particles.count())
                    .map(|p| particles.position_of(p))
                    .collect();
                let momenta = (0..particles.count())
                    .map(|p| particles.momentum_of(p))
                    .collect();

                let particle_snapshots = vec![ParticleSpeciesSnapshot {
                    name: "dark matter particles".to_string(),
                    positions,
                    momenta,
                    mass_particle: particles.mass_particle,
                }];

                let n = fields.grid.n_cells;
                let mut field_snapshots = Vec::new();
                if let Some(ref alpha) = fields.alpha {
                    let density_field = alpha.density(fields.params.mass_alpha);
                    let density: Vec<f64> = density_field.data.iter().copied().collect();
                    field_snapshots.push(FieldSpeciesSnapshot {
                        name: "dark matter field".to_string(),
                        density,
                    });
                }

                Self {
                    step,
                    scale_factor,
                    particles: particle_snapshots,
                    fields: field_snapshots,
                    n_cells: n,
                    diagnostics: DiagnosticsSummary {
                        scale_factor,
                        mass_total: particles.total_mass(),
                        momentum_magnitude: particles.total_momentum().norm(),
                        energy_kinetic: particles.kinetic_energy(scale_factor),
                        energy_potential: 0.0,
                        angular_momentum_magnitude: 0.0,
                    },
                }
            }
        }
    }

    /// Capture a single particle species snapshot.
    fn capture_particles(
        name: &str,
        particles: &Particles,
        step: usize,
        scale_factor: f64,
    ) -> Self {
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
            particles: vec![ParticleSpeciesSnapshot {
                name: name.to_string(),
                positions,
                momenta,
                mass_particle: particles.mass_particle,
            }],
            fields: Vec::new(),
            n_cells: 0,
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
            particles: vec![ParticleSpeciesSnapshot {
                name: "dark_matter".to_string(),
                positions,
                momenta,
                mass_particle: particles.mass_particle,
            }],
            fields: Vec::new(),
            n_cells: 0,
            diagnostics: DiagnosticsSummary::from_diagnostics(&diag),
        }
    }

    /// Access the first particle species' positions (legacy convenience).
    pub fn positions(&self) -> Option<&[Vector<3>]> {
        self.particles.first().map(|s| s.positions.as_slice())
    }

    /// Access the first particle species' momenta (legacy convenience).
    pub fn momenta(&self) -> Option<&[Vector<3>]> {
        self.particles.first().map(|s| s.momenta.as_slice())
    }

    /// Reconstruct a Particles object from the first particle species.
    pub fn to_particles(&self) -> Option<Particles> {
        let species = self.particles.first()?;
        let n = species.positions.len();
        let mut particles = Particles::zeros(n, species.mass_particle);
        for (p, (pos, mom)) in species
            .positions
            .iter()
            .zip(species.momenta.iter())
            .enumerate()
        {
            particles.set_position(p, pos);
            particles.set_momentum(p, mom);
        }

        Some(particles)
    }

    // ========================================================================
    // Serialization
    // ========================================================================

    /// Convert to serializable disk format.
    pub fn to_disk(&self) -> SnapshotOnDisk {
        let particle_species = self
            .particles
            .iter()
            .map(|s| ParticleSpeciesOnDisk {
                name: s.name.clone(),
                positions: s.positions.iter().map(components_from_vector).collect(),
                momenta: s.momenta.iter().map(components_from_vector).collect(),
                mass_particle: s.mass_particle,
            })
            .collect();

        let field_species = self
            .fields
            .iter()
            .map(|s| FieldSpeciesOnDisk {
                name: s.name.clone(),
                density: s.density.clone(),
            })
            .collect();

        SnapshotOnDisk {
            step: self.step,
            scale_factor: self.scale_factor,
            particles: particle_species,
            fields: field_species,
            n_cells: self.n_cells,
            diagnostics: self.diagnostics.clone(),
        }
    }

    /// Reconstruct from disk format.
    pub fn from_disk(disk: SnapshotOnDisk) -> Self {
        let particles = disk
            .particles
            .into_iter()
            .map(|s| ParticleSpeciesSnapshot {
                name: s.name,
                positions: s.positions.iter().map(|c| vector_from_array(*c)).collect(),
                momenta: s.momenta.iter().map(|c| vector_from_array(*c)).collect(),
                mass_particle: s.mass_particle,
            })
            .collect();

        let fields = disk
            .fields
            .into_iter()
            .map(|s| FieldSpeciesSnapshot {
                name: s.name,
                density: s.density,
            })
            .collect();

        Self {
            step: disk.step,
            scale_factor: disk.scale_factor,
            particles,
            fields,
            n_cells: disk.n_cells,
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
    /// Zero diagnostics at a given scale factor.
    pub fn zero(scale_factor: f64) -> Self {
        Self {
            scale_factor,
            mass_total: 0.0,
            momentum_magnitude: 0.0,
            energy_kinetic: 0.0,
            energy_potential: 0.0,
            angular_momentum_magnitude: 0.0,
        }
    }

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
    pub particles: Vec<ParticleSpeciesOnDisk>,
    pub fields: Vec<FieldSpeciesOnDisk>,
    pub n_cells: usize,
    pub diagnostics: DiagnosticsSummary,
}

/// On-disk particle species.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParticleSpeciesOnDisk {
    pub name: String,
    pub positions: Vec<[f64; 3]>,
    pub momenta: Vec<[f64; 3]>,
    pub mass_particle: f64,
}

/// On-disk field species density.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldSpeciesOnDisk {
    pub name: String,
    pub density: Vec<f64>,
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
