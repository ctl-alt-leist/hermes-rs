/// Particle storage in structure-of-arrays layout.
///
/// Positions and momenta are stored as `Array2<f64>` with shape [3, N_p]
/// (component-major) so that each Cartesian component is contiguous in
/// memory. This layout is cache-friendly for the CIC kernel, which
/// processes one dimension at a time.
use ndarray::Array2;

use crate::grid::Grid;

/// Dark matter (or single-species) particle ensemble.
#[derive(Debug, Clone)]
pub struct Particles {
    /// Comoving positions, shape [3, N_p]. Each column is one particle.
    pub position: Array2<f64>,
    /// Canonical momenta p = a² m dx/dt, shape [3, N_p].
    pub momentum: Array2<f64>,
    /// Mass per particle (M_☉). Fixed for the species.
    pub mass_particle: f64,
}

impl Particles {
    /// Create particles with zero position and momentum.
    pub fn zeros(n_particles: usize, mass_particle: f64) -> Self {
        Self {
            position: Array2::zeros((3, n_particles)),
            momentum: Array2::zeros((3, n_particles)),
            mass_particle,
        }
    }

    /// Number of particles.
    pub fn count(&self) -> usize {
        self.position.ncols()
    }

    /// Total mass of all particles (M_☉).
    pub fn total_mass(&self) -> f64 {
        self.count() as f64 * self.mass_particle
    }

    /// Get the position of particle `n` as [x, y, z].
    pub fn position_of(&self, n: usize) -> [f64; 3] {
        [
            self.position[[0, n]],
            self.position[[1, n]],
            self.position[[2, n]],
        ]
    }

    /// Get the momentum of particle `n` as [px, py, pz].
    pub fn momentum_of(&self, n: usize) -> [f64; 3] {
        [
            self.momentum[[0, n]],
            self.momentum[[1, n]],
            self.momentum[[2, n]],
        ]
    }

    /// Set the position of particle `n`.
    pub fn set_position(&mut self, n: usize, pos: [f64; 3]) {
        self.position[[0, n]] = pos[0];
        self.position[[1, n]] = pos[1];
        self.position[[2, n]] = pos[2];
    }

    /// Set the momentum of particle `n`.
    pub fn set_momentum(&mut self, n: usize, mom: [f64; 3]) {
        self.momentum[[0, n]] = mom[0];
        self.momentum[[1, n]] = mom[1];
        self.momentum[[2, n]] = mom[2];
    }

    /// Wrap all particle positions into the periodic box [0, box_length).
    pub fn wrap_positions(&mut self, grid: &Grid) {
        for n in 0..self.count() {
            for d in 0..3 {
                self.position[[d, n]] = grid.wrap_position(self.position[[d, n]]);
            }
        }
    }

    /// Place particles on a uniform lattice with `n_per_side`³ particles.
    ///
    /// Positions are at the center of each Lagrangian cell. The particle
    /// mass is set from the total mass budget: mass_particle = density_mean × box_volume / N_p.
    pub fn on_lattice(n_per_side: usize, grid: &Grid, density_mean: f64) -> Self {
        let n_total = n_per_side * n_per_side * n_per_side;
        let mass_particle = density_mean * grid.box_volume() / n_total as f64;
        let spacing = grid.box_length / n_per_side as f64;

        let mut particles = Self::zeros(n_total, mass_particle);

        let mut n = 0;
        for m0 in 0..n_per_side {
            for m1 in 0..n_per_side {
                for m2 in 0..n_per_side {
                    particles.position[[0, n]] = (m0 as f64 + 0.5) * spacing;
                    particles.position[[1, n]] = (m1 as f64 + 0.5) * spacing;
                    particles.position[[2, n]] = (m2 as f64 + 0.5) * spacing;
                    n += 1;
                }
            }
        }

        particles
    }
}
