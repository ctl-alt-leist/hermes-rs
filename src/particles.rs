//! Particle storage with morphis geometric algebra interface.
//!
//! Positions and momenta are geometric objects — grade-1 vectors in the
//! Euclidean 3-metric. The primary interface returns morphis `Vector<3>`
//! values and uses morphis operations (norms, wedge products) for all
//! derived quantities.
//!
//! Internal storage is `Array2<f64>` with shape [3, N_p] for cache-friendly
//! access in CIC kernels. The `_components` methods provide raw access at
//! the flat-array boundary; all other access is through morphis.

use morphis::ops::wedge;
use morphis::vector::Vector;
use ndarray::Array2;

use crate::algebra::{components_from_vector, euclidean_3, vector_from_components};
use crate::grid::Grid;

/// Dark matter (or single-species) particle ensemble.
#[derive(Debug, Clone)]
pub struct Particles {
    /// Comoving positions, shape [3, N_p].
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

    // ========================================================================
    // Morphis-native access — primary interface
    // ========================================================================

    /// Position of particle `n` as a morphis grade-1 vector.
    pub fn position_of(&self, n: usize) -> Vector<3> {
        vector_from_components(
            self.position[[0, n]],
            self.position[[1, n]],
            self.position[[2, n]],
        )
    }

    /// Momentum of particle `n` as a morphis grade-1 vector.
    pub fn momentum_of(&self, n: usize) -> Vector<3> {
        vector_from_components(
            self.momentum[[0, n]],
            self.momentum[[1, n]],
            self.momentum[[2, n]],
        )
    }

    /// Set the position of particle `n` from a morphis grade-1 vector.
    pub fn set_position(&mut self, n: usize, v: &Vector<3>) {
        let c = components_from_vector(v);
        self.position[[0, n]] = c[0];
        self.position[[1, n]] = c[1];
        self.position[[2, n]] = c[2];
    }

    /// Set the momentum of particle `n` from a morphis grade-1 vector.
    pub fn set_momentum(&mut self, n: usize, v: &Vector<3>) {
        let c = components_from_vector(v);
        self.momentum[[0, n]] = c[0];
        self.momentum[[1, n]] = c[1];
        self.momentum[[2, n]] = c[2];
    }

    // ========================================================================
    // Fast-path component access — for CIC/FFT hot loops only
    // ========================================================================

    /// Position of particle `n` as raw Cartesian components.
    pub fn position_components(&self, n: usize) -> [f64; 3] {
        [
            self.position[[0, n]],
            self.position[[1, n]],
            self.position[[2, n]],
        ]
    }

    /// Momentum of particle `n` as raw Cartesian components.
    pub fn momentum_components(&self, n: usize) -> [f64; 3] {
        [
            self.momentum[[0, n]],
            self.momentum[[1, n]],
            self.momentum[[2, n]],
        ]
    }

    /// Set position of particle `n` from raw Cartesian components.
    pub fn set_position_components(&mut self, n: usize, pos: [f64; 3]) {
        self.position[[0, n]] = pos[0];
        self.position[[1, n]] = pos[1];
        self.position[[2, n]] = pos[2];
    }

    /// Set momentum of particle `n` from raw Cartesian components.
    pub fn set_momentum_components(&mut self, n: usize, mom: [f64; 3]) {
        self.momentum[[0, n]] = mom[0];
        self.momentum[[1, n]] = mom[1];
        self.momentum[[2, n]] = mom[2];
    }

    // ========================================================================
    // Morphis-native derived quantities
    // ========================================================================

    /// Total momentum as a morphis grade-1 vector: P = Σ p_n.
    pub fn total_momentum(&self) -> Vector<3> {
        let mut total = Vector::<3>::zero(1, euclidean_3());
        for n in 0..self.count() {
            let p = self.momentum_of(n);
            total = &total + &p;
        }

        total
    }

    /// Angular momentum bivector of particle `n`: L_n = x_n ∧ p_n.
    ///
    /// Returns a grade-2 morphis vector (bivector) encoding the oriented
    /// plane of rotation. This is the natural representation — not a
    /// pseudovector cross product.
    pub fn angular_momentum(&self, n: usize) -> Vector<3> {
        let x = self.position_of(n);
        let p = self.momentum_of(n);

        wedge(&x, &p)
    }

    /// Total angular momentum bivector: L = Σ x_n ∧ p_n.
    pub fn total_angular_momentum(&self) -> Vector<3> {
        let mut total = Vector::<3>::zero(2, euclidean_3());
        for n in 0..self.count() {
            let angular_momentum_n = self.angular_momentum(n);
            total = &total + &angular_momentum_n;
        }

        total
    }

    /// Total kinetic energy: E_k = Σ |p_n|² / (2 m a²).
    ///
    /// Uses the morphis norm (metric-aware) to compute |p|².
    pub fn kinetic_energy(&self, scale_factor: f64) -> f64 {
        let denominator = 2.0 * self.mass_particle * scale_factor * scale_factor;
        let mut energy = 0.0;
        for n in 0..self.count() {
            let p = self.momentum_of(n);
            energy += p.norm_squared() / denominator;
        }

        energy
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
                    particles.set_position_components(
                        n,
                        [
                            (m0 as f64 + 0.5) * spacing,
                            (m1 as f64 + 0.5) * spacing,
                            (m2 as f64 + 0.5) * spacing,
                        ],
                    );
                    n += 1;
                }
            }
        }

        particles
    }
}
