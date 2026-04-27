//! Cloud-in-Cell (CIC) mass assignment and force interpolation.
//!
//! CIC uses a trilinear (8-cell) kernel to distribute particle masses onto
//! the grid and to interpolate grid forces back to particle positions. Using
//! the same kernel for both operations is the discrete analog of Newton's
//! third law and ensures momentum conservation to machine precision.
//!
//! The CIC inner loops operate on flat arrays for cache performance. The
//! morphis boundary lives at the public interface: `assign_density` returns
//! a metric-aware `ScalarField`, and `interpolate_force` returns
//! `ParticleForces` with morphis vector access.

use morphis::vector::Vector;
use ndarray::Array2;

use crate::algebra::{euclidean_3, vector_from_components};
use crate::field::{ScalarField, VectorField};
use crate::grid::Grid;
use crate::particles::Particles;

// ============================================================================
// ParticleForces — morphis-aware per-particle force result
// ============================================================================

/// Per-particle force vectors from CIC interpolation.
///
/// Internal storage is flat `Array2<f64>` for hot-path access.
/// The primary interface returns morphis grade-1 vectors.
#[derive(Debug, Clone)]
pub struct ParticleForces {
    /// Raw force components, shape [3, N_p].
    pub data: Array2<f64>,
}

impl ParticleForces {
    /// Number of particles.
    pub fn count(&self) -> usize {
        self.data.ncols()
    }

    /// Force on particle `n` as a morphis grade-1 vector.
    pub fn force_on(&self, n: usize) -> Vector<3> {
        vector_from_components(self.data[[0, n]], self.data[[1, n]], self.data[[2, n]])
    }

    /// Force on particle `n` as raw Cartesian components (CIC/integrator hot path).
    pub fn force_components(&self, n: usize) -> [f64; 3] {
        [self.data[[0, n]], self.data[[1, n]], self.data[[2, n]]]
    }

    /// Total force summed over all particles as a morphis grade-1 vector.
    ///
    /// Should vanish for a self-consistent periodic system (Newton's third law).
    pub fn total_force(&self) -> Vector<3> {
        let mut total = Vector::<3>::zero(1, euclidean_3());
        for n in 0..self.count() {
            let force_n = self.force_on(n);
            total = &total + &force_n;
        }

        total
    }
}

// ============================================================================
// Mass assignment
// ============================================================================

/// Deposit particle masses onto the grid using CIC interpolation.
///
/// Returns a grade-0 density field ρ(x) in M_☉ / kpc³. Each particle's
/// mass is distributed across the 8 surrounding cells using trilinear weights.
pub fn assign_density(particles: &Particles, grid: &Grid) -> ScalarField {
    let mut field = ScalarField::zeros(grid);
    let n = grid.n_cells;
    let h = grid.cell_length;
    let h_inv = 1.0 / h;

    for p in 0..particles.count() {
        let pos = particles.position_components(p);
        let cell = [
            pos[0] * h_inv - 0.5,
            pos[1] * h_inv - 0.5,
            pos[2] * h_inv - 0.5,
        ];

        let base = [
            cell[0].floor() as isize,
            cell[1].floor() as isize,
            cell[2].floor() as isize,
        ];

        let frac = [
            cell[0] - base[0] as f64,
            cell[1] - base[1] as f64,
            cell[2] - base[2] as f64,
        ];

        let weight = [
            [1.0 - frac[0], frac[0]],
            [1.0 - frac[1], frac[1]],
            [1.0 - frac[2], frac[2]],
        ];

        for (a, &weight_a) in weight[0].iter().enumerate() {
            let g0 = ((base[0] + a as isize) % n as isize + n as isize) as usize % n;
            for (b, &weight_b) in weight[1].iter().enumerate() {
                let g1 = ((base[1] + b as isize) % n as isize + n as isize) as usize % n;
                for (c, &weight_c) in weight[2].iter().enumerate() {
                    let g2 = ((base[2] + c as isize) % n as isize + n as isize) as usize % n;
                    field.data[[g0, g1, g2]] +=
                        particles.mass_particle * weight_a * weight_b * weight_c;
                }
            }
        }
    }

    let volume_inv = 1.0 / grid.cell_volume();
    field.data *= volume_inv;

    field
}

// ============================================================================
// Force interpolation
// ============================================================================

/// Interpolate a vector force field at particle positions using CIC weights.
///
/// Returns morphis-aware `ParticleForces`. Uses the same trilinear kernel
/// as `assign_density` to ensure momentum conservation by kernel symmetry.
pub fn interpolate_force(
    force: &VectorField,
    particles: &Particles,
    grid: &Grid,
) -> ParticleForces {
    let n_particles = particles.count();
    let n = grid.n_cells;
    let h = grid.cell_length;
    let h_inv = 1.0 / h;
    let mut result = Array2::zeros((3, n_particles));

    for p in 0..n_particles {
        let pos = particles.position_components(p);
        let cell = [
            pos[0] * h_inv - 0.5,
            pos[1] * h_inv - 0.5,
            pos[2] * h_inv - 0.5,
        ];

        let base = [
            cell[0].floor() as isize,
            cell[1].floor() as isize,
            cell[2].floor() as isize,
        ];

        let frac = [
            cell[0] - base[0] as f64,
            cell[1] - base[1] as f64,
            cell[2] - base[2] as f64,
        ];

        let weight = [
            [1.0 - frac[0], frac[0]],
            [1.0 - frac[1], frac[1]],
            [1.0 - frac[2], frac[2]],
        ];

        for (a, &weight_a) in weight[0].iter().enumerate() {
            let g0 = ((base[0] + a as isize) % n as isize + n as isize) as usize % n;
            for (b, &weight_b) in weight[1].iter().enumerate() {
                let g1 = ((base[1] + b as isize) % n as isize + n as isize) as usize % n;
                for (c, &weight_c) in weight[2].iter().enumerate() {
                    let g2 = ((base[2] + c as isize) % n as isize + n as isize) as usize % n;
                    let w = weight_a * weight_b * weight_c;
                    for d in 0..3 {
                        result[[d, p]] += force.data[d][[g0, g1, g2]] * w;
                    }
                }
            }
        }
    }

    ParticleForces { data: result }
}
