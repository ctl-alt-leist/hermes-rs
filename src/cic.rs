//! Cloud-in-Cell (CIC) mass assignment and force interpolation.
//!
//! CIC uses a trilinear (8-cell) kernel to distribute particle masses onto
//! the grid and to interpolate grid forces back to particle positions. Using
//! the same kernel for both operations is the discrete analog of Newton's
//! third law and ensures momentum conservation to machine precision.

use ndarray::Array2;

use crate::field::{ScalarField, VectorField};
use crate::grid::Grid;
use crate::particles::Particles;

/// Deposit particle masses onto the grid using CIC interpolation.
///
/// Returns the density field ρ(x) in M_☉ / kpc³. Each particle's mass is
/// distributed across the 8 surrounding cells using trilinear weights.
pub fn assign_density(particles: &Particles, grid: &Grid) -> ScalarField {
    let mut field = ScalarField::zeros(grid);
    let n = grid.n_cells;
    let h = grid.cell_length;
    let h_inv = 1.0 / h;

    for p in 0..particles.count() {
        let x = particles.position[[0, p]];
        let y = particles.position[[1, p]];
        let z = particles.position[[2, p]];

        // Cell index of the lower-left corner of the CIC stencil.
        // Particle at cell center (i+0.5)*h has equal weight to cell i and i+1.
        let xc = x * h_inv - 0.5;
        let yc = y * h_inv - 0.5;
        let zc = z * h_inv - 0.5;

        let m0 = xc.floor() as isize;
        let m1 = yc.floor() as isize;
        let m2 = zc.floor() as isize;

        let dx = xc - m0 as f64;
        let dy = yc - m1 as f64;
        let dz = zc - m2 as f64;

        let wx = [1.0 - dx, dx];
        let wy = [1.0 - dy, dy];
        let wz = [1.0 - dz, dz];

        for (a, &weight_x) in wx.iter().enumerate() {
            let g0 = ((m0 + a as isize) % n as isize + n as isize) as usize % n;
            for (b, &weight_y) in wy.iter().enumerate() {
                let g1 = ((m1 + b as isize) % n as isize + n as isize) as usize % n;
                for (c, &weight_z) in wz.iter().enumerate() {
                    let g2 = ((m2 + c as isize) % n as isize + n as isize) as usize % n;
                    let weight = weight_x * weight_y * weight_z;
                    field.data[[g0, g1, g2]] += particles.mass_particle * weight;
                }
            }
        }
    }

    // Convert from mass per cell to density (M_☉ / kpc³).
    let volume_inv = 1.0 / grid.cell_volume();
    field.data *= volume_inv;

    field
}

/// Interpolate a vector force field at particle positions using CIC weights.
///
/// Returns per-particle forces as `Array2<f64>` with shape [3, N_p].
/// Uses the same trilinear kernel as `assign_density` to ensure
/// momentum conservation by kernel symmetry (Newton's third law).
pub fn interpolate_force(force: &VectorField, particles: &Particles, grid: &Grid) -> Array2<f64> {
    let n_particles = particles.count();
    let n = grid.n_cells;
    let h = grid.cell_length;
    let h_inv = 1.0 / h;
    let mut result = Array2::zeros((3, n_particles));

    for p in 0..n_particles {
        let x = particles.position[[0, p]];
        let y = particles.position[[1, p]];
        let z = particles.position[[2, p]];

        let xc = x * h_inv - 0.5;
        let yc = y * h_inv - 0.5;
        let zc = z * h_inv - 0.5;

        let m0 = xc.floor() as isize;
        let m1 = yc.floor() as isize;
        let m2 = zc.floor() as isize;

        let dx = xc - m0 as f64;
        let dy = yc - m1 as f64;
        let dz = zc - m2 as f64;

        let wx = [1.0 - dx, dx];
        let wy = [1.0 - dy, dy];
        let wz = [1.0 - dz, dz];

        for (a, &weight_x) in wx.iter().enumerate() {
            let g0 = ((m0 + a as isize) % n as isize + n as isize) as usize % n;
            for (b, &weight_y) in wy.iter().enumerate() {
                let g1 = ((m1 + b as isize) % n as isize + n as isize) as usize % n;
                for (c, &weight_z) in wz.iter().enumerate() {
                    let g2 = ((m2 + c as isize) % n as isize + n as isize) as usize % n;
                    let weight = weight_x * weight_y * weight_z;
                    for d in 0..3 {
                        result[[d, p]] += force.data[d][[g0, g1, g2]] * weight;
                    }
                }
            }
        }
    }

    result
}
