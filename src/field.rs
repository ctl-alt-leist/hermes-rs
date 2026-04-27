/// Grid-based scalar and vector fields.
///
/// These are derived (per-step) quantities living on the grid — density,
/// potential, force — not persistent state. Thin wrappers over ndarray
/// providing field-specific semantics.
use std::ops::{Add, Mul, Sub};

use ndarray::Array3;

use crate::grid::Grid;

// ============================================================================
// ScalarField
// ============================================================================

/// A scalar field on the grid (e.g. density, potential, energy density).
#[derive(Debug, Clone)]
pub struct ScalarField {
    pub data: Array3<f64>,
}

impl ScalarField {
    /// Create a zero-valued scalar field for the given grid.
    pub fn zeros(grid: &Grid) -> Self {
        let n = grid.n_cells;

        Self {
            data: Array3::zeros((n, n, n)),
        }
    }

    /// Create a scalar field from an existing array.
    pub fn from_array(data: Array3<f64>) -> Self {
        Self { data }
    }

    /// Total sum of all values (useful for mass conservation checks).
    pub fn sum(&self) -> f64 {
        self.data.sum()
    }

    /// Scale all values by a constant factor.
    pub fn scale(&mut self, factor: f64) {
        self.data *= factor;
    }

    /// Element-wise access.
    pub fn get(&self, m0: usize, m1: usize, m2: usize) -> f64 {
        self.data[[m0, m1, m2]]
    }

    /// Mutable element-wise access.
    pub fn get_mut(&mut self, m0: usize, m1: usize, m2: usize) -> &mut f64 {
        &mut self.data[[m0, m1, m2]]
    }
}

impl Add for &ScalarField {
    type Output = ScalarField;

    fn add(self, other: &ScalarField) -> ScalarField {
        ScalarField {
            data: &self.data + &other.data,
        }
    }
}

impl Sub for &ScalarField {
    type Output = ScalarField;

    fn sub(self, other: &ScalarField) -> ScalarField {
        ScalarField {
            data: &self.data - &other.data,
        }
    }
}

impl Mul<f64> for &ScalarField {
    type Output = ScalarField;

    fn mul(self, scalar: f64) -> ScalarField {
        ScalarField {
            data: &self.data * scalar,
        }
    }
}

// ============================================================================
// VectorField
// ============================================================================

/// A 3D vector field on the grid (e.g. force, velocity, momentum density).
///
/// Stored as three separate scalar arrays, one per Cartesian component.
#[derive(Debug, Clone)]
pub struct VectorField {
    pub data: [Array3<f64>; 3],
}

impl VectorField {
    /// Create a zero-valued vector field for the given grid.
    pub fn zeros(grid: &Grid) -> Self {
        let n = grid.n_cells;

        Self {
            data: [
                Array3::zeros((n, n, n)),
                Array3::zeros((n, n, n)),
                Array3::zeros((n, n, n)),
            ],
        }
    }

    /// Access one Cartesian component (0, 1, or 2).
    pub fn component(&self, d: usize) -> &Array3<f64> {
        &self.data[d]
    }

    /// Mutable access to one Cartesian component.
    pub fn component_mut(&mut self, d: usize) -> &mut Array3<f64> {
        &mut self.data[d]
    }

    /// Element-wise access: the d-th component at cell (m0, m1, m2).
    pub fn get(&self, d: usize, m0: usize, m1: usize, m2: usize) -> f64 {
        self.data[d][[m0, m1, m2]]
    }

    /// Scale all components by a constant factor.
    pub fn scale(&mut self, factor: f64) {
        for component in &mut self.data {
            *component *= factor;
        }
    }
}

impl Add for &VectorField {
    type Output = VectorField;

    fn add(self, other: &VectorField) -> VectorField {
        VectorField {
            data: [
                &self.data[0] + &other.data[0],
                &self.data[1] + &other.data[1],
                &self.data[2] + &other.data[2],
            ],
        }
    }
}

impl Sub for &VectorField {
    type Output = VectorField;

    fn sub(self, other: &VectorField) -> VectorField {
        VectorField {
            data: [
                &self.data[0] - &other.data[0],
                &self.data[1] - &other.data[1],
                &self.data[2] - &other.data[2],
            ],
        }
    }
}

impl Mul<f64> for &VectorField {
    type Output = VectorField;

    fn mul(self, scalar: f64) -> VectorField {
        VectorField {
            data: [
                &self.data[0] * scalar,
                &self.data[1] * scalar,
                &self.data[2] * scalar,
            ],
        }
    }
}
