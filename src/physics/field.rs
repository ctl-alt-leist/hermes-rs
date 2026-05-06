//! Grid-based scalar and vector fields as graded geometric objects.
//!
//! Fields are derived (per-step) quantities living on the grid — density,
//! potential, force — not persistent state. Each field carries a morphis
//! `Metric<3>` and provides access to its values as morphis vectors of
//! the appropriate grade.

use std::ops::{Add, Mul, Sub};

use morphis::metric::Metric;
use morphis::ops::geometric;
use morphis::vector::Vector;
use ndarray::Array3;

use crate::algebra::{euclidean_3, scalar_from_f64, vector_from_components};
use crate::physics::grid::Grid;

// ============================================================================
// ScalarField — grade-0 field
// ============================================================================

/// A grade-0 (scalar) field on the grid (e.g. density, potential, energy).
///
/// Each value is a morphis grade-0 scalar in the Euclidean 3-metric.
#[derive(Debug, Clone)]
pub struct ScalarField {
    pub data: Array3<f64>,
    /// Metric defining the geometric context.
    pub metric: Metric<3>,
}

impl ScalarField {
    /// Create a zero-valued scalar field for the given grid.
    pub fn zeros(grid: &Grid) -> Self {
        let n = grid.n_cells;

        Self {
            data: Array3::zeros((n, n, n)),
            metric: euclidean_3(),
        }
    }

    /// Create a scalar field from an existing array.
    pub fn from_array(data: Array3<f64>) -> Self {
        Self {
            data,
            metric: euclidean_3(),
        }
    }

    /// Total sum of all values.
    pub fn sum(&self) -> f64 {
        self.data.sum()
    }

    /// Scale all values by a constant factor.
    pub fn scale(&mut self, factor: f64) {
        self.data *= factor;
    }

    /// Raw value at cell (m0, m1, m2).
    pub fn get(&self, m0: usize, m1: usize, m2: usize) -> f64 {
        self.data[[m0, m1, m2]]
    }

    /// Mutable raw access at cell (m0, m1, m2).
    pub fn get_mut(&mut self, m0: usize, m1: usize, m2: usize) -> &mut f64 {
        &mut self.data[[m0, m1, m2]]
    }

    /// Extract the value at a cell as a morphis grade-0 scalar.
    pub fn scalar_at(&self, m0: usize, m1: usize, m2: usize) -> Vector<3> {
        scalar_from_f64(self.data[[m0, m1, m2]])
    }
}

impl Add for &ScalarField {
    type Output = ScalarField;

    fn add(self, other: &ScalarField) -> ScalarField {
        ScalarField {
            data: &self.data + &other.data,
            metric: self.metric,
        }
    }
}

impl Sub for &ScalarField {
    type Output = ScalarField;

    fn sub(self, other: &ScalarField) -> ScalarField {
        ScalarField {
            data: &self.data - &other.data,
            metric: self.metric,
        }
    }
}

impl Mul<f64> for &ScalarField {
    type Output = ScalarField;

    fn mul(self, scalar: f64) -> ScalarField {
        ScalarField {
            data: &self.data * scalar,
            metric: self.metric,
        }
    }
}

// ============================================================================
// VectorField — grade-1 field
// ============================================================================

/// A grade-1 (vector) field on the grid (e.g. force, velocity, momentum density).
///
/// Each value is a morphis grade-1 vector in the Euclidean 3-metric.
/// Internal storage is three component arrays for FFT compatibility.
#[derive(Debug, Clone)]
pub struct VectorField {
    pub data: [Array3<f64>; 3],
    /// Metric defining the geometric context.
    pub metric: Metric<3>,
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
            metric: euclidean_3(),
        }
    }

    /// Access one component array (for FFT/CIC kernels).
    pub fn component(&self, d: usize) -> &Array3<f64> {
        &self.data[d]
    }

    /// Mutable access to one component array (for FFT/CIC kernels).
    pub fn component_mut(&mut self, d: usize) -> &mut Array3<f64> {
        &mut self.data[d]
    }

    /// Raw component access at a cell.
    pub fn get(&self, d: usize, m0: usize, m1: usize, m2: usize) -> f64 {
        self.data[d][[m0, m1, m2]]
    }

    /// Scale all components by a constant factor.
    pub fn scale(&mut self, factor: f64) {
        for comp in &mut self.data {
            *comp *= factor;
        }
    }

    /// Extract the vector at a cell as a morphis grade-1 vector.
    pub fn vector_at(&self, m0: usize, m1: usize, m2: usize) -> Vector<3> {
        vector_from_components(
            self.data[0][[m0, m1, m2]],
            self.data[1][[m0, m1, m2]],
            self.data[2][[m0, m1, m2]],
        )
    }

    /// Write a morphis grade-1 vector into a cell.
    pub fn set_vector_at(&mut self, m0: usize, m1: usize, m2: usize, v: &Vector<3>) {
        self.data[0][[m0, m1, m2]] = v.component(&[1]);
        self.data[1][[m0, m1, m2]] = v.component(&[2]);
        self.data[2][[m0, m1, m2]] = v.component(&[3]);
    }

    /// Dot product of two vector fields at a cell, via the geometric product.
    ///
    /// Returns the scalar part of u * v, which for grade-1 vectors is
    /// the metric inner product g(u, v).
    pub fn dot_at(&self, m0: usize, m1: usize, m2: usize, other: &VectorField) -> f64 {
        let u = self.vector_at(m0, m1, m2);
        let v = other.vector_at(m0, m1, m2);

        geometric(&u, &v).scalar_part()
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
            metric: self.metric,
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
            metric: self.metric,
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
            metric: self.metric,
        }
    }
}
