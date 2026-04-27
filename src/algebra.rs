//! Shared geometric algebra utilities for hermes.
//!
//! Establishes the Euclidean 3-metric used by all physical quantities in
//! hermes and provides conversion functions at the boundary between morphis
//! geometric objects and flat array storage (CIC kernels, FFT workspace).
//!
//! All positions, momenta, forces, and velocities in hermes live in flat
//! 3D Euclidean space with metric signature (+1, +1, +1).

use morphis::metric::{self, Metric};
use morphis::vector::Vector;
use ndarray::{ArrayD, IxDyn};

/// The Euclidean 3-metric shared by all hermes physical quantities.
pub fn euclidean_3() -> Metric<3> {
    metric::euclidean::<3>()
}

/// Build a grade-1 `Vector<3>` from three Cartesian components.
///
/// For use at the boundary between flat storage and morphis algebra.
/// Not for tight loops — each call allocates an `ArrayD`.
pub fn vector_from_components(x: f64, y: f64, z: f64) -> Vector<3> {
    let mut data = ArrayD::zeros(IxDyn(&[3]));
    data[IxDyn(&[0])] = x;
    data[IxDyn(&[1])] = y;
    data[IxDyn(&[2])] = z;

    Vector::new(data, 1, euclidean_3())
}

/// Build a grade-1 `Vector<3>` from a 3-element array.
pub fn vector_from_array(components: [f64; 3]) -> Vector<3> {
    vector_from_components(components[0], components[1], components[2])
}

/// Extract Cartesian components from a grade-1 `Vector<3>`.
pub fn components_from_vector(v: &Vector<3>) -> [f64; 3] {
    [v.component(&[0]), v.component(&[1]), v.component(&[2])]
}

/// Build a grade-0 scalar `Vector<3>` from an `f64`.
pub fn scalar_from_f64(value: f64) -> Vector<3> {
    Vector::scalar(value, euclidean_3())
}
