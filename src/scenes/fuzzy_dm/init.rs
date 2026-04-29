//! Fuzzy dark matter initial conditions.
//!
//! Initializes a Gaussian wavepacket in the even subalgebra:
//! psi(x) = sqrt(rho(x) / m) * exp(I * S(x) / ell)
//!
//! The density profile is a Gaussian centered in the box, and the
//! phase S encodes zero initial velocity (stationary soliton-like IC).

use morphis::even_field::EvenField;
use morphis::metric;

use crate::core::content::FieldParams;
use crate::physics::grid::Grid as HermesGrid;

/// Create a Gaussian wavepacket centered in the box.
///
/// The density profile is:
///   rho(x) = rho_0 * exp(-|x - x_center|^2 / (2 sigma^2))
///
/// The wavefunction is:
///   psi = sqrt(rho / m) (scalar part only, zero phase = stationary)
pub fn gaussian_wavepacket(hermes_grid: &HermesGrid, params: &FieldParams) -> EvenField<3> {
    let n = hermes_grid.n_cells;
    let box_length = hermes_grid.box_length;
    let _cell_length = hermes_grid.cell_length;
    let center = box_length / 2.0;
    let sigma = box_length * 0.1; // width = 10% of box
    let mass = params.mass_alpha;

    let morphis_grid = morphis::grid::Grid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();

    EvenField::from_fn(&morphis_grid, g, |x| {
        let dx = x[0] - center;
        let dy = x[1] - center;
        let dz = x[2] - center;
        let r2 = dx * dx + dy * dy + dz * dz;

        let rho = (-r2 / (2.0 * sigma * sigma)).exp();
        let amplitude = (rho / mass).sqrt();

        // (scalar, pseudoscalar) = (amplitude, 0) — zero phase, stationary
        (amplitude, 0.0)
    })
}
