use std::f64::consts::PI;

use hermes_rs::field::ScalarField;
use hermes_rs::grid::Grid;
use hermes_rs::poisson::PoissonSolver;

const DENSITY_MEAN: f64 = 1e-7; // M_☉ / kpc³
const SCALE_FACTOR: f64 = 1.0;

// ============================================================================
// Zero force from uniform density
// ============================================================================

#[test]
fn uniform_density_gives_zero_force() {
    let grid = Grid::new(16, 100.0);
    let mut solver = PoissonSolver::new(&grid);

    // Overdensity δ = 0 everywhere (uniform density).
    let overdensity = ScalarField::zeros(&grid);
    let force = solver.solve(&overdensity, DENSITY_MEAN, SCALE_FACTOR);

    for d in 0..3 {
        let max_force = force.data[d]
            .iter()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_force < 1e-10,
            "uniform density should give zero force, component {d} has max {max_force}"
        );
    }
}

// ============================================================================
// Single sinusoidal mode
// ============================================================================

#[test]
fn single_mode_potential() {
    let n = 32;
    let box_length = 100.0;
    let grid = Grid::new(n, box_length);
    let mut solver = PoissonSolver::new(&grid);

    // Place a single Fourier mode: δ(x) = A sin(2π x / L).
    // The potential should be: ϕ(x) = -(4πG ρ̄ a²) A / k² sin(2π x / L)
    // and the x-force: Fx = -dϕ/dx = (4πG ρ̄ a²) A k / k² cos(...)
    //                               = (4πG ρ̄ a²) A / k cos(2π x / L)
    let amplitude = 0.1;
    let k_fund = 2.0 * PI / box_length;
    let h = grid.cell_length;

    let mut overdensity = ScalarField::zeros(&grid);
    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * h;
        let overdensity_x = amplitude * (k_fund * x).sin();
        for m1 in 0..n {
            for m2 in 0..n {
                *overdensity.get_mut(m0, m1, m2) = overdensity_x;
            }
        }
    }

    let force = solver.solve(&overdensity, DENSITY_MEAN, SCALE_FACTOR);

    // The force in x should follow cos(2π x / L), and y,z forces should be ~0.
    // Use the discrete k² for the expected amplitude.
    let kx_h_half = PI / n as f64;
    let k2_discrete = (2.0 / h) * (2.0 / h) * kx_h_half.sin().powi(2);
    let prefactor = 4.0 * PI * hermes_rs::constants::G * DENSITY_MEAN * SCALE_FACTOR * SCALE_FACTOR;
    let expected_amplitude = prefactor * amplitude * k_fund / k2_discrete;

    // Check Fx at a few points.
    let m1 = n / 2;
    let m2 = n / 2;
    let mut max_err: f64 = 0.0;
    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * h;
        let expected = expected_amplitude * (k_fund * x).cos();
        let computed = force.get(0, m0, m1, m2);
        let err = (computed - expected).abs();
        max_err = max_err.max(err);
    }

    let rel_err = max_err / expected_amplitude;
    assert!(
        rel_err < 1e-4,
        "single mode Fx: relative error {rel_err} exceeds tolerance"
    );

    // Fy and Fz should be negligible.
    for d in 1..3 {
        let max_force = force.data[d]
            .iter()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_force < expected_amplitude * 1e-10,
            "F{} should be ~0 for x-only mode, max = {max_force}",
            ["x", "y", "z"][d]
        );
    }
}

// ============================================================================
// Isotropy: same mode along each axis produces same force magnitude
// ============================================================================

#[test]
fn force_isotropic() {
    let n = 16;
    let box_length = 100.0;
    let grid = Grid::new(n, box_length);
    let h = grid.cell_length;
    let k_fund = 2.0 * PI / box_length;
    let amplitude = 0.05;

    let mut amplitudes_by_axis = [0.0_f64; 3];

    for axis in 0..3 {
        let mut solver = PoissonSolver::new(&grid);
        let mut overdensity = ScalarField::zeros(&grid);

        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    let coord = match axis {
                        0 => (m0 as f64 + 0.5) * h,
                        1 => (m1 as f64 + 0.5) * h,
                        _ => (m2 as f64 + 0.5) * h,
                    };
                    *overdensity.get_mut(m0, m1, m2) = amplitude * (k_fund * coord).sin();
                }
            }
        }

        let force = solver.solve(&overdensity, DENSITY_MEAN, SCALE_FACTOR);

        // Max absolute force in the direction of the mode.
        amplitudes_by_axis[axis] = force.data[axis]
            .iter()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
    }

    // All three axes should produce the same peak force amplitude.
    let mean = amplitudes_by_axis.iter().sum::<f64>() / 3.0;
    for (axis, &amp) in amplitudes_by_axis.iter().enumerate() {
        let rel_err = (amp - mean).abs() / mean;
        assert!(
            rel_err < 1e-10,
            "axis {axis} amplitude {amp} differs from mean {mean} by {rel_err}"
        );
    }
}

// ============================================================================
// Force antisymmetry: F(-δ) = -F(δ)
// ============================================================================

#[test]
fn force_antisymmetric() {
    let n = 16;
    let box_length = 100.0;
    let grid = Grid::new(n, box_length);
    let h = grid.cell_length;
    let k_fund = 2.0 * PI / box_length;

    let mut solver = PoissonSolver::new(&grid);

    let mut overdensity_pos = ScalarField::zeros(&grid);
    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * h;
        let overdensity_x = 0.1 * (k_fund * x).sin();
        for m1 in 0..n {
            for m2 in 0..n {
                *overdensity_pos.get_mut(m0, m1, m2) = overdensity_x;
            }
        }
    }

    let overdensity_neg = &overdensity_pos * -1.0;

    let force_pos = solver.solve(&overdensity_pos, DENSITY_MEAN, SCALE_FACTOR);
    let force_neg = solver.solve(&overdensity_neg, DENSITY_MEAN, SCALE_FACTOR);

    // F(-δ) should equal -F(δ).
    for d in 0..3 {
        for ((m0, m1, m2), &force_cell_pos) in force_pos.data[d].indexed_iter() {
            let force_cell_neg = force_neg.data[d][[m0, m1, m2]];
            let err = (force_cell_pos + force_cell_neg).abs();
            assert!(
                err < 1e-10,
                "F{d}[{m0},{m1},{m2}]: F(+δ) = {force_cell_pos}, F(-δ) = {force_cell_neg}, sum = {err}"
            );
        }
    }
}

// ============================================================================
// Zero mode is projected out
// ============================================================================

#[test]
fn nonzero_mean_overdensity_still_works() {
    let grid = Grid::new(16, 100.0);
    let mut solver = PoissonSolver::new(&grid);
    let h = grid.cell_length;
    let k_fund = 2.0 * PI / grid.box_length;

    // δ = 0.5 + 0.1 sin(2π x / L).
    // The constant part (0.5) has zero wavevector and should be projected
    // out by G(0,0,0) = 0, leaving only the sinusoidal response.
    let mut overdensity = ScalarField::zeros(&grid);
    for m0 in 0..grid.n_cells {
        let x = (m0 as f64 + 0.5) * h;
        let overdensity_x = 0.5 + 0.1 * (k_fund * x).sin();
        for m1 in 0..grid.n_cells {
            for m2 in 0..grid.n_cells {
                *overdensity.get_mut(m0, m1, m2) = overdensity_x;
            }
        }
    }

    let force = solver.solve(&overdensity, DENSITY_MEAN, SCALE_FACTOR);

    // Fx should be nonzero (from the sine), Fy and Fz should be ~0.
    let max_fx = force.data[0]
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    let max_fy = force.data[1]
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    let max_fz = force.data[2]
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);

    assert!(max_fx > 1e-15, "Fx should be nonzero from the sine mode");
    assert!(max_fy < 1e-10, "Fy should be ~0, got {max_fy}");
    assert!(max_fz < 1e-10, "Fz should be ~0, got {max_fz}");
}
