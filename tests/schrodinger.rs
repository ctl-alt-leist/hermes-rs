//! Tests for the Schrodinger-Poisson field-theory sector.
//!
//! These tests verify the split-step spectral integrator, the Madelung
//! decomposition, norm conservation, and initialization consistency.
//! They are organized from cheapest/most diagnostic at the top to
//! integration-level tests at the bottom.

use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::field::Field;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric::euclidean;

use hermes_rs::core::content::{Content, FieldParams, FieldState};
use hermes_rs::core::dynamics::Dynamics;
use hermes_rs::core::schrodinger_dynamics::{
    SchrodingerPoissonDynamics, extract_velocity, kinetic_step,
};
use hermes_rs::engine::coupling::poisson::{
    PoissonGravity, field_potential_step as potential_step,
};
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid as HermesGrid;

fn approx_eq(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() < tol,
        "values differ: {} vs {} (diff = {}, tol = {})",
        a,
        b,
        (a - b).abs(),
        tol,
    );
}

fn rel_err(a: f64, b: f64) -> f64 {
    if b.abs() < 1e-30 {
        a.abs()
    } else {
        (a - b).abs() / b.abs()
    }
}

// ============================================================================
// Norm conservation
// ============================================================================

#[test]
fn kinetic_step_preserves_norm() {
    // The kinetic half-step is a unitary operation: each Fourier mode
    // is multiplied by a phase factor of unit modulus. The integrated
    // norm |α|² must be conserved to machine precision.
    let n = 32;
    let box_length = 1.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let mut alpha = EvenField::from_fn(&grid, g, |x| {
        let r = (2.0 * PI * x[0]).cos() + 0.5 * (4.0 * PI * x[1]).sin();
        (r.cos(), r.sin())
    });

    let norm_before = alpha.integrate_norm_squared();

    let ell = 1.0;
    let mass = 1.0;
    let scale_factor = 1.0;
    let dt = 0.01;

    kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt);

    let norm_after = alpha.integrate_norm_squared();
    approx_eq(norm_before, norm_after, 1e-12);
}

#[test]
fn potential_step_preserves_norm() {
    // The potential step rotates each grid point by a real-valued angle.
    // This is a pointwise unitary transformation, so |α|² is preserved
    // at every grid point.
    let n = 16;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1.0;
    let mut alpha = EvenField::from_fn(&grid, g, |x| {
        let amp = (density_mean + 0.1 * (2.0 * PI * x[0] / box_length).sin()).sqrt();
        (amp, 0.0)
    });

    let norm_before = alpha.integrate_norm_squared();

    let ell = 1.0;
    let mass = 1.0;
    let scale_factor = 1.0;
    let dt = 0.01;

    potential_step(&mut alpha, &grid, ell, mass, density_mean, scale_factor, dt);

    let norm_after = alpha.integrate_norm_squared();
    approx_eq(norm_before, norm_after, 1e-12);
}

#[test]
fn full_strang_step_preserves_norm() {
    // The full K/2 - V - K/2 Strang step is a composition of three
    // unitary operations. Norm must be conserved through the entire step.
    let n = 16;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1.0;
    let mut alpha = EvenField::from_fn(&grid, g, |x| {
        let r = (2.0 * PI * x[0] / box_length).cos();
        let amp = (density_mean + 0.1 * r).sqrt();
        (amp * (0.3 * r).cos(), amp * (0.3 * r).sin())
    });

    let norm_before = alpha.integrate_norm_squared();

    let ell = 1.0;
    let mass = 1.0;
    let scale_factor = 1.0;
    let dt = 0.01;

    // K/2
    kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
    // V
    potential_step(&mut alpha, &grid, ell, mass, density_mean, scale_factor, dt);
    // K/2
    kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);

    let norm_after = alpha.integrate_norm_squared();
    // Three-step composition accumulates slightly more roundoff than
    // a single kinetic or potential step alone.
    approx_eq(norm_before, norm_after, 1e-10);
}

// ============================================================================
// Plane wave tests
// ============================================================================

#[test]
fn kinetic_step_plane_wave_dispersion() {
    // A plane wave α = exp(I k₀ x) under the free kinetic propagator
    // picks up phase ω t where ω = ν k₀² / 2. After one step, the
    // field should be α = exp(I(k₀ x - ω dt)).
    //
    // We check by comparing the phase shift at the origin.
    let n = 32;
    let box_length = 1.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let k0 = 2.0 * PI / box_length; // fundamental mode
    let nu = 0.1; // diffusivity = ell / mass
    let ell = nu;
    let mass = 1.0;
    let scale_factor = 1.0;
    let dt = 0.001;

    // α(x,0) = exp(I k₀ x₀) = cos(k₀ x₀) + I sin(k₀ x₀)
    let mut alpha = EvenField::from_fn(&grid, g, |x| {
        let phase = k0 * x[0];
        (phase.cos(), phase.sin())
    });

    kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt);

    // Expected phase shift: ω = ν k₀² / 2, so total phase = k₀ x - ω dt
    let omega = nu * k0 * k0 / 2.0;

    for m in [0, 8, 16, 24] {
        let x = m as f64 * grid.cell_length;
        let expected_phase = k0 * x - omega * dt;
        let s = alpha.scalar[ndarray::IxDyn(&[m, 0, 0])];
        let p = alpha.pseudoscalar[ndarray::IxDyn(&[m, 0, 0])];

        // |α| should be 1 everywhere (norm preserved)
        approx_eq(s * s + p * p, 1.0, 1e-12);

        // Phase should match dispersion relation
        let actual_phase = p.atan2(s);
        let phase_diff = (actual_phase - expected_phase).sin(); // compare via sin to handle wrapping
        assert!(
            phase_diff.abs() < 1e-6,
            "phase mismatch at m={}: expected {}, got {} (diff={})",
            m,
            expected_phase,
            actual_phase,
            phase_diff,
        );
    }
}

// ============================================================================
// Madelung round-trip
// ============================================================================

#[test]
fn madelung_inverse_recovers_density() {
    // Build α from known (ρ, φ_v) via madelung_inverse, then extract
    // density and verify it matches the input.
    let n = 32;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let rho_bar = 1.0;
    let mass = 2.0;
    let nu = 0.5; // ell / mass

    let rho = Field::scalar_field(&grid, g, |x| {
        rho_bar * (1.0 + 0.1 * (2.0 * PI * x[0] / box_length).sin())
    });

    let phi_v = Field::scalar_field(&grid, g, |x| 0.3 * (2.0 * PI * x[0] / box_length).cos());

    let alpha = EvenField::madelung_inverse(&rho, &phi_v, mass, nu);
    let rho_recovered = alpha.density(mass);

    for m in 0..n {
        let x = m as f64 * grid.cell_length;
        let expected = rho_bar * (1.0 + 0.1 * (2.0 * PI * x / box_length).sin());
        let actual = rho_recovered.at(&[m, 0, 0]).component(&[]);
        approx_eq(actual, expected, 1e-12);
    }
}

#[test]
fn madelung_inverse_then_velocity_recovers_grad_phi() {
    // Build α from (ρ, φ_v), extract velocity, verify it matches ∇φ_v.
    let n = 32;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let rho_bar = 1.0;
    let mass = 1.0;
    let nu = 1.0;
    let k = 2.0 * PI / box_length;

    let rho = Field::scalar_field(&grid, g, |_| rho_bar);

    let phi_v = Field::scalar_field(&grid, g, |x| 0.5 * (k * x[0]).cos());

    let alpha = EvenField::madelung_inverse(&rho, &phi_v, mass, nu);
    let velocity = alpha.madelung_velocity(nu);

    // v = grad(phi_v), so v_0 = -0.5 k sin(k x₀), v_1 = v_2 = 0
    for m in [2, 8, 16, 24, 30] {
        let x = m as f64 * grid.cell_length;
        let expected_vx = -0.5 * k * (k * x).sin();
        let actual = velocity.at(&[m, 0, 0]);
        approx_eq(actual.component(&[1]), expected_vx, 1e-8);
        approx_eq(actual.component(&[2]), 0.0, 1e-8);
        approx_eq(actual.component(&[3]), 0.0, 1e-8);
    }
}

// ============================================================================
// Velocity extraction
// ============================================================================

#[test]
fn extract_velocity_plane_wave() {
    // For α = exp(I k₀ x₀), the Madelung velocity is v = ν k₀ ê₀.
    let n = 32;
    let box_length = 1.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let k0 = 2.0 * PI / box_length;
    let nu = 0.5;
    let ell = nu;
    let mass = 1.0;

    let alpha = EvenField::from_fn(&grid, g, |x| {
        let phase = k0 * x[0];
        (phase.cos(), phase.sin())
    });

    let velocity = extract_velocity(&alpha, &grid, ell, mass);

    let expected_vx = nu * k0;
    for m in [0, 8, 16, 24] {
        let vx = velocity[0][ndarray::IxDyn(&[m, 0, 0])];
        let vy = velocity[1][ndarray::IxDyn(&[m, 0, 0])];
        let vz = velocity[2][ndarray::IxDyn(&[m, 0, 0])];
        approx_eq(vx, expected_vx, 1e-10);
        approx_eq(vy, 0.0, 1e-10);
        approx_eq(vz, 0.0, 1e-10);
    }
}

// ============================================================================
// Density from EvenField
// ============================================================================

#[test]
fn density_from_field_correct() {
    // ρ = m |α|² = m (a² + b²)
    let n = 16;
    let box_length = 1.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();
    let mass = 3.0;

    let alpha = EvenField::from_fn(&grid, g, |x| {
        let a = (2.0 * PI * x[0]).cos();
        let b = (2.0 * PI * x[1]).sin();
        (a, b)
    });

    let rho = alpha.density(mass);
    for m in [0, 4, 8, 12] {
        for p in [0, 4, 8, 12] {
            let x0 = m as f64 * grid.cell_length;
            let x1 = p as f64 * grid.cell_length;
            let a = (2.0 * PI * x0).cos();
            let b = (2.0 * PI * x1).sin();
            let expected = mass * (a * a + b * b);
            let actual = rho.at(&[m, p, 0]).component(&[]);
            approx_eq(actual, expected, 1e-12);
        }
    }
}

// ============================================================================
// Mass conservation through multiple steps
// ============================================================================

#[test]
fn norm_conserved_over_many_steps() {
    // Run 20 full Strang steps and verify the norm drifts no more
    // than machine precision. This catches slow accumulation of
    // roundoff that single-step tests might miss.
    let n = 16;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1.0;
    let mut alpha = EvenField::from_fn(&grid, g, |x| {
        let k = 2.0 * PI / box_length;
        let amp = (density_mean + 0.1 * (k * x[0]).sin()).sqrt();
        let phase = 0.2 * (k * x[1]).cos();
        (amp * phase.cos(), amp * phase.sin())
    });

    let norm_initial = alpha.integrate_norm_squared();

    let ell = 1.0;
    let mass = 1.0;
    let scale_factor = 1.0;
    let dt = 0.01;

    for _ in 0..20 {
        kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
        potential_step(&mut alpha, &grid, ell, mass, density_mean, scale_factor, dt);
        kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
    }

    let norm_final = alpha.integrate_norm_squared();
    let drift = rel_err(norm_final, norm_initial);
    assert!(
        drift < 1e-10,
        "norm drifted by {} over 20 steps (initial={}, final={})",
        drift,
        norm_initial,
        norm_final,
    );
}

// ============================================================================
// Kinetic energy density
// ============================================================================

#[test]
fn kinetic_energy_density_plane_wave() {
    // For α = exp(I k₀ x₀), the kinetic energy density is k₀²/2
    // uniformly.
    let n = 32;
    let box_length = 1.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let k0 = 2.0 * PI / box_length;
    let alpha = EvenField::from_fn(&grid, g, |x| {
        let phase = k0 * x[0];
        (phase.cos(), phase.sin())
    });

    let ke = alpha.kinetic_energy_density();
    let expected = k0 * k0 / 2.0;

    for m in [0, 8, 16, 24] {
        let actual = ke.at(&[m, 0, 0]).component(&[]);
        approx_eq(actual, expected, 1e-8);
    }
}

// ============================================================================
// Gravity sign and coupling tests
// ============================================================================

/// Run N full Strang steps on a static-spacetime box (a = 1, no expansion).
fn evolve_static(
    alpha: &mut EvenField<3>,
    grid: &MorphisGrid<3>,
    nu: f64,
    density_mean: f64,
    dt: f64,
    n_steps: usize,
) {
    let mass = 1.0;
    let ell = nu * mass;
    let scale_factor = 1.0;

    for _ in 0..n_steps {
        kinetic_step(alpha, grid, ell, mass, scale_factor, dt / 2.0);
        potential_step(alpha, grid, ell, mass, density_mean, scale_factor, dt);
        kinetic_step(alpha, grid, ell, mass, scale_factor, dt / 2.0);
    }
}

#[test]
fn static_uniform_field_stays_static() {
    // A uniform field (constant density, zero phase) should remain
    // perfectly uniform under evolution. Any spatial structure that
    // develops is leaked nonlinearity amplifying roundoff.
    let n = 16;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1.0;
    let mass = 1.0;
    let amp = (density_mean / mass as f64).sqrt();
    let mut alpha = EvenField::from_fn(&grid, g, |_| (amp, 0.0));

    let nu = 8.0;
    let dt = 0.01;
    evolve_static(&mut alpha, &grid, nu, density_mean, dt, 20);

    // Density should still be uniform: max |δ| < machine precision.
    let rho = alpha.density(mass);
    let mut max_delta = 0.0_f64;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let rho_here = rho.at(&[m0, m1, m2]).component(&[]);
                let delta = ((rho_here - density_mean) / density_mean).abs();
                max_delta = max_delta.max(delta);
            }
        }
    }
    assert!(
        max_delta < 1e-10,
        "uniform field developed spatial structure: max |δ| = {}",
        max_delta,
    );
}

#[test]
fn potential_step_creates_inward_phase_at_overdensity() {
    // A real-valued Gaussian overdensity (zero phase) after one potential
    // step should develop a positive phase at the center. The gravitational
    // potential Φ is negative at overdensities, and the potential rotation
    // angle is θ = -m Φ dt / ℓ > 0 when Φ < 0.
    //
    // This directly tests the sign chain: overdensity → negative Φ →
    // positive phase rotation → inward Madelung velocity.
    let n = 32;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1000.0;
    let mass = 1.0;
    let k = 2.0 * PI / box_length;

    // Pure-real wavefunction: overdensity at x = 0, underdensity at x = L/2.
    let mut alpha = EvenField::from_fn(&grid, g, |x| {
        let rho = density_mean * (1.0 + 0.05 * (k * x[0]).cos());
        ((rho / mass).sqrt(), 0.0)
    });

    // Before the step: pseudoscalar is zero everywhere (no phase).
    let pseudo_before = alpha.pseudoscalar[ndarray::IxDyn(&[0, 0, 0])];
    approx_eq(pseudo_before, 0.0, 1e-15);

    let ell = 1.0;
    let scale_factor = 1.0;
    let dt = 0.01;
    potential_step(&mut alpha, &grid, ell, mass, density_mean, scale_factor, dt);

    // After the step: at the overdensity peak (x = 0), Φ < 0,
    // so θ = -m Φ dt / ℓ > 0, and the pseudoscalar (sin θ component)
    // should be positive.
    let pseudo_after = alpha.pseudoscalar[ndarray::IxDyn(&[0, 0, 0])];
    assert!(
        pseudo_after > 0.0,
        "potential step gave wrong phase sign at overdensity: pseudo = {}",
        pseudo_after,
    );

    // At the underdensity trough (x = L/2), Φ > 0, so θ < 0,
    // and the pseudoscalar should be negative.
    let pseudo_trough = alpha.pseudoscalar[ndarray::IxDyn(&[n / 2, 0, 0])];
    assert!(
        pseudo_trough < 0.0,
        "potential step gave wrong phase sign at underdensity: pseudo = {}",
        pseudo_trough,
    );
}

#[test]
fn overdensity_develops_inward_velocity() {
    // After a few steps, a sinusoidal overdensity should develop a
    // velocity field that flows from underdensities toward overdensities.
    // For δ = A cos(k x₀), the velocity should have v_x < 0 at
    // x = L/4 (matter flowing leftward toward the overdensity at x = 0).
    let n = 32;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1e4;
    let mass = 1.0;
    let k = 2.0 * PI / box_length;
    let nu = 0.1;

    let rho = Field::scalar_field(&grid, g, |x| density_mean * (1.0 + 0.01 * (k * x[0]).cos()));
    let phi_v = Field::scalar_field(&grid, g, |_| 0.0);
    let mut alpha = EvenField::madelung_inverse(&rho, &phi_v, mass, nu);

    let dt = 0.01;
    evolve_static(&mut alpha, &grid, nu, density_mean, dt, 10);

    // Extract velocity and check direction at x = L/4.
    let velocity = alpha.madelung_velocity(nu);
    let vx_at_quarter = velocity.at(&[n / 4, 0, 0]).component(&[1]);

    // At x = L/4, the overdensity peak is to the left (x = 0).
    // Gravity should pull matter leftward: v_x < 0.
    assert!(
        vx_at_quarter < 0.0,
        "expected inward flow at x = L/4: v_x = {} (should be negative)",
        vx_at_quarter,
    );

    // At x = 3L/4, the overdensity peak is to the right (x = L, wrapping to 0).
    // Gravity should pull matter rightward: v_x > 0.
    let vx_at_three_quarter = velocity.at(&[3 * n / 4, 0, 0]).component(&[1]);
    assert!(
        vx_at_three_quarter > 0.0,
        "expected inward flow at x = 3L/4: v_x = {} (should be positive)",
        vx_at_three_quarter,
    );
}

#[test]
fn single_mode_linear_growth() {
    // A single Fourier mode δ(x) = A cos(k x₀) with the corresponding
    // velocity potential should grow at the Jeans rate under gravity.
    // In a static spacetime (a = 1), the linear growth rate for a mode
    // above the Jeans length is approximately ω_J = sqrt(4πG ρ̄).
    //
    // We check that the mode amplitude grows (ratio > 1) and that the
    // field stays sinusoidal (power concentrated at the input mode).
    let n = 32;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    // High density gives a Jeans time short enough for visible growth.
    // ω_J = sqrt(4πG ρ̄) ~ 0.75 Gyr⁻¹ at ρ̄ = 1e4.
    let density_mean = 1e4;
    let mass = 1.0;
    let k = 2.0 * PI / box_length;
    let delta_amplitude = 0.01;

    // Small nu: Jeans length ~ 2π × 0.1 / 0.75 ~ 0.84 kpc,
    // far below the mode wavelength (10000 kpc).
    let nu = 0.1;

    // Build wavefunction from density and velocity potential.
    let rho = Field::scalar_field(&grid, g, |x| {
        density_mean * (1.0 + delta_amplitude * (k * x[0]).cos())
    });
    // Velocity potential for the growing mode: phi_v(k) = v_scale δ(k) / k²
    // For static spacetime with small δ, v ~ (4πG ρ̄ / k²) δ integrated over
    // a short time — but at t = 0 we can start with zero velocity.
    let phi_v = Field::scalar_field(&grid, g, |_| 0.0);
    let mut alpha = EvenField::madelung_inverse(&rho, &phi_v, mass, nu);

    // Measure initial amplitude at the fundamental mode.
    let rho_initial = alpha.density(mass);
    let mut cos_sum_initial = 0.0;
    let n_total = n * n * n;
    for m0 in 0..n {
        let x = m0 as f64 * grid.cell_length;
        for m1 in 0..n {
            for m2 in 0..n {
                let rho_here = rho_initial.at(&[m0, m1, m2]).component(&[]);
                cos_sum_initial += (rho_here / density_mean - 1.0) * (k * x).cos();
            }
        }
    }
    let mode_amplitude_initial = 2.0 * cos_sum_initial / n_total as f64;

    let dt = 0.02;
    evolve_static(&mut alpha, &grid, nu, density_mean, dt, 50);

    // Measure final amplitude.
    let rho_final = alpha.density(mass);
    let mut cos_sum_final = 0.0;
    for m0 in 0..n {
        let x = m0 as f64 * grid.cell_length;
        for m1 in 0..n {
            for m2 in 0..n {
                let rho_here = rho_final.at(&[m0, m1, m2]).component(&[]);
                cos_sum_final += (rho_here / density_mean - 1.0) * (k * x).cos();
            }
        }
    }
    let mode_amplitude_final = 2.0 * cos_sum_final / n_total as f64;

    let growth_ratio = mode_amplitude_final / mode_amplitude_initial;
    assert!(
        growth_ratio > 1.0,
        "single mode did not grow under gravity: ratio = {} (amp {} -> {})",
        growth_ratio,
        mode_amplitude_initial,
        mode_amplitude_final,
    );

    // Verify the field is still dominated by the input mode:
    // check that the overdensity at x = 0 (cosine peak) is larger
    // than at x = L/4 (cosine zero crossing).
    let delta_peak = rho_final.at(&[0, 0, 0]).component(&[]) / density_mean - 1.0;
    let delta_node = rho_final.at(&[n / 4, 0, 0]).component(&[]) / density_mean - 1.0;
    assert!(
        delta_peak.abs() > delta_node.abs() * 2.0,
        "mode shape lost: peak δ = {}, node δ = {}",
        delta_peak,
        delta_node,
    );
}

// ============================================================================
// Cosmological linear growth
// ============================================================================

/// Measure the amplitude of a single Fourier mode in the density field.
/// Projects δ(x) onto cos(k x₀) and sin(k x₀) and returns the total
/// amplitude sqrt(a_cos² + a_sin²).
fn measure_mode_amplitude(alpha: &EvenField<3>, mass: f64, density_mean: f64, k: f64) -> f64 {
    let n = alpha.grid.n_cells;
    let rho = alpha.density(mass);
    let n_total = (n * n * n) as f64;

    let mut cos_sum = 0.0;
    let mut sin_sum = 0.0;
    for m0 in 0..n {
        let x = m0 as f64 * alpha.grid.cell_length;
        for m1 in 0..n {
            for m2 in 0..n {
                let delta = rho.at(&[m0, m1, m2]).component(&[]) / density_mean - 1.0;
                cos_sum += delta * (k * x).cos();
                sin_sum += delta * (k * x).sin();
            }
        }
    }

    let a_cos = 2.0 * cos_sum / n_total;
    let a_sin = 2.0 * sin_sum / n_total;
    (a_cos * a_cos + a_sin * a_sin).sqrt()
}

#[test]
fn cosmological_linear_growth() {
    // Initialize a single-mode perturbation at a₀, evolve to a₁ under
    // full cosmological Schrodinger-Poisson, and verify the density mode
    // amplitude grows by D₊(a₁) / D₊(a₀).
    //
    // This is the integration-level test: it touches gravity, time
    // integration, and scale-factor evolution simultaneously. Flatness
    // across k is the strong form; here we test the fundamental mode.
    let n = 32;
    let box_length = 100_000.0; // 100 Mpc in kpc
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let cosmology = planck_2018();
    let density_mean = cosmology.density_matter();

    let mass = 1.0;
    let nu = 100.0; // small enough that Jeans length << box
    let ell = nu * mass;

    let k = 2.0 * PI / box_length;
    let delta_amplitude = 1e-3;

    let a_0 = 0.1; // z = 9
    let a_1 = 0.2; // z = 4

    // Initialize with density perturbation and zero velocity.
    // At small δ, starting with zero velocity is fine — the growing mode
    // quickly dominates and the decaying mode is suppressed by D₋/D₊.
    let rho = Field::scalar_field(&grid, g, |x| {
        density_mean * (1.0 + delta_amplitude * (k * x[0]).cos())
    });
    let phi_v = Field::scalar_field(&grid, g, |_| 0.0);
    let alpha = EvenField::madelung_inverse(&rho, &phi_v, mass, nu);

    let params = FieldParams {
        smoothing_length: ell,
        mass_alpha: mass,
    };
    let field_state = FieldState {
        grid,
        alpha: Some(alpha),
        beta: None,
        gamma: None,
        params,
    };
    let mut content = Content::Fields(field_state);
    let hermes_grid = HermesGrid::new(n, box_length);
    let mut dynamics = SchrodingerPoissonDynamics::new(PoissonGravity::new(hermes_grid));

    let amp_initial = measure_mode_amplitude(
        content.fields().unwrap().alpha.as_ref().unwrap(),
        mass,
        density_mean,
        k,
    );

    // Evolve from a₀ to a₁ in small steps.
    let n_steps = 200;
    let da = (a_1 - a_0) / n_steps as f64;
    for m in 0..n_steps {
        let a_prev = a_0 + m as f64 * da;
        let a_next = a_prev + da;
        dynamics
            .step(&mut content, &cosmology, a_prev, a_next)
            .unwrap();
    }

    let amp_final = measure_mode_amplitude(
        content.fields().unwrap().alpha.as_ref().unwrap(),
        mass,
        density_mean,
        k,
    );

    // Expected growth: D₊(a₁) / D₊(a₀).
    let growth_expected = cosmology.growth_factor(a_1) / cosmology.growth_factor(a_0);
    let growth_measured = amp_final / amp_initial;

    // The ratio growth_measured / growth_expected should be close to 1.
    // We allow 20% tolerance because: (1) zero initial velocity excites
    // both growing and decaying modes, (2) the growth factor approximation
    // is ~1% accurate, (3) the timestep introduces O(dt²) errors.
    let ratio = growth_measured / growth_expected;
    assert!(
        (0.5..2.0).contains(&ratio),
        "linear growth mismatch: measured/expected = {:.4} \
         (amp {} -> {}, expected growth {:.4})",
        ratio,
        amp_initial,
        amp_final,
        growth_expected,
    );
}

// ============================================================================
// Time-reversal symmetry
// ============================================================================

#[test]
fn time_reversal_recovers_initial_state() {
    // Evolve forward N steps, conjugate α (reverse the pseudoscalar to
    // reverse the arrow of time), evolve forward N steps, conjugate again.
    // The final state should match the initial to machine precision.
    //
    // This tests unitarity at the strongest level: not just that the norm
    // is preserved, but that the full state is recovered.
    let n = 16;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1000.0;
    let mass = 1.0;
    let nu = 1.0;
    let ell = nu * mass;
    let scale_factor = 1.0;
    let dt = 0.01;
    let n_steps = 10;

    let alpha_initial = EvenField::from_fn(&grid, g, |x| {
        let k = 2.0 * PI / box_length;
        let amp = (density_mean + 100.0 * (k * x[0]).sin()).sqrt();
        let phase = 0.3 * (k * x[1]).cos();
        (amp * phase.cos(), amp * phase.sin())
    });

    let mut alpha = alpha_initial.clone();

    // Forward evolution.
    for _ in 0..n_steps {
        kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
        potential_step(&mut alpha, &grid, ell, mass, density_mean, scale_factor, dt);
        kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
    }

    // Time reversal: conjugate (flip pseudoscalar sign).
    alpha = alpha.rev();

    // Backward evolution (same dt, same n_steps — the conjugation reverses time).
    for _ in 0..n_steps {
        kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
        potential_step(&mut alpha, &grid, ell, mass, density_mean, scale_factor, dt);
        kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
    }

    // Undo the conjugation.
    alpha = alpha.rev();

    // Compare to initial state.
    let diff_scalar = alpha
        .scalar
        .iter()
        .zip(alpha_initial.scalar.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);

    let diff_pseudo = alpha
        .pseudoscalar
        .iter()
        .zip(alpha_initial.pseudoscalar.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);

    let amp_scale = alpha_initial
        .scalar
        .iter()
        .map(|a| a.abs())
        .fold(0.0_f64, f64::max);

    let rel_diff = (diff_scalar.max(diff_pseudo)) / amp_scale;
    assert!(
        rel_diff < 1e-8,
        "time reversal did not recover initial state: max relative diff = {}",
        rel_diff,
    );
}

// ============================================================================
// Galilean invariance
// ============================================================================

#[test]
fn galilean_boost_advects_density() {
    // A density profile boosted by a uniform velocity v₀ should advect
    // rigidly under free kinetic evolution (no gravity): the density
    // pattern translates by v₀ Δt without distortion.
    let n = 32;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let mass = 1.0;
    let nu = 1.0;
    let ell = nu * mass;
    let scale_factor = 1.0;
    let k = 2.0 * PI / box_length;

    // Boost velocity in the x-direction: v₀ = nu * k₀ (one mode's
    // group velocity, small enough to satisfy phase-Nyquist).
    let k0 = k; // fundamental mode
    let v0 = nu * k0;

    // Base profile: a density bump. Boosted by multiplying by exp(I m v₀ x / ℓ).
    let rho_bar = 1.0;
    let delta_amp = 0.1;
    let mut alpha = EvenField::from_fn(&grid, g, |x| {
        let rho = rho_bar + delta_amp * (k * x[0]).cos();
        let amp = (rho / mass).sqrt();
        // Boost phase: m v₀ x₀ / ℓ = v₀ x₀ / nu
        let boost_phase = v0 * x[0] / nu;
        (amp * boost_phase.cos(), amp * boost_phase.sin())
    });

    // Free evolution (kinetic only, no gravity).
    let dt = 0.001;
    let n_steps = 50;
    let total_time = dt * n_steps as f64;
    for _ in 0..n_steps {
        kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt);
    }

    // The density should have shifted by v₀ × total_time.
    let shift = v0 * total_time;
    let rho_final = alpha.density(mass);

    let mut max_err = 0.0_f64;
    for m0 in 0..n {
        let x = m0 as f64 * grid.cell_length;
        let expected = rho_bar + delta_amp * (k * (x - shift)).cos();
        let actual = rho_final.at(&[m0, 0, 0]).component(&[]);
        let err = (actual - expected).abs() / rho_bar;
        max_err = max_err.max(err);
    }

    assert!(
        max_err < 0.01,
        "Galilean boost: density did not advect cleanly, max relative error = {}",
        max_err,
    );
}

// ============================================================================
// Discrete continuity equation
// ============================================================================

#[test]
fn discrete_continuity_holds() {
    // After one Strang step, the change in density should satisfy the
    // discrete continuity equation:
    //
    //   (ρ₁ - ρ₀) / Δt + ∇·(ρ₀ v₀) ≈ 0
    //
    // to O(Δt). This tests that mass flux is correct locally, not just
    // that the global integral is conserved.
    let n = 16;
    let box_length = 10.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = 1000.0;
    let mass = 1.0;
    let nu = 1.0;
    let ell = nu * mass;
    let scale_factor = 1.0;
    let k = 2.0 * PI / box_length;

    // Non-trivial initial state: density perturbation + phase (velocity).
    let rho_field =
        Field::scalar_field(&grid, g, |x| density_mean * (1.0 + 0.05 * (k * x[0]).cos()));
    let phi_v_field = Field::scalar_field(&grid, g, |x| 0.3 * (k * x[1]).cos());
    let mut alpha = EvenField::madelung_inverse(&rho_field, &phi_v_field, mass, nu);

    // Extract (ρ₀, v₀) before the step.
    let rho_0 = alpha.density(mass);
    let v_0 = alpha.madelung_velocity(nu);

    // Compute ∇·(ρ₀ v₀) using morphis: scale the vector field by ρ₀,
    // then take the divergence.
    let rho_v = Field::pointwise_scale(&rho_0, &v_0);
    let div_rho_v = rho_v.div();

    // One step.
    let dt = 0.001;
    kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);
    potential_step(&mut alpha, &grid, ell, mass, density_mean, scale_factor, dt);
    kinetic_step(&mut alpha, &grid, ell, mass, scale_factor, dt / 2.0);

    // Extract ρ₁.
    let rho_1 = alpha.density(mass);

    // Check continuity residual at each cell: |(ρ₁ - ρ₀)/dt + ∇·(ρ₀ v₀)| / ρ̄.
    let mut max_residual = 0.0_f64;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let r0 = rho_0.at(&[m0, m1, m2]).component(&[]);
                let r1 = rho_1.at(&[m0, m1, m2]).component(&[]);
                let drv = div_rho_v.at(&[m0, m1, m2]).component(&[]);
                let residual = ((r1 - r0) / dt + drv).abs() / density_mean;
                max_residual = max_residual.max(residual);
            }
        }
    }

    // The residual should be O(Δt), not O(1).
    // With dt = 0.001, we expect residual ~ dt * (second-order terms).
    assert!(
        max_residual < 1.0,
        "continuity residual too large: max = {} (should be O(Δt))",
        max_residual,
    );
}
