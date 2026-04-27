use hermes_rs::algebra::components_from_vector;
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::initial::{power_spectrum, transfer_function, zeldovich_init};

// ============================================================================
// Transfer function
// ============================================================================

#[test]
fn transfer_function_unity_at_low_k() {
    let cosmology = planck_2018();

    // T(k) → 1 for k → 0.
    let transfer = transfer_function(1e-4, &cosmology);
    assert!(
        (transfer - 1.0).abs() < 0.1,
        "T(k→0) should be ≈1, got {transfer}"
    );
}

#[test]
fn transfer_function_decreases_at_high_k() {
    let cosmology = planck_2018();

    let transfer_low = transfer_function(0.01, &cosmology);
    let transfer_high = transfer_function(1.0, &cosmology);

    assert!(
        transfer_high < transfer_low,
        "T(k) should decrease with k: T(0.01) = {transfer_low}, T(1.0) = {transfer_high}"
    );
}

#[test]
fn transfer_function_positive() {
    let cosmology = planck_2018();

    for &k in &[1e-4, 1e-3, 0.01, 0.1, 1.0, 10.0] {
        let transfer = transfer_function(k, &cosmology);
        assert!(
            transfer > 0.0,
            "T(k={k}) should be positive, got {transfer}"
        );
    }
}

// ============================================================================
// Power spectrum
// ============================================================================

#[test]
fn power_spectrum_positive() {
    let cosmology = planck_2018();

    for &k in &[0.001, 0.01, 0.1, 1.0] {
        let power = power_spectrum(k, &cosmology);
        assert!(power > 0.0, "P(k={k}) should be positive, got {power}");
    }
}

#[test]
fn power_spectrum_peaks_near_keq() {
    let cosmology = planck_2018();

    // The power spectrum P(k) = k^{n_s} T(k)² should peak near the
    // matter-radiation equality scale k_eq ~ 0.01 h/Mpc.
    let power_low = power_spectrum(0.001, &cosmology);
    let power_peak = power_spectrum(0.02, &cosmology);
    let power_high = power_spectrum(1.0, &cosmology);

    assert!(
        power_peak > power_low && power_peak > power_high,
        "P(k) should peak near k_eq: P(0.001) = {power_low}, P(0.02) = {power_peak}, P(1.0) = {power_high}"
    );
}

// ============================================================================
// Zel'dovich initialization
// ============================================================================

#[test]
fn zeldovich_deterministic() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();
    let scale_factor_initial = 0.02;

    let particles_1 = zeldovich_init(16, &grid, &cosmology, scale_factor_initial, 42).unwrap();
    let particles_2 = zeldovich_init(16, &grid, &cosmology, scale_factor_initial, 42).unwrap();

    for p in 0..particles_1.count() {
        let pos_1 = components_from_vector(&particles_1.position_of(p));
        let pos_2 = components_from_vector(&particles_2.position_of(p));
        for d in 0..3 {
            assert!(
                (pos_1[d] - pos_2[d]).abs() < 1e-15,
                "same seed should produce identical particles"
            );
        }
    }
}

#[test]
fn zeldovich_different_seeds_differ() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();

    let particles_1 = zeldovich_init(16, &grid, &cosmology, 0.02, 42).unwrap();
    let particles_2 = zeldovich_init(16, &grid, &cosmology, 0.02, 99).unwrap();

    let mut any_different = false;
    for p in 0..particles_1.count() {
        let pos_1 = components_from_vector(&particles_1.position_of(p));
        let pos_2 = components_from_vector(&particles_2.position_of(p));
        for d in 0..3 {
            if (pos_1[d] - pos_2[d]).abs() > 1e-10 {
                any_different = true;
            }
        }
    }
    assert!(
        any_different,
        "different seeds should produce different ICs"
    );
}

#[test]
fn zeldovich_mass_conservation() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();

    let particles = zeldovich_init(16, &grid, &cosmology, 0.02, 42).unwrap();

    let expected = cosmology.density_matter() * grid.box_volume();
    let computed = particles.total_mass();
    let rel_err = (computed - expected).abs() / expected;

    assert!(
        rel_err < 1e-12,
        "mass conservation: expected {expected}, got {computed}, rel_err {rel_err}"
    );
}

#[test]
fn zeldovich_positions_in_box() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();

    let particles = zeldovich_init(16, &grid, &cosmology, 0.02, 42).unwrap();

    for p in 0..particles.count() {
        let pos = components_from_vector(&particles.position_of(p));
        for d in 0..3 {
            assert!(
                pos[d] >= 0.0 && pos[d] < grid.box_length,
                "particle {p} component {d} = {} out of bounds",
                pos[d]
            );
        }
    }
}

#[test]
fn zeldovich_momenta_nonzero() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();

    let particles = zeldovich_init(16, &grid, &cosmology, 0.02, 42).unwrap();

    // Zel'dovich ICs should give nonzero momenta (velocity ∝ displacement).
    let total_momentum_norm = particles.total_momentum().norm();
    let any_nonzero = (0..particles.count()).any(|p| particles.momentum_of(p).norm() > 1e-20);

    assert!(any_nonzero, "Zel'dovich ICs should have nonzero momenta");

    // Total momentum should be approximately zero by symmetry of the random field,
    // but not exactly (finite realization).
    let typical_momentum = particles.momentum_of(0).norm();
    let n_particles = particles.count() as f64;
    let expected_rms = typical_momentum * n_particles.sqrt();

    // Total momentum should be much smaller than the sum of magnitudes.
    assert!(
        total_momentum_norm < expected_rms,
        "total momentum {total_momentum_norm} should be less than RMS {expected_rms}"
    );
}

#[test]
fn zeldovich_positions_are_grade_1_vectors() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();

    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();

    let pos = particles.position_of(0);
    let mom = particles.momentum_of(0);

    assert_eq!(pos.grade(), 1, "position should be grade-1 vector");
    assert_eq!(mom.grade(), 1, "momentum should be grade-1 vector");
}

#[test]
fn zeldovich_angular_momentum_is_bivector() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();

    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();

    // Angular momentum L = x ∧ p should be a grade-2 bivector.
    let angular_momentum = particles.angular_momentum(0);
    assert_eq!(
        angular_momentum.grade(),
        2,
        "angular momentum should be grade-2 bivector"
    );
}
