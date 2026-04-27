use hermes_rs::constants::{CRITICAL_DENSITY_FACTOR, H100};
use hermes_rs::cosmology::{Cosmology, einstein_de_sitter, planck_2018};

// ============================================================================
// Einstein-de Sitter analytic solutions
//
// EdS: Ω_m = 1, Ω_Λ = 0, flat. Every quantity has a closed form:
//   E(a) = a^{-3/2}
//   t(a) = (2/3) a^{3/2} / H₀
//   D₊(a) = a
//   f(a) = 1
//   kick_factor(a0, a1)  = 2(√a1 - √a0) / H₀
//   drift_factor(a0, a1) = 2(1/√a0 - 1/√a1) / H₀
// ============================================================================

#[test]
fn eds_expansion_rate() {
    let eds = einstein_de_sitter();
    for &a in &[0.1_f64, 0.25, 0.5, 0.75, 1.0] {
        let expected = a.powf(-1.5);
        let computed = eds.expansion_rate(a);
        let rel_err = (computed - expected).abs() / expected;
        assert!(
            rel_err < 1e-12,
            "E({a}): expected {expected}, got {computed}, rel_err {rel_err}"
        );
    }
}

#[test]
fn eds_cosmic_time() {
    let eds = einstein_de_sitter();
    let h0 = eds.hubble_constant();

    for &a in &[0.1_f64, 0.25, 0.5, 0.75, 1.0] {
        let expected = (2.0 / 3.0) * a.powf(1.5) / h0;
        let computed = eds.cosmic_time(a);
        let rel_err = (computed - expected).abs() / expected;
        assert!(
            rel_err < 1e-4,
            "t({a}): expected {expected}, got {computed}, rel_err {rel_err}"
        );
    }
}

#[test]
fn eds_growth_factor() {
    let eds = einstein_de_sitter();

    for &a in &[0.1, 0.25, 0.5, 0.75, 1.0] {
        let computed = eds.growth_factor(a);
        let rel_err = (computed - a).abs() / a;
        assert!(
            rel_err < 1e-3,
            "D₊({a}): expected {a}, got {computed}, rel_err {rel_err}"
        );
    }
}

#[test]
fn eds_growth_rate() {
    let eds = einstein_de_sitter();

    for &a in &[0.1, 0.25, 0.5, 0.75, 1.0] {
        let computed = eds.growth_rate(a);
        let err = (computed - 1.0).abs();
        assert!(
            err < 1e-3,
            "f({a}): expected 1.0, got {computed}, err {err}"
        );
    }
}

#[test]
fn eds_kick_factor() {
    let eds = einstein_de_sitter();
    let h0 = eds.hubble_constant();

    let a0: f64 = 0.1;
    let a1: f64 = 0.5;
    let expected = 2.0 * (a1.sqrt() - a0.sqrt()) / h0;
    let computed = eds.kick_factor(a0, a1);
    let rel_err = (computed - expected).abs() / expected;
    assert!(
        rel_err < 1e-4,
        "kick_factor({a0}, {a1}): expected {expected}, got {computed}, rel_err {rel_err}"
    );
}

#[test]
fn eds_drift_factor() {
    let eds = einstein_de_sitter();
    let h0 = eds.hubble_constant();

    let a0: f64 = 0.1;
    let a1: f64 = 0.5;
    let expected = 2.0 * (1.0 / a0.sqrt() - 1.0 / a1.sqrt()) / h0;
    let computed = eds.drift_factor(a0, a1);
    let rel_err = (computed - expected).abs() / expected;
    assert!(
        rel_err < 1e-4,
        "drift_factor({a0}, {a1}): expected {expected}, got {computed}, rel_err {rel_err}"
    );
}

// ============================================================================
// Planck 2018 reference values
// ============================================================================

#[test]
fn planck_age_of_universe() {
    let cosmo = planck_2018();
    let age = cosmo.cosmic_time(1.0);

    assert!(
        (age - 13.8).abs() < 0.2,
        "age of universe: expected ≈13.8 Gyr, got {age}"
    );
}

#[test]
fn planck_hubble_constant() {
    let cosmo = planck_2018();
    let h0 = cosmo.hubble_constant();
    let expected = 0.674 * H100;
    let rel_err = (h0 - expected).abs() / expected;

    assert!(rel_err < 1e-10, "H₀: expected {expected}, got {h0}");
}

#[test]
fn planck_critical_density() {
    let cosmo = planck_2018();
    let rho_c = cosmo.density_critical(1.0);
    let expected = CRITICAL_DENSITY_FACTOR * cosmo.hubble_constant().powi(2);
    let rel_err = (rho_c - expected).abs() / expected;

    assert!(
        rel_err < 1e-10,
        "ρ_c(a=1): expected {expected}, got {rho_c}"
    );
}

#[test]
fn planck_growth_factor_normalization() {
    let cosmo = planck_2018();
    let d_1 = cosmo.growth_factor(1.0);

    assert!(
        (d_1 - 1.0).abs() < 1e-10,
        "D₊(1) should be 1.0 by normalization, got {d_1}"
    );
}

#[test]
fn planck_growth_factor_monotone() {
    let cosmo = planck_2018();
    let scale_factors = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0];
    let growth: Vec<f64> = scale_factors
        .iter()
        .map(|&a| cosmo.growth_factor(a))
        .collect();

    for n in 1..growth.len() {
        assert!(
            growth[n] > growth[n - 1],
            "D₊ must be monotonically increasing: D₊({}) = {} ≤ D₊({}) = {}",
            scale_factors[n],
            growth[n],
            scale_factors[n - 1],
            growth[n - 1]
        );
    }
}

#[test]
fn planck_growth_rate_reasonable() {
    let cosmo = planck_2018();

    // At a=1 for Planck 2018, f ≈ 0.52 (matter-dominated would give 1.0,
    // but Λ suppresses growth).
    let f_1 = cosmo.growth_rate(1.0);
    assert!(
        f_1 > 0.4 && f_1 < 0.7,
        "f(a=1) should be ≈0.52 for Planck 2018, got {f_1}"
    );

    // At high redshift (matter-dominated era), f → 1.
    let f_early = cosmo.growth_rate(0.01);
    assert!(
        (f_early - 1.0).abs() < 0.05,
        "f(a=0.01) should be ≈1.0 in the matter era, got {f_early}"
    );
}

// ============================================================================
// Derived quantities
// ============================================================================

#[test]
fn derived_omega_cdm() {
    let cosmo = planck_2018();
    let expected = cosmo.omega_m - cosmo.omega_b;
    let computed = cosmo.omega_cdm();

    assert!(
        (computed - expected).abs() < 1e-15,
        "omega_cdm: expected {expected}, got {computed}"
    );
}

#[test]
fn derived_baryon_fraction() {
    let cosmo = planck_2018();
    let expected = cosmo.omega_b / cosmo.omega_m;
    let computed = cosmo.baryon_fraction();

    assert!(
        (computed - expected).abs() < 1e-15,
        "baryon_fraction: expected {expected}, got {computed}"
    );
}

#[test]
fn derived_density_matter_comoving() {
    let cosmo = planck_2018();
    let rho_m = cosmo.density_matter();
    let expected = cosmo.omega_m * cosmo.density_critical(1.0);

    assert!(
        (rho_m - expected).abs() / expected < 1e-10,
        "density_matter: expected {expected}, got {rho_m}"
    );
}

#[test]
fn planck_is_flat() {
    let cosmo = planck_2018();
    assert!(cosmo.is_flat(), "Planck 2018 cosmology should be flat");
}

// ============================================================================
// Validation
// ============================================================================

#[test]
fn validation_passes_for_planck() {
    let cosmo = planck_2018();
    assert!(cosmo.validate().is_ok());
}

#[test]
fn validation_passes_for_eds() {
    let eds = einstein_de_sitter();
    assert!(eds.validate().is_ok());
}

#[test]
fn validation_rejects_negative_omega_m() {
    let cosmo = Cosmology {
        omega_m: -0.1,
        ..planck_2018()
    };
    assert!(cosmo.validate().is_err());
}

#[test]
fn validation_rejects_omega_b_exceeding_omega_m() {
    let cosmo = Cosmology {
        omega_b: 0.5,
        omega_m: 0.3,
        ..planck_2018()
    };
    assert!(cosmo.validate().is_err());
}

#[test]
fn validation_rejects_nonunit_density_sum() {
    let cosmo = Cosmology {
        omega_lambda: 0.0,
        ..planck_2018()
    };
    assert!(cosmo.validate().is_err());
}
