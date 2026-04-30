use hermes_rs::config::{build_configuration, load_defaults};

// ============================================================================
// Default loading
// ============================================================================

#[test]
fn load_defaults_succeeds() {
    let config = load_defaults().expect("defaults should load");

    assert!((config.cosmology.hubble - 0.674).abs() < 1e-10);
    assert!((config.cosmology.omega_m - 0.315).abs() < 1e-10);
    assert_eq!(config.simulation.n_grid, 64);
    assert_eq!(config.simulation.n_particles, 32);
    assert_eq!(config.time.n_steps, 300);
    assert_eq!(config.output.write_interval, 1);
    assert_eq!(config.output.diagnostic_interval, 10);
}

#[test]
fn defaults_cosmology_validates() {
    let config = load_defaults().expect("defaults should load");
    assert!(config.cosmology.validate().is_ok());
}

// ============================================================================
// Deep merge
// ============================================================================

#[test]
fn partial_override_merges_correctly() {
    let override_val: toml::Value = toml::from_str(
        r#"
        [cosmology]
        omega_m = 0.30
        omega_lambda = 0.6999085

        [simulation]
        n_grid = 128
        "#,
    )
    .unwrap();

    let config = build_configuration(None, Some(&override_val)).expect("merge should succeed");

    // Overridden values
    assert!((config.cosmology.omega_m - 0.30).abs() < 1e-10);
    assert_eq!(config.simulation.n_grid, 128);

    // Non-overridden values remain at defaults
    assert!((config.cosmology.hubble - 0.674).abs() < 1e-10);
    assert_eq!(config.simulation.n_particles, 32);
    assert_eq!(config.time.n_steps, 300);
}

#[test]
fn override_rejects_invalid_cosmology() {
    let override_val: toml::Value = toml::from_str(
        r#"
        [cosmology]
        omega_m = -0.1
        "#,
    )
    .unwrap();

    let result = build_configuration(None, Some(&override_val));
    assert!(result.is_err(), "negative omega_m should be rejected");
}

#[test]
fn two_tier_override() {
    let file_val: toml::Value = toml::from_str(
        r#"
        [simulation]
        n_grid = 128
        n_particles = 128
        "#,
    )
    .unwrap();

    let override_val: toml::Value = toml::from_str(
        r#"
        [simulation]
        n_particles = 256
        "#,
    )
    .unwrap();

    let config =
        build_configuration(Some(&file_val), Some(&override_val)).expect("merge should succeed");

    // file_val sets n_grid = 128, override_val overrides n_particles to 256
    assert_eq!(config.simulation.n_grid, 128);
    assert_eq!(config.simulation.n_particles, 256);
}
