use hermes_rs::config::build_configuration;
use hermes_rs::physics::simulation::Simulation;

/// Small configuration for fast tests: 8³ grid, 8³ particles, 5 steps.
fn small_config() -> hermes_rs::config::Configuration {
    let overrides: toml::Value = toml::from_str(
        r#"
        [simulation]
        n_cells = 8
        n_particles = 8
        box_length = 100000.0

        [time]
        scale_factor_initial = 0.02
        scale_factor_final = 0.05
        n_steps = 5
        stepping = "log_a"

        [output]
        snapshot_interval = 1
        "#,
    )
    .unwrap();

    build_configuration(None, Some(&overrides)).unwrap()
}

// ============================================================================
// Smoke test
// ============================================================================

#[test]
fn simulation_runs_without_panic() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    let n_steps = sim.run(&mut []).unwrap();
    assert_eq!(n_steps, 5);
}

#[test]
fn simulation_records_diagnostics() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    sim.run(&mut []).unwrap();

    // With snapshot_interval = 1 and 5 steps, we get initial + 5 snapshots = 6.
    assert_eq!(
        sim.diagnostics_history.len(),
        6,
        "should have initial + 5 step diagnostics"
    );
}

#[test]
fn simulation_reaches_final_scale_factor() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    sim.run(&mut []).unwrap();

    assert!(
        (sim.scale_factor - 0.05).abs() < 1e-6,
        "should reach final scale factor 0.05, got {}",
        sim.scale_factor
    );
}

// ============================================================================
// Conservation
// ============================================================================

#[test]
fn mass_conserved_throughout_run() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    sim.run(&mut []).unwrap();

    let mass_initial = sim.diagnostics_history[0].mass_total;
    for (n, diag) in sim.diagnostics_history.iter().enumerate() {
        let rel_err = (diag.mass_total - mass_initial).abs() / mass_initial;
        assert!(
            rel_err < 1e-12,
            "mass not conserved at step {n}: initial {mass_initial}, current {}",
            diag.mass_total
        );
    }
}

#[test]
fn momentum_stays_small() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    sim.run(&mut []).unwrap();

    // Total comoving momentum should remain small (not exactly zero due to
    // finite Fourier-space sampling, but bounded).
    for diag in &sim.diagnostics_history {
        let momentum_norm = diag.momentum_magnitude();
        // This is a weak bound — just checking it doesn't blow up.
        assert!(
            momentum_norm < 1e10,
            "total momentum blew up: {}",
            momentum_norm
        );
    }
}

// ============================================================================
// Diagnostics content
// ============================================================================

#[test]
fn diagnostics_have_correct_scale_factors() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    sim.run(&mut []).unwrap();

    // First diagnostics at initial scale factor.
    assert!(
        (sim.diagnostics_history[0].scale_factor - 0.02).abs() < 1e-6,
        "initial diagnostics should be at a = 0.02"
    );

    // Last diagnostics at final scale factor.
    let last = sim.diagnostics_history.last().unwrap();
    assert!(
        (last.scale_factor - 0.05).abs() < 1e-6,
        "final diagnostics should be at a = 0.05"
    );

    // Scale factors should be monotonically increasing.
    for n in 1..sim.diagnostics_history.len() {
        assert!(
            sim.diagnostics_history[n].scale_factor > sim.diagnostics_history[n - 1].scale_factor,
            "scale factors should increase monotonically"
        );
    }
}

#[test]
fn latest_diagnostics_is_final() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    sim.run(&mut []).unwrap();

    let latest = sim.latest_diagnostics().unwrap();
    assert!(
        (latest.scale_factor - 0.05).abs() < 1e-6,
        "latest diagnostics should be at final scale factor"
    );
}

#[test]
fn angular_momentum_grade_2_throughout() {
    let config = small_config();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    sim.run(&mut []).unwrap();

    for diag in &sim.diagnostics_history {
        assert_eq!(
            diag.angular_momentum.grade(),
            2,
            "angular momentum should always be grade-2 bivector"
        );
    }
}
