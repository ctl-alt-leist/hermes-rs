use hermes_rs::config::{FieldGrade, SpacetimeBackground, load_scene_config};
use std::path::Path;

// ============================================================================
// Scene config parsing
// ============================================================================

#[test]
fn cosmic_web_pm_parses() {
    let config = load_scene_config(Path::new("scenes/cosmic-web-pm.toml")).unwrap();

    assert_eq!(
        config.ontology.spacetime.background,
        SpacetimeBackground::Flrw
    );
    assert!(config.ontology.has_particles());
    assert!(!config.ontology.has_fields());
    assert!(config.ontology.has_gravity());

    let dm = &config.ontology.particles["dark_matter"];
    assert_eq!(dm.count, 262144);
    assert_eq!(dm.mass, 1e10);

    assert_eq!(config.simulation.grid.n_cells, 64);
    assert_eq!(config.simulation.grid.box_length, 100000.0);
    assert_eq!(config.simulation.time.n_steps, 300);
    assert_eq!(config.simulation.time.stepping, "log");
}

#[test]
fn cosmic_web_field_parses() {
    let config = load_scene_config(Path::new("scenes/cosmic-web-ft.toml")).unwrap();

    assert_eq!(
        config.ontology.spacetime.background,
        SpacetimeBackground::Flrw
    );
    assert!(!config.ontology.has_particles());
    assert!(config.ontology.has_fields());
    assert!(config.ontology.has_gravity());

    let alpha = &config.ontology.fields["alpha"];
    assert_eq!(alpha.grade, FieldGrade::Multi(vec![0, 3]));
    assert_eq!(alpha.mass, Some(1e10));
    assert_eq!(alpha.length_scale, Some(2000.0));
    assert_eq!(alpha.free.as_deref(), Some("schrodinger"));

    assert_eq!(config.simulation.grid.box_length, 10000.0);
}

#[test]
fn galaxy_group_pm_parses() {
    let config = load_scene_config(Path::new("scenes/galaxy-group-pm.toml")).unwrap();

    let dm = &config.ontology.particles["dark_matter"];
    assert_eq!(dm.count, 32768);
    assert_eq!(config.simulation.time.stepping, "linear");
    assert_eq!(config.simulation.time.n_steps, 600);
}

#[test]
fn static_wave_packet_parses() {
    let config = load_scene_config(Path::new("scenes/static-wave-packet.toml")).unwrap();

    assert_eq!(
        config.ontology.spacetime.background,
        SpacetimeBackground::Static
    );
    assert!(!config.ontology.has_gravity());
    assert!(config.ontology.has_fields());

    let alpha = &config.ontology.fields["alpha"];
    assert_eq!(alpha.free.as_deref(), Some("schrodinger"));

    assert!(config.simulation.time.time_range.is_some());
    assert_eq!(config.simulation.time.time_range.unwrap(), [0.0, 10.0]);
}

#[test]
fn multi_species_parses() {
    let config = load_scene_config(Path::new("scenes/multi-species.toml")).unwrap();

    assert!(config.ontology.has_particles());
    assert!(config.ontology.has_fields());
    assert!(config.ontology.has_gravity());
    assert_eq!(config.ontology.lagrangian.electromagnetic.len(), 2);

    assert!(config.ontology.fields.contains_key("alpha"));
    assert!(config.ontology.fields.contains_key("beta"));
    assert!(config.ontology.fields.contains_key("gamma"));

    let gamma = &config.ontology.fields["gamma"];
    assert_eq!(gamma.grade, FieldGrade::Single(2));
    assert_eq!(gamma.free.as_deref(), Some("maxwell"));
}

#[test]
fn propagating_gravity_parses() {
    let config = load_scene_config(Path::new("scenes/propagating-gravity.toml")).unwrap();

    assert!(config.ontology.fields.contains_key("phi"));
    let phi = &config.ontology.fields["phi"];
    assert_eq!(phi.grade, FieldGrade::Single(0));
    assert_eq!(phi.free.as_deref(), Some("wave"));
    assert_eq!(phi.speed, Some(3e5));
}

// ============================================================================
// Validation
// ============================================================================

#[test]
fn flrw_requires_cosmological_parameters() {
    let toml_str = r#"
        [ontology.spacetime]
        background = "flrw"

        [ontology.particles.test]
        count = 8
        mass = 1.0

        [ontology.lagrangian]
        gravity = true

        [simulation.grid]
        n_cells = 8
        box_length = 1.0

        [simulation.time]
        scale_factor_range = [0.1, 1.0]
        n_steps = 10
    "#;

    let config: hermes_rs::config::EngineConfig = toml::from_str(toml_str).unwrap();
    let result = config.validate();
    assert!(result.is_err());
}

#[test]
fn static_requires_time_range() {
    let toml_str = r#"
        [ontology.spacetime]
        background = "static"

        [ontology.fields.alpha]
        grade = [0, 3]
        mass = 1.0
        length_scale = 1.0
        free = "schrodinger"

        [simulation.grid]
        n_cells = 8
        box_length = 1.0

        [simulation.time]
        scale_factor_range = [0.1, 1.0]
        n_steps = 10
    "#;

    let config: hermes_rs::config::EngineConfig = toml::from_str(toml_str).unwrap();
    let result = config.validate();
    assert!(result.is_err());
}

#[test]
fn schrodinger_requires_length_scale() {
    let toml_str = r#"
        [ontology.spacetime]
        background = "static"

        [ontology.fields.alpha]
        grade = [0, 3]
        mass = 1.0
        free = "schrodinger"

        [simulation.grid]
        n_cells = 8
        box_length = 1.0

        [simulation.time]
        time_range = [0.0, 1.0]
        n_steps = 10
    "#;

    let config: hermes_rs::config::EngineConfig = toml::from_str(toml_str).unwrap();
    let result = config.validate();
    assert!(result.is_err());
}

#[test]
fn empty_ontology_fails_validation() {
    let toml_str = r#"
        [ontology.spacetime]
        background = "static"

        [simulation.grid]
        n_cells = 8
        box_length = 1.0

        [simulation.time]
        time_range = [0.0, 1.0]
        n_steps = 10
    "#;

    let config: hermes_rs::config::EngineConfig = toml::from_str(toml_str).unwrap();
    let result = config.validate();
    assert!(result.is_err());
}

// ============================================================================
// Deep merge
// ============================================================================

#[test]
fn scene_overrides_base_defaults() {
    let config = load_scene_config(Path::new("scenes/cosmic-web-pm.toml")).unwrap();

    // Scene overrides base n_cells (32 → 64)
    assert_eq!(config.simulation.grid.n_cells, 64);

    // Base default for snapshot interval persists
    assert_eq!(config.output.snapshots.interval, 1);
}

#[test]
fn output_display_defaults_populated() {
    let config = load_scene_config(Path::new("scenes/cosmic-web-pm.toml")).unwrap();

    // Display defaults from base.toml should be present
    assert_eq!(config.output.display.point_size, 5.0);
    assert_eq!(config.output.display.camera_distance, 1.9);
    assert_eq!(config.output.display.gif_resolution, 512);
}
