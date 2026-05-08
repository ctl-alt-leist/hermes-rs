/// Tests for the EngineConfig system.
///
/// These tests verify config parsing, validation, and deep merge
/// using inline TOML strings. No test depends on specific scene
/// file contents — scenes are user-tunable and not relied upon.
use hermes_rs::config::{EngineConfig, FieldGrade, SpacetimeBackground, load_scene_config};
use std::path::Path;

// ============================================================================
// Config parsing
// ============================================================================

#[test]
fn flrw_particle_config_parses() {
    let config: EngineConfig = toml::from_str(
        r#"
        [ontology.spacetime]
        background = "flrw"
        hubble = 67.4
        omega_m = 0.315
        omega_v = 0.685

        [ontology.particles.dark_matter]
        count = 8000
        mass = 1e10

        [ontology.lagrangian]
        gravity = true

        [simulation.grid]
        n_cells = 32
        box_length = 100000.0

        [simulation.time]
        scale_factor_range = [0.02, 1.0]
        n_steps = 100
    "#,
    )
    .unwrap();

    assert_eq!(
        config.ontology.spacetime.background,
        SpacetimeBackground::Flrw
    );
    assert!(config.ontology.has_particles());
    assert!(!config.ontology.has_fields());
    assert!(config.ontology.has_gravity());

    let dm = &config.ontology.particles["dark_matter"];
    assert_eq!(dm.count, 8000);
    assert_eq!(dm.mass, 1e10);
}

#[test]
fn flrw_field_config_parses() {
    let config: EngineConfig = toml::from_str(
        r#"
        [ontology.spacetime]
        background = "flrw"
        hubble = 67.4
        omega_m = 0.315
        omega_v = 0.685

        [ontology.fields.alpha]
        grade = [0, 3]
        mass = 1e10
        length_scale = 2000.0
        free = "schrodinger"

        [ontology.lagrangian]
        gravity = true

        [simulation.grid]
        n_cells = 64
        box_length = 10000.0

        [simulation.time]
        scale_factor_range = [0.1, 1.0]
        n_steps = 200
    "#,
    )
    .unwrap();

    assert!(!config.ontology.has_particles());
    assert!(config.ontology.has_fields());
    assert!(config.ontology.has_gravity());

    let alpha = &config.ontology.fields["alpha"];
    assert_eq!(alpha.grade, FieldGrade::Multi(vec![0, 3]));
    assert_eq!(alpha.mass, Some(1e10));
    assert_eq!(alpha.length_scale, Some(2000.0));
    assert_eq!(alpha.free.as_deref(), Some("schrodinger"));
}

#[test]
fn static_field_config_parses() {
    let config: EngineConfig = toml::from_str(
        r#"
        [ontology.spacetime]
        background = "static"

        [ontology.fields.alpha]
        grade = [0, 3]
        mass = 1.0
        length_scale = 1.0
        free = "schrodinger"

        [simulation.grid]
        n_cells = 16
        box_length = 1.0

        [simulation.time]
        time_range = [0.0, 10.0]
        n_steps = 100
    "#,
    )
    .unwrap();

    assert_eq!(
        config.ontology.spacetime.background,
        SpacetimeBackground::Static
    );
    assert!(!config.ontology.has_gravity());
    assert!(config.ontology.has_fields());
    assert!(config.simulation.time.time_range.is_some());
    assert_eq!(config.simulation.time.time_range.unwrap(), [0.0, 10.0]);
}

#[test]
fn multi_species_config_parses() {
    let config: EngineConfig = toml::from_str(
        r#"
        [ontology.spacetime]
        background = "flrw"
        hubble = 67.4
        omega_m = 0.315
        omega_v = 0.685

        [ontology.particles.dark_matter]
        count = 1000
        mass = 1e10

        [ontology.fields.alpha]
        grade = [0, 3]
        mass = 1e10
        length_scale = 2000.0
        free = "schrodinger"

        [ontology.fields.beta]
        grade = [0, 3]
        mass = 1e10
        length_scale = 2000.0
        free = "schrodinger"
        self_interaction = 1e6

        [ontology.fields.gamma]
        grade = 2
        free = "maxwell"
        speed = 3e5

        [ontology.lagrangian]
        gravity = true
        electromagnetic = ["beta", "gamma"]

        [simulation.grid]
        n_cells = 32
        box_length = 10000.0

        [simulation.time]
        scale_factor_range = [0.1, 1.0]
        n_steps = 100
    "#,
    )
    .unwrap();

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

    let beta = &config.ontology.fields["beta"];
    assert_eq!(beta.self_interaction, Some(1e6));
}

#[test]
fn propagating_gravity_config_parses() {
    let config: EngineConfig = toml::from_str(
        r#"
        [ontology.spacetime]
        background = "flrw"
        hubble = 67.4
        omega_m = 0.315
        omega_v = 0.685

        [ontology.fields.alpha]
        grade = [0, 3]
        mass = 1e10
        length_scale = 2000.0
        free = "schrodinger"

        [ontology.fields.phi]
        grade = 0
        free = "wave"
        speed = 3e5

        [ontology.lagrangian]
        gravity = true

        [simulation.grid]
        n_cells = 32
        box_length = 10000.0

        [simulation.time]
        scale_factor_range = [0.1, 1.0]
        n_steps = 100
    "#,
    )
    .unwrap();

    let phi = &config.ontology.fields["phi"];
    assert_eq!(phi.grade, FieldGrade::Single(0));
    assert_eq!(phi.free.as_deref(), Some("wave"));
    assert_eq!(phi.speed, Some(3e5));
}

#[test]
fn per_species_display_config_parses() {
    let config: EngineConfig = toml::from_str(
        r#"
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
        time_range = [0.0, 1.0]
        n_steps = 10

        [output.display]
        blob_size = 28.0

        [output.display.species.alpha]
        colormap = "cool"
        colormap_range = [0.1, 5.0]
    "#,
    )
    .unwrap();

    assert_eq!(config.output.display.blob_size, 28.0);

    let alpha_vis = &config.output.display.species["alpha"];
    assert_eq!(alpha_vis.colormap, "cool");
    assert_eq!(alpha_vis.colormap_range, Some([0.1, 5.0]));
}

// ============================================================================
// Validation
// ============================================================================

#[test]
fn flrw_requires_cosmological_parameters() {
    let config: EngineConfig = toml::from_str(
        r#"
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
    "#,
    )
    .unwrap();

    assert!(config.validate().is_err());
}

#[test]
fn static_requires_time_range() {
    let config: EngineConfig = toml::from_str(
        r#"
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
    "#,
    )
    .unwrap();

    assert!(config.validate().is_err());
}

#[test]
fn schrodinger_requires_length_scale() {
    let config: EngineConfig = toml::from_str(
        r#"
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
    "#,
    )
    .unwrap();

    assert!(config.validate().is_err());
}

#[test]
fn empty_ontology_fails_validation() {
    let config: EngineConfig = toml::from_str(
        r#"
        [ontology.spacetime]
        background = "static"

        [simulation.grid]
        n_cells = 8
        box_length = 1.0

        [simulation.time]
        time_range = [0.0, 1.0]
        n_steps = 10
    "#,
    )
    .unwrap();

    assert!(config.validate().is_err());
}

// ============================================================================
// Deep merge
// ============================================================================

#[test]
fn scene_overrides_base_defaults() {
    // Use a real scene file to verify the merge pipeline works,
    // but only check structural properties — not specific values.
    let config = load_scene_config(Path::new("scenes/cosmic-web-pm.toml")).unwrap();

    // Scene should override the base n_cells (32).
    assert_ne!(config.simulation.grid.n_cells, 32);

    // Display defaults from base.toml should be present.
    assert_eq!(config.output.display.point_size, 5.0);
    assert_eq!(config.output.display.camera_distance, 1.9);
    assert_eq!(config.output.display.gif_resolution, 512);
}
