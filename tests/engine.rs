/// Tests for the composable physics engine.
///
/// These tests verify that the sector-based Engine with merged Strang
/// splitting preserves conservation laws and produces correct dynamics.
use std::collections::BTreeMap;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric::euclidean;

use hermes_rs::engine::Engine;
use hermes_rs::engine::sector::Sector;
use hermes_rs::engine::sector::schrodinger::SchrodingerSector;
use hermes_rs::engine::solver::GravitySolver;
use hermes_rs::engine::state::{FieldEntry, SimulationState};
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;

// ============================================================================
// Helpers
// ============================================================================

fn make_state(
    n: usize,
    box_length: f64,
    alpha: EvenField<3>,
    ell: f64,
    mass: f64,
) -> (SimulationState, MorphisGrid<3>) {
    let grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: ell,
            mass,
            self_interaction: None,
        },
    );

    let state = SimulationState {
        particles: BTreeMap::new(),
        fields,
        grid,
        morphis_grid,
        time: 0.0,
        step: 0,
    };

    (state, morphis_grid)
}

fn field_norm(entry: &FieldEntry) -> f64 {
    entry
        .data
        .scalar
        .iter()
        .zip(entry.data.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum()
}

// ============================================================================
// Free evolution tests
// ============================================================================

#[test]
fn schrodinger_sector_preserves_norm() {
    let n = 16;
    let box_length = 1.0;
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();
    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (1.0, 0.0));
    let norm_before: f64 = alpha
        .scalar
        .iter()
        .zip(alpha.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum();

    let (state, _) = make_state(n, box_length, alpha, 1.0, 1.0);

    let sectors: Vec<Box<dyn Sector>> = vec![Box::new(SchrodingerSector::new("alpha".to_string()))];

    let mut engine = Engine::new(state, sectors, None, None);

    let dt = 0.01;
    for _ in 0..10 {
        engine.step(1.0, dt).unwrap();
    }
    engine.finalize(1.0, dt).unwrap();

    let norm_after = field_norm(&engine.state.fields["alpha"]);

    let relative_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        relative_error < 1e-12,
        "norm not preserved: relative error = {relative_error}"
    );
    assert_eq!(engine.state.step, 10);
}

#[test]
fn simulation_state_tracks_species() {
    let n = 8;
    let box_length = 1.0;
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();
    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (1.0, 0.0));

    let (state, _) = make_state(n, box_length, alpha, 1.0, 1.0);

    assert!(!state.has_particles());
    assert!(state.has_fields());
    assert_eq!(state.total_particle_count(), 0);
    assert!(state.fields.contains_key("alpha"));
}

// ============================================================================
// Gravity-coupled engine tests
// ============================================================================

#[test]
fn engine_gravity_coupled_field_preserves_norm() {
    let n = 16;
    let box_length = 10000.0;
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let cosmology = planck_2018();
    let density_mean = cosmology.density_matter();
    let mass = 1e10;
    let length_scale = 2000.0;
    let ell = length_scale * mass;

    let uniform = (density_mean / mass).sqrt();
    let alpha = EvenField::from_fn(&morphis_grid, g, |pos| {
        let r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
        let perturbation = 1.0 + 0.01 * (-r2 / (2000.0 * 2000.0)).exp();
        (uniform * perturbation, 0.0)
    });

    let norm_before = alpha
        .scalar
        .iter()
        .zip(alpha.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum::<f64>();

    let (state, morphis_grid) = make_state(n, box_length, alpha, ell, mass);

    let sectors: Vec<Box<dyn Sector>> = vec![Box::new(SchrodingerSector::new("alpha".to_string()))];
    let solver = GravitySolver::new(morphis_grid);

    let mut engine = Engine::new(state, sectors, Some(solver), Some(cosmology));

    let scale_factor = 0.5;
    let dt = 0.01;

    for _ in 0..20 {
        engine.step(scale_factor, dt).unwrap();
    }
    engine.finalize(scale_factor, dt).unwrap();

    let norm_after = field_norm(&engine.state.fields["alpha"]);

    let relative_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        relative_error < 1e-12,
        "norm not preserved under gravity: relative error = {relative_error}"
    );
}

#[test]
fn engine_gravity_coupled_field_runs_many_steps() {
    let n = 16;
    let box_length = 10000.0;
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let cosmology = planck_2018();
    let density_mean = cosmology.density_matter();
    let mass = 1e10;
    let ell = 2000.0 * mass;

    let uniform = (density_mean / mass).sqrt();
    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (uniform, 0.0));

    let norm_before: f64 = alpha
        .scalar
        .iter()
        .zip(alpha.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum();

    let (state, morphis_grid) = make_state(n, box_length, alpha, ell, mass);

    let sectors: Vec<Box<dyn Sector>> = vec![Box::new(SchrodingerSector::new("alpha".to_string()))];
    let solver = GravitySolver::new(morphis_grid);

    let mut engine = Engine::new(state, sectors, Some(solver), Some(cosmology));

    let dt = 0.005;
    for _ in 0..100 {
        engine.step(0.5, dt).unwrap();
    }
    engine.finalize(0.5, dt).unwrap();

    let norm_after = field_norm(&engine.state.fields["alpha"]);

    let relative_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        relative_error < 1e-12,
        "norm drifted over 100 steps: relative error = {relative_error}"
    );
    assert_eq!(engine.state.step, 100);
}

// ============================================================================
// Two-sector tests
// ============================================================================

#[test]
fn two_sector_gravity_preserves_mass_independently() {
    use hermes_rs::engine::sector::gross_pitaevskii::GrossPitaevskiiSector;

    let n = 16;
    let box_length = 10000.0;
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let cosmology = planck_2018();
    let density_mean = cosmology.density_matter();
    let mass = 1e10;
    let ell = 2000.0 * mass;

    // Alpha: 84% of matter (dark matter), uniform with perturbation.
    let fraction_alpha = 0.844;
    let uniform_alpha = (fraction_alpha * density_mean / mass).sqrt();
    let alpha = EvenField::from_fn(&morphis_grid, g, |pos| {
        let r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
        let perturbation = 1.0 + 0.01 * (-r2 / (2000.0 * 2000.0)).exp();
        (uniform_alpha * perturbation, 0.0)
    });

    // Beta: 16% of matter (baryons), uniform with offset perturbation.
    let fraction_beta = 0.156;
    let uniform_beta = (fraction_beta * density_mean / mass).sqrt();
    let beta = EvenField::from_fn(&morphis_grid, g, |pos| {
        let dx = pos[0] - box_length / 4.0;
        let r2 = dx * dx + pos[1] * pos[1] + pos[2] * pos[2];
        let perturbation = 1.0 + 0.02 * (-r2 / (1500.0 * 1500.0)).exp();
        (uniform_beta * perturbation, 0.0)
    });

    let norm_alpha_before = field_norm_raw(&alpha);
    let norm_beta_before = field_norm_raw(&beta);

    let grid = Grid::new(n, box_length);

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: ell,
            mass,
            self_interaction: None,
        },
    );
    fields.insert(
        "beta".to_string(),
        FieldEntry {
            data: beta,
            smoothing_length: ell,
            mass,
            self_interaction: Some(1e6),
        },
    );

    let state = SimulationState {
        particles: BTreeMap::new(),
        fields,
        grid,
        morphis_grid,
        time: 0.0,
        step: 0,
    };

    let sectors: Vec<Box<dyn Sector>> = vec![
        Box::new(SchrodingerSector::new("alpha".to_string())),
        Box::new(GrossPitaevskiiSector::new("beta".to_string())),
    ];
    let solver = GravitySolver::new(morphis_grid);

    let mut engine = Engine::new(state, sectors, Some(solver), Some(cosmology));

    let dt = 0.005;
    for _ in 0..50 {
        engine.step(0.5, dt).unwrap();
    }
    engine.finalize(0.5, dt).unwrap();

    // Each sector's mass must be independently conserved.
    let norm_alpha_after = field_norm(&engine.state.fields["alpha"]);
    let norm_beta_after = field_norm(&engine.state.fields["beta"]);

    let error_alpha = (norm_alpha_after - norm_alpha_before).abs() / norm_alpha_before;
    let error_beta = (norm_beta_after - norm_beta_before).abs() / norm_beta_before;

    assert!(
        error_alpha < 1e-12,
        "α norm drifted: relative error = {error_alpha}"
    );
    assert!(
        error_beta < 1e-12,
        "β norm drifted: relative error = {error_beta}"
    );
    assert_eq!(engine.state.step, 50);
}

#[test]
fn two_sector_cross_sourcing_differs_from_isolation() {
    use hermes_rs::engine::sector::gross_pitaevskii::GrossPitaevskiiSector;

    let n = 16;
    let box_length = 10000.0;
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let cosmology = planck_2018();
    let density_mean = cosmology.density_matter();
    let mass = 1e10;
    let ell = 2000.0 * mass;

    // A lopsided perturbation in beta that should gravitationally
    // affect alpha differently than alpha's own self-gravity alone.
    let uniform_alpha = (0.844 * density_mean / mass).sqrt();
    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (uniform_alpha, 0.0));

    let uniform_beta = (0.156 * density_mean / mass).sqrt();
    let beta = EvenField::from_fn(&morphis_grid, g, |pos| {
        let dx = pos[0] - box_length / 4.0;
        let r2 = dx * dx + pos[1] * pos[1] + pos[2] * pos[2];
        let perturbation = 1.0 + 0.1 * (-r2 / (1000.0 * 1000.0)).exp();
        (uniform_beta * perturbation, 0.0)
    });

    // Run alpha alone under its own gravity.
    let grid = Grid::new(n, box_length);

    let mut fields_alone = BTreeMap::new();
    fields_alone.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha.clone(),
            smoothing_length: ell,
            mass,
            self_interaction: None,
        },
    );

    let state_alone = SimulationState {
        particles: BTreeMap::new(),
        fields: fields_alone,
        grid: grid.clone(),
        morphis_grid,
        time: 0.0,
        step: 0,
    };

    let sectors_alone: Vec<Box<dyn Sector>> =
        vec![Box::new(SchrodingerSector::new("alpha".to_string()))];
    let solver_alone = GravitySolver::new(morphis_grid);

    let mut engine_alone = Engine::new(
        state_alone,
        sectors_alone,
        Some(solver_alone),
        Some(cosmology.clone()),
    );

    let dt = 0.005;
    let n_steps = 20;

    for _ in 0..n_steps {
        engine_alone.step(0.5, dt).unwrap();
    }
    engine_alone.finalize(0.5, dt).unwrap();

    // Run alpha + beta together under shared gravity.
    let mut fields_joint = BTreeMap::new();
    fields_joint.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: ell,
            mass,
            self_interaction: None,
        },
    );
    fields_joint.insert(
        "beta".to_string(),
        FieldEntry {
            data: beta,
            smoothing_length: ell,
            mass,
            self_interaction: Some(1e6),
        },
    );

    let state_joint = SimulationState {
        particles: BTreeMap::new(),
        fields: fields_joint,
        grid,
        morphis_grid,
        time: 0.0,
        step: 0,
    };

    let sectors_joint: Vec<Box<dyn Sector>> = vec![
        Box::new(SchrodingerSector::new("alpha".to_string())),
        Box::new(GrossPitaevskiiSector::new("beta".to_string())),
    ];
    let solver_joint = GravitySolver::new(morphis_grid);

    let mut engine_joint = Engine::new(
        state_joint,
        sectors_joint,
        Some(solver_joint),
        Some(cosmology),
    );

    for _ in 0..n_steps {
        engine_joint.step(0.5, dt).unwrap();
    }
    engine_joint.finalize(0.5, dt).unwrap();

    // Alpha's final state should differ between the two runs,
    // because in the joint run beta's density perturbation contributes
    // to the gravitational potential that alpha feels.
    let alpha_alone = &engine_alone.state.fields["alpha"].data;
    let alpha_joint = &engine_joint.state.fields["alpha"].data;

    let diff: f64 = alpha_alone
        .scalar
        .iter()
        .zip(alpha_joint.scalar.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();

    assert!(
        diff > 1e-15,
        "α should evolve differently with β present, but diff = {diff}"
    );
}

/// Norm computed directly from an EvenField (before it's in a FieldEntry).
fn field_norm_raw(field: &EvenField<3>) -> f64 {
    field
        .scalar
        .iter()
        .zip(field.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum()
}

// ============================================================================
// Snapshot capture tests
// ============================================================================

#[test]
fn snapshot_captures_all_field_species() {
    let n = 8;
    let box_length = 1.0;
    let grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (1.0, 0.0));
    let beta = EvenField::from_fn(&morphis_grid, g, |_| (0.5, 0.0));

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: 1.0,
            mass: 1.0,
            self_interaction: None,
        },
    );
    fields.insert(
        "beta".to_string(),
        FieldEntry {
            data: beta,
            smoothing_length: 1.0,
            mass: 1.0,
            self_interaction: Some(1.0),
        },
    );

    let state = SimulationState {
        particles: BTreeMap::new(),
        fields,
        grid,
        morphis_grid,
        time: 0.0,
        step: 5,
    };

    let snapshot = hermes_rs::io::snapshot::Snapshot::capture_from_state(&state, 5, 1.0);

    assert_eq!(snapshot.fields.len(), 2);
    assert!(snapshot.fields.iter().any(|f| f.name == "alpha"));
    assert!(snapshot.fields.iter().any(|f| f.name == "beta"));
    assert!(snapshot.particles.is_empty());
    assert_eq!(snapshot.n_cells, n);
}

#[test]
fn snapshot_captures_all_particle_species() {
    use hermes_rs::physics::particles::Particles;

    let n = 8;
    let box_length = 100.0;
    let grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);

    let dm = Particles::zeros(10, 1e10);
    let stars = Particles::zeros(5, 1e8);

    let mut particles = BTreeMap::new();
    particles.insert("dark_matter".to_string(), dm);
    particles.insert("stars".to_string(), stars);

    let state = SimulationState {
        particles,
        fields: BTreeMap::new(),
        grid,
        morphis_grid,
        time: 0.0,
        step: 3,
    };

    let snapshot = hermes_rs::io::snapshot::Snapshot::capture_from_state(&state, 3, 0.5);

    assert_eq!(snapshot.particles.len(), 2);
    assert!(snapshot.particles.iter().any(|p| p.name == "dark_matter"));
    assert!(snapshot.particles.iter().any(|p| p.name == "stars"));
    assert_eq!(snapshot.particle_count(), 15);
    assert!(snapshot.fields.is_empty());
}

#[test]
fn snapshot_captures_mixed_content() {
    use hermes_rs::physics::particles::Particles;

    let n = 8;
    let box_length = 100.0;
    let grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (1.0, 0.0));
    let dm = Particles::zeros(10, 1e10);

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: 1.0,
            mass: 1.0,
            self_interaction: None,
        },
    );

    let mut particles = BTreeMap::new();
    particles.insert("dark_matter".to_string(), dm);

    let state = SimulationState {
        particles,
        fields,
        grid,
        morphis_grid,
        time: 0.0,
        step: 0,
    };

    let snapshot = hermes_rs::io::snapshot::Snapshot::capture_from_state(&state, 0, 1.0);

    assert!(snapshot.has_particles());
    assert!(snapshot.has_fields());
    assert_eq!(snapshot.particles.len(), 1);
    assert_eq!(snapshot.fields.len(), 1);
    assert_eq!(snapshot.particle_count(), 10);
    assert_eq!(snapshot.n_cells, n);
}
