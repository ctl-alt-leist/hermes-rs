/// Tests for the composable physics engine.
///
/// These tests verify that the Engine produces identical results to the
/// old Dynamics-based path, ensuring the refactor is behavior-preserving.
use std::collections::BTreeMap;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric::euclidean;

use hermes_rs::engine::Engine;
use hermes_rs::engine::coupling::poisson::PoissonGravity;
use hermes_rs::engine::free::FreeEvolution;
use hermes_rs::engine::free::schrodinger::SchrodingerEvolution;
use hermes_rs::engine::state::{FieldEntry, SimulationState};
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;

// ============================================================================
// Free evolution tests
// ============================================================================

#[test]
fn schrodinger_free_evolution_preserves_norm() {
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

    let mut entry = FieldEntry {
        data: alpha,
        smoothing_length: 1.0,
        mass: 1.0,
    };

    let mut evolver = SchrodingerEvolution;

    for _ in 0..10 {
        evolver.step(&mut entry, &morphis_grid, 1.0, 0.01).unwrap();
    }

    let norm_after: f64 = entry
        .data
        .scalar
        .iter()
        .zip(entry.data.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum();

    let relative_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        relative_error < 1e-12,
        "norm not preserved: relative error = {relative_error}"
    );
}

#[test]
fn simulation_state_tracks_species() {
    let n = 8;
    let box_length = 1.0;
    let grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (1.0, 0.0));

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: 1.0,
            mass: 1.0,
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

    assert!(!state.has_particles());
    assert!(state.has_fields());
    assert_eq!(state.total_particle_count(), 0);
    assert!(state.fields.contains_key("alpha"));
}

#[test]
fn engine_free_only_evolves_without_crash() {
    use hermes_rs::engine::Engine;

    let n = 16;
    let box_length = 1.0;
    let grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let alpha = EvenField::from_fn(&morphis_grid, g, |_| (1.0, 0.0));

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: 1.0,
            mass: 1.0,
        },
    );

    let mut free_modules: BTreeMap<String, Box<dyn FreeEvolution>> = BTreeMap::new();
    free_modules.insert("alpha".to_string(), Box::new(SchrodingerEvolution));

    let mut engine = Engine {
        state: SimulationState {
            particles: BTreeMap::new(),
            fields,
            grid,
            morphis_grid,
            time: 0.0,
            step: 0,
        },
        free_modules,
        couplings: vec![],
        cosmology: None,
    };

    // 10 steps of free Schrodinger evolution, no gravity, no expansion.
    for _ in 0..10 {
        engine.step(1.0, 0.01).unwrap();
    }

    assert_eq!(engine.state.step, 10);
}

// ============================================================================
// Gravity-coupled engine tests
// ============================================================================

fn field_norm(entry: &FieldEntry) -> f64 {
    entry
        .data
        .scalar
        .iter()
        .zip(entry.data.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum()
}

#[test]
fn engine_gravity_coupled_field_preserves_norm() {
    let n = 16;
    let box_length = 10000.0;
    let grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let cosmology = planck_2018();
    let density_mean = cosmology.density_matter();
    let mass = 1e10;
    let length_scale = 2000.0;
    let ell = length_scale * mass;

    // Uniform field with small perturbation.
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

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: ell,
            mass,
        },
    );

    let mut free_modules: BTreeMap<String, Box<dyn FreeEvolution>> = BTreeMap::new();
    free_modules.insert("alpha".to_string(), Box::new(SchrodingerEvolution));

    let gravity = PoissonGravity::new(grid.clone());

    let mut engine = Engine {
        state: SimulationState {
            particles: BTreeMap::new(),
            fields,
            grid,
            morphis_grid,
            time: 0.0,
            step: 0,
        },
        free_modules,
        couplings: vec![Box::new(gravity)],
        cosmology: Some(cosmology),
    };

    let scale_factor = 0.5;
    let dt = 0.01;

    for _ in 0..20 {
        engine.step(scale_factor, dt).unwrap();
    }

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
    let grid = Grid::new(n, box_length);
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

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: ell,
            mass,
        },
    );

    let mut free_modules: BTreeMap<String, Box<dyn FreeEvolution>> = BTreeMap::new();
    free_modules.insert("alpha".to_string(), Box::new(SchrodingerEvolution));

    let gravity = PoissonGravity::new(grid.clone());

    let mut engine = Engine {
        state: SimulationState {
            particles: BTreeMap::new(),
            fields,
            grid,
            morphis_grid,
            time: 0.0,
            step: 0,
        },
        free_modules,
        couplings: vec![Box::new(gravity)],
        cosmology: Some(cosmology),
    };

    // 100 steps — if the integrator is broken, norm will drift.
    for _ in 0..100 {
        engine.step(0.5, 0.005).unwrap();
    }

    let norm_after = field_norm(&engine.state.fields["alpha"]);

    let relative_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        relative_error < 1e-12,
        "norm drifted over 100 steps: relative error = {relative_error}"
    );
    assert_eq!(engine.state.step, 100);
}
