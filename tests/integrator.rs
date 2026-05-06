use hermes_rs::algebra::vector_from_components;
use hermes_rs::physics::cic::ParticleForces;
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::initial::zeldovich_init;
use hermes_rs::physics::integrator::{drift, kick, midpoint, scale_factor_schedule, step_kdk};
use hermes_rs::physics::particles::Particles;
use hermes_rs::physics::poisson::PoissonSolver;
use morphis::vector::Vector;

// ============================================================================
// Kick
// ============================================================================

#[test]
fn kick_adds_force_times_factor() {
    let mut particles = Particles::zeros(1, 1.0);
    particles.set_momentum(0, &vector_from_components(1.0, 0.0, 0.0));

    // Constant force via ParticleForces.
    let mut force_data = ndarray::Array2::zeros((3, 1));
    force_data[[0, 0]] = 0.0;
    force_data[[1, 0]] = 10.0;
    force_data[[2, 0]] = 0.0;
    let forces = ParticleForces { data: force_data };

    let kick_factor = 0.5;
    kick(&mut particles, &forces, kick_factor);

    let momentum = particles.momentum_of(0);
    assert!((momentum.component(&[1]) - 1.0).abs() < 1e-12);
    assert!((momentum.component(&[2]) - 5.0).abs() < 1e-12);
    assert!((momentum.component(&[3]) - 0.0).abs() < 1e-12);
}

#[test]
fn kick_preserves_total_momentum_with_zero_net_force() {
    let grid = Grid::new(8, 80.0);
    let mut particles = Particles::on_lattice(8, &grid, 1e-7);

    // Give all particles the same momentum.
    for p in 0..particles.count() {
        particles.set_momentum(p, &vector_from_components(1.0, 2.0, 3.0));
    }

    let momentum_before = particles.total_momentum();

    // Zero net force.
    let forces = ParticleForces {
        data: ndarray::Array2::zeros((3, particles.count())),
    };

    kick(&mut particles, &forces, 1.0);

    let momentum_after = particles.total_momentum();
    let diff = &momentum_after - &momentum_before;

    assert!(
        diff.is_zero(1e-12),
        "zero force kick should preserve total momentum"
    );
}

// ============================================================================
// Drift
// ============================================================================

#[test]
fn drift_moves_particles_by_momentum() {
    let grid = Grid::new(8, 80.0);
    let mut particles = Particles::zeros(1, 2.0);

    particles.set_position(0, &vector_from_components(10.0, 20.0, 30.0));
    particles.set_momentum(0, &vector_from_components(4.0, 6.0, 8.0));

    let drift_factor = 0.5;
    drift(&mut particles, drift_factor, &grid);

    // x_new = x + (p / m) * drift_factor = (10, 20, 30) + (4, 6, 8) / 2 * 0.5
    //       = (10 + 1, 20 + 1.5, 30 + 2) = (11, 21.5, 32)
    let pos = particles.position_of(0);
    assert!((pos.component(&[1]) - 11.0).abs() < 1e-12);
    assert!((pos.component(&[2]) - 21.5).abs() < 1e-12);
    assert!((pos.component(&[3]) - 32.0).abs() < 1e-12);
}

#[test]
fn drift_wraps_periodic() {
    let grid = Grid::new(8, 100.0);
    let mut particles = Particles::zeros(1, 1.0);

    particles.set_position(0, &vector_from_components(99.0, 50.0, 50.0));
    particles.set_momentum(0, &vector_from_components(5.0, 0.0, 0.0));

    drift(&mut particles, 1.0, &grid);

    // 99 + 5 = 104, wraps to 4.
    let pos = particles.position_of(0);
    assert!(
        (pos.component(&[1]) - 4.0).abs() < 1e-10,
        "x should wrap: got {}",
        pos.component(&[1])
    );
}

#[test]
fn drift_with_zero_momentum_does_nothing() {
    let grid = Grid::new(8, 80.0);
    let mut particles = Particles::zeros(1, 1.0);
    particles.set_position(0, &vector_from_components(10.0, 20.0, 30.0));

    let pos_before = particles.position_of(0);
    drift(&mut particles, 1.0, &grid);
    let pos_after = particles.position_of(0);

    let diff = &pos_after - &pos_before;
    assert!(
        diff.is_zero(1e-15),
        "zero momentum drift should be stationary"
    );
}

// ============================================================================
// Time-reversal symmetry
// ============================================================================

#[test]
fn time_reversal_symmetry() {
    let grid = Grid::new(8, 80.0);
    let n_particles = 4;
    let mut particles = Particles::zeros(n_particles, 1.0);

    // Set up simple initial state with known positions and momenta.
    particles.set_position(0, &vector_from_components(10.0, 10.0, 10.0));
    particles.set_position(1, &vector_from_components(30.0, 30.0, 30.0));
    particles.set_position(2, &vector_from_components(50.0, 50.0, 50.0));
    particles.set_position(3, &vector_from_components(70.0, 70.0, 70.0));

    particles.set_momentum(0, &vector_from_components(1.0, 0.5, 0.0));
    particles.set_momentum(1, &vector_from_components(-0.5, 1.0, 0.3));
    particles.set_momentum(2, &vector_from_components(0.0, -1.0, 0.5));
    particles.set_momentum(3, &vector_from_components(-0.5, 0.5, -0.8));

    // Save initial state.
    let positions_initial: Vec<Vector<3>> =
        (0..n_particles).map(|p| particles.position_of(p)).collect();
    let momenta_initial: Vec<Vector<3>> =
        (0..n_particles).map(|p| particles.momentum_of(p)).collect();

    // Drift forward 10 steps.
    let n_steps = 10;
    let drift_factor = 0.1;
    for _ in 0..n_steps {
        drift(&mut particles, drift_factor, &grid);
    }

    // Negate all momenta (time reversal).
    for p in 0..n_particles {
        let neg_momentum = &particles.momentum_of(p) * -1.0;
        particles.set_momentum(p, &neg_momentum);
    }

    // Drift backward 10 steps.
    for _ in 0..n_steps {
        drift(&mut particles, drift_factor, &grid);
    }

    // Negate momenta again to restore original direction.
    for p in 0..n_particles {
        let neg_momentum = &particles.momentum_of(p) * -1.0;
        particles.set_momentum(p, &neg_momentum);
    }

    // Should return to initial state (within floating-point tolerance).
    for p in 0..n_particles {
        let pos = particles.position_of(p);
        let mom = particles.momentum_of(p);
        let pos_diff = &pos - &positions_initial[p];
        let mom_diff = &mom - &momenta_initial[p];

        assert!(
            pos_diff.is_zero(1e-10),
            "particle {p} position not restored after time reversal"
        );
        assert!(
            mom_diff.is_zero(1e-10),
            "particle {p} momentum not restored after time reversal"
        );
    }
}

// ============================================================================
// Scale factor schedule
// ============================================================================

#[test]
fn schedule_log_a_endpoints() {
    let schedule = scale_factor_schedule(0.02, 1.0, 100, "log_a");

    assert_eq!(schedule.len(), 101);
    assert!((schedule[0] - 0.02).abs() < 1e-12);
    assert!((schedule[100] - 1.0).abs() < 1e-12);
}

#[test]
fn schedule_log_a_monotone() {
    let schedule = scale_factor_schedule(0.02, 1.0, 50, "log_a");

    for n in 1..schedule.len() {
        assert!(
            schedule[n] > schedule[n - 1],
            "schedule must be monotonically increasing"
        );
    }
}

#[test]
fn schedule_linear_a_endpoints() {
    let schedule = scale_factor_schedule(0.1, 1.0, 10, "linear_a");

    assert_eq!(schedule.len(), 11);
    assert!((schedule[0] - 0.1).abs() < 1e-12);
    assert!((schedule[10] - 1.0).abs() < 1e-12);
}

#[test]
fn schedule_linear_a_uniform_spacing() {
    let schedule = scale_factor_schedule(0.0, 1.0, 10, "linear_a");
    let da = 0.1;

    for n in 0..=10 {
        let expected = n as f64 * da;
        assert!(
            (schedule[n] - expected).abs() < 1e-12,
            "step {n}: expected {expected}, got {}",
            schedule[n]
        );
    }
}

#[test]
fn midpoint_geometric_mean() {
    let mid = midpoint(0.04, 0.16);

    // sqrt(0.04 * 0.16) = sqrt(0.0064) = 0.08
    assert!(
        (mid - 0.08).abs() < 1e-12,
        "midpoint should be geometric mean: got {mid}"
    );
}

// ============================================================================
// Full KDK step (smoke test with real gravity)
// ============================================================================

#[test]
fn kdk_step_conserves_mass() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let mut particles = zeldovich_init(16, &grid, &cosmology, 0.02, 42).unwrap();
    let mass_before = particles.total_mass();

    let scale_factor_prev = 0.02;
    let scale_factor_next = 0.025;
    let scale_factor_mid = midpoint(scale_factor_prev, scale_factor_next);

    step_kdk(
        &mut particles,
        &mut solver,
        &grid,
        &cosmology,
        scale_factor_prev,
        scale_factor_mid,
        scale_factor_next,
        None,
    )
    .unwrap();

    let mass_after = particles.total_mass();

    assert!(
        (mass_after - mass_before).abs() / mass_before < 1e-12,
        "KDK step should conserve mass: before {mass_before}, after {mass_after}"
    );
}

#[test]
fn kdk_step_returns_forces_for_reuse() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let mut particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();

    let scale_factor_prev = 0.02;
    let scale_factor_next = 0.025;
    let scale_factor_mid = midpoint(scale_factor_prev, scale_factor_next);

    let forces = step_kdk(
        &mut particles,
        &mut solver,
        &grid,
        &cosmology,
        scale_factor_prev,
        scale_factor_mid,
        scale_factor_next,
        None,
    )
    .unwrap();

    assert_eq!(forces.count(), particles.count());

    // Force vectors should be grade-1.
    let force_0 = forces.force_on(0);
    assert_eq!(force_0.grade(), 1);
}

#[test]
fn kdk_two_steps_with_force_caching() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let mut particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();

    let schedule = scale_factor_schedule(0.02, 0.03, 2, "log_a");

    // Step 1: no cached forces.
    let mid_1 = midpoint(schedule[0], schedule[1]);
    let forces_1 = step_kdk(
        &mut particles,
        &mut solver,
        &grid,
        &cosmology,
        schedule[0],
        mid_1,
        schedule[1],
        None,
    )
    .unwrap();

    // Step 2: reuse forces from step 1.
    let mid_2 = midpoint(schedule[1], schedule[2]);
    let forces_2 = step_kdk(
        &mut particles,
        &mut solver,
        &grid,
        &cosmology,
        schedule[1],
        mid_2,
        schedule[2],
        Some(&forces_1),
    )
    .unwrap();

    assert_eq!(forces_2.count(), particles.count());
}
