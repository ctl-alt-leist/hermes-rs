use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::diagnostics::Diagnostics;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::initial::zeldovich_init;
use hermes_rs::physics::particles::Particles;
use hermes_rs::physics::poisson::PoissonSolver;

// ============================================================================
// Stationary particles
// ============================================================================

#[test]
fn stationary_particles_zero_momentum() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = Particles::on_lattice(8, &grid, cosmology.density_matter());
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 1.0);

    assert!(
        diag.momentum_total.is_zero(1e-15),
        "stationary particles should have zero total momentum"
    );
}

#[test]
fn stationary_particles_zero_kinetic_energy() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = Particles::on_lattice(8, &grid, cosmology.density_matter());
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 1.0);

    assert!(
        diag.energy_kinetic.abs() < 1e-15,
        "stationary particles should have zero kinetic energy, got {}",
        diag.energy_kinetic
    );
}

#[test]
fn stationary_particles_zero_angular_momentum() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = Particles::on_lattice(8, &grid, cosmology.density_matter());
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 1.0);

    assert!(
        diag.angular_momentum.is_zero(1e-15),
        "stationary particles should have zero angular momentum"
    );
}

// ============================================================================
// Mass conservation
// ============================================================================

#[test]
fn mass_equals_particle_count_times_mass() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = Particles::on_lattice(8, &grid, cosmology.density_matter());
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 1.0);

    let expected = particles.total_mass();
    assert!(
        (diag.mass_total - expected).abs() / expected < 1e-12,
        "mass_total should be N_p × m_p"
    );
}

// ============================================================================
// Grade checks on morphis objects
// ============================================================================

#[test]
fn momentum_is_grade_1() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 0.02);

    assert_eq!(diag.momentum_total.grade(), 1);
}

#[test]
fn angular_momentum_is_grade_2() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 0.02);

    assert_eq!(diag.angular_momentum.grade(), 2);
}

// ============================================================================
// Energy
// ============================================================================

#[test]
fn zeldovich_has_nonzero_kinetic_energy() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 0.02);

    assert!(
        diag.energy_kinetic > 0.0,
        "Zel'dovich ICs should have nonzero kinetic energy"
    );
}

#[test]
fn potential_energy_is_negative() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = zeldovich_init(16, &grid, &cosmology, 0.02, 42).unwrap();
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 0.02);

    assert!(
        diag.energy_potential < 0.0,
        "gravitational potential energy should be negative, got {}",
        diag.energy_potential
    );
}

// ============================================================================
// Derived quantities
// ============================================================================

#[test]
fn energy_total_is_sum() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 0.02);

    let expected = diag.energy_kinetic + diag.energy_potential;
    assert!(
        (diag.energy_total() - expected).abs() < 1e-15,
        "energy_total should be kinetic + potential"
    );
}

#[test]
fn momentum_magnitude_consistent() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);

    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();
    let diag = Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, 0.02);

    let norm_direct = diag.momentum_total.norm();
    assert!(
        (diag.momentum_magnitude() - norm_direct).abs() < 1e-15,
        "momentum_magnitude should match norm()"
    );
}
