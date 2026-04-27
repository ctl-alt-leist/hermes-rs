use hermes_rs::algebra::{components_from_vector, vector_from_components};
use hermes_rs::grid::Grid;
use hermes_rs::particles::Particles;

// ============================================================================
// Construction
// ============================================================================

#[test]
fn zeros_construction() {
    let particles = Particles::zeros(1000, 1e10);

    assert_eq!(particles.count(), 1000);
    assert!((particles.mass_particle - 1e10).abs() < 1e-5);
    assert_eq!(particles.position.shape(), &[3, 1000]);
    assert_eq!(particles.momentum.shape(), &[3, 1000]);
}

#[test]
fn total_mass() {
    let particles = Particles::zeros(100, 1e10);
    let expected = 100.0 * 1e10;

    assert!(
        (particles.total_mass() - expected).abs() / expected < 1e-12,
        "total_mass: expected {expected}, got {}",
        particles.total_mass()
    );
}

// ============================================================================
// Lattice placement
// ============================================================================

#[test]
fn lattice_particle_count() {
    let grid = Grid::new(64, 100_000.0);
    let density_mean = 1e-7;
    let particles = Particles::on_lattice(8, &grid, density_mean);

    assert_eq!(particles.count(), 8 * 8 * 8);
}

#[test]
fn lattice_mass_conservation() {
    let grid = Grid::new(64, 100_000.0);
    let density_mean = 1e-7;
    let particles = Particles::on_lattice(8, &grid, density_mean);

    let expected_total = density_mean * grid.box_volume();
    let computed_total = particles.total_mass();
    let rel_err = (computed_total - expected_total).abs() / expected_total;

    assert!(
        rel_err < 1e-12,
        "lattice total mass: expected {expected_total}, got {computed_total}"
    );
}

#[test]
fn lattice_positions_in_box() {
    let grid = Grid::new(64, 100_000.0);
    let particles = Particles::on_lattice(8, &grid, 1e-7);

    for n in 0..particles.count() {
        let pos = components_from_vector(&particles.position_of(n));
        for d in 0..3 {
            assert!(
                pos[d] >= 0.0 && pos[d] < grid.box_length,
                "particle {n} component {d} = {} is out of bounds [0, {})",
                pos[d],
                grid.box_length
            );
        }
    }
}

#[test]
fn lattice_uniform_spacing() {
    let grid = Grid::new(64, 100.0);
    let particles = Particles::on_lattice(4, &grid, 1e-7);
    let spacing = 100.0 / 4.0;

    let pos = particles.position_of(0);
    let c = components_from_vector(&pos);
    for d in 0..3 {
        assert!(
            (c[d] - spacing / 2.0).abs() < 1e-10,
            "first particle component {d} = {}, expected {}",
            c[d],
            spacing / 2.0
        );
    }
}

// ============================================================================
// Morphis-native position/momentum access
// ============================================================================

#[test]
fn set_and_get_position_morphis() {
    let mut particles = Particles::zeros(10, 1e10);

    let pos = vector_from_components(1.0, 2.0, 3.0);
    particles.set_position(5, &pos);
    let got = particles.position_of(5);

    assert_eq!(got.grade(), 1);
    assert!((got.component(&[0]) - 1.0).abs() < 1e-15);
    assert!((got.component(&[1]) - 2.0).abs() < 1e-15);
    assert!((got.component(&[2]) - 3.0).abs() < 1e-15);
}

#[test]
fn set_and_get_momentum_morphis() {
    let mut particles = Particles::zeros(10, 1e10);

    let mom = vector_from_components(4.0, 5.0, 6.0);
    particles.set_momentum(3, &mom);
    let got = particles.momentum_of(3);

    assert_eq!(got.grade(), 1);
    assert!((got.component(&[0]) - 4.0).abs() < 1e-15);
    assert!((got.component(&[1]) - 5.0).abs() < 1e-15);
    assert!((got.component(&[2]) - 6.0).abs() < 1e-15);
}

#[test]
fn position_of_returns_grade_1_vector() {
    let grid = Grid::new(8, 80.0);
    let particles = Particles::on_lattice(4, &grid, 1e-7);
    let pos = particles.position_of(0);

    assert_eq!(pos.grade(), 1);
    assert!(
        pos.norm() > 0.0,
        "lattice particle should have nonzero position"
    );
}

// ============================================================================
// Morphis-native derived quantities
// ============================================================================

#[test]
fn total_momentum_zero_for_stationary() {
    let particles = Particles::zeros(100, 1e10);
    let total = particles.total_momentum();

    assert!(
        total.is_zero(1e-15),
        "total momentum of stationary particles should be zero"
    );
}

#[test]
fn total_momentum_sums_correctly() {
    let mut particles = Particles::zeros(2, 1e10);
    particles.set_momentum(0, &vector_from_components(1.0, 0.0, 0.0));
    particles.set_momentum(1, &vector_from_components(0.0, 2.0, 3.0));

    let total = particles.total_momentum();
    assert!((total.component(&[0]) - 1.0).abs() < 1e-12);
    assert!((total.component(&[1]) - 2.0).abs() < 1e-12);
    assert!((total.component(&[2]) - 3.0).abs() < 1e-12);
}

#[test]
fn angular_momentum_is_bivector() {
    let mut particles = Particles::zeros(1, 1e10);
    particles.set_position(0, &vector_from_components(1.0, 0.0, 0.0));
    particles.set_momentum(0, &vector_from_components(0.0, 1.0, 0.0));

    let angular_momentum = particles.angular_momentum(0);

    // x ∧ p for x along e0, p along e1 → bivector in the e0∧e1 plane.
    assert_eq!(angular_momentum.grade(), 2);
    assert!((angular_momentum.component(&[0, 1]) - 1.0).abs() < 1e-12);
}

#[test]
fn total_angular_momentum_zero_for_symmetric_lattice() {
    let grid = Grid::new(8, 80.0);
    let particles = Particles::on_lattice(4, &grid, 1e-7);

    // Stationary particles on a symmetric lattice: L = Σ x ∧ 0 = 0.
    let total_angular_momentum = particles.total_angular_momentum();
    assert!(
        total_angular_momentum.is_zero(1e-15),
        "angular momentum should be zero for stationary particles"
    );
}

#[test]
fn kinetic_energy_matches_norm_squared() {
    let mut particles = Particles::zeros(1, 2.0);
    particles.set_momentum(0, &vector_from_components(3.0, 4.0, 0.0));

    let scale_factor = 1.0;
    let energy = particles.kinetic_energy(scale_factor);

    // E_k = |p|² / (2 m a²) = 25 / (2 * 2 * 1) = 6.25
    assert!(
        (energy - 6.25).abs() < 1e-12,
        "kinetic energy: expected 6.25, got {energy}"
    );
}

// ============================================================================
// Periodic wrapping
// ============================================================================

#[test]
fn wrap_positions_identity() {
    let grid = Grid::new(64, 100.0);
    let mut particles = Particles::on_lattice(4, &grid, 1e-7);
    let original = particles.position.clone();

    particles.wrap_positions(&grid);

    for n in 0..particles.count() {
        for d in 0..3 {
            assert!(
                (particles.position[[d, n]] - original[[d, n]]).abs() < 1e-15,
                "wrapping changed an in-bounds position"
            );
        }
    }
}

#[test]
fn wrap_positions_overflow() {
    let grid = Grid::new(64, 100.0);
    let mut particles = Particles::zeros(1, 1e10);

    particles.set_position_components(0, [105.0, -3.0, 250.0]);
    particles.wrap_positions(&grid);

    let pos = components_from_vector(&particles.position_of(0));
    assert!(
        (pos[0] - 5.0).abs() < 1e-10,
        "x: expected 5.0, got {}",
        pos[0]
    );
    assert!(
        (pos[1] - 97.0).abs() < 1e-10,
        "y: expected 97.0, got {}",
        pos[1]
    );
    assert!(
        (pos[2] - 50.0).abs() < 1e-10,
        "z: expected 50.0, got {}",
        pos[2]
    );
}
