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
        let pos = particles.position_of(n);
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

    // First particle should be at (spacing/2, spacing/2, spacing/2)
    let pos = particles.position_of(0);
    for d in 0..3 {
        assert!(
            (pos[d] - spacing / 2.0).abs() < 1e-10,
            "first particle component {d} = {}, expected {}",
            pos[d],
            spacing / 2.0
        );
    }
}

// ============================================================================
// Position/momentum access
// ============================================================================

#[test]
fn set_and_get_position() {
    let mut particles = Particles::zeros(10, 1e10);

    particles.set_position(5, [1.0, 2.0, 3.0]);
    let pos = particles.position_of(5);

    assert!((pos[0] - 1.0).abs() < 1e-15);
    assert!((pos[1] - 2.0).abs() < 1e-15);
    assert!((pos[2] - 3.0).abs() < 1e-15);
}

#[test]
fn set_and_get_momentum() {
    let mut particles = Particles::zeros(10, 1e10);

    particles.set_momentum(3, [4.0, 5.0, 6.0]);
    let mom = particles.momentum_of(3);

    assert!((mom[0] - 4.0).abs() < 1e-15);
    assert!((mom[1] - 5.0).abs() < 1e-15);
    assert!((mom[2] - 6.0).abs() < 1e-15);
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

    // All positions already in [0, 100), so wrapping should change nothing.
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

    particles.set_position(0, [105.0, -3.0, 250.0]);
    particles.wrap_positions(&grid);

    let pos = particles.position_of(0);
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
