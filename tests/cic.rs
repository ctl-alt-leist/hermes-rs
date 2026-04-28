use hermes_rs::algebra::vector_from_array;
use hermes_rs::physics::cic::{assign_density, interpolate_force};
use hermes_rs::physics::field::VectorField;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::particles::Particles;

// ============================================================================
// Mass conservation
// ============================================================================

#[test]
fn total_mass_conserved() {
    let grid = Grid::new(16, 100.0);
    let particles = Particles::on_lattice(8, &grid, 1e-7);

    let density = assign_density(&particles, &grid);
    let total_deposited = density.sum() * grid.cell_volume();
    let total_particle = particles.total_mass();
    let rel_err = (total_deposited - total_particle).abs() / total_particle;

    assert!(
        rel_err < 1e-12,
        "total mass not conserved: deposited {total_deposited}, particles {total_particle}, rel_err {rel_err}"
    );
}

#[test]
fn total_mass_conserved_random_positions() {
    let grid = Grid::new(16, 100.0);
    let n_particles = 500;
    let mass_particle = 1e8;
    let mut particles = Particles::zeros(n_particles, mass_particle);

    // Deterministic pseudo-random positions using a simple LCG.
    let mut seed: u64 = 12345;
    for p in 0..n_particles {
        for d in 0..3 {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let frac = (seed >> 33) as f64 / (1u64 << 31) as f64;
            particles.position[[d, p]] = frac * grid.box_length;
        }
    }

    let density = assign_density(&particles, &grid);
    let total_deposited = density.sum() * grid.cell_volume();
    let total_particle = particles.total_mass();
    let rel_err = (total_deposited - total_particle).abs() / total_particle;

    assert!(
        rel_err < 1e-12,
        "mass not conserved for random positions: rel_err {rel_err}"
    );
}

// ============================================================================
// Single particle placement
// ============================================================================

#[test]
fn particle_at_cell_center() {
    let grid = Grid::new(8, 80.0);
    let mut particles = Particles::zeros(1, 1e10);

    // Place particle at the center of cell (3, 3, 3).
    let center = grid.cell_center_vector(3, 3, 3);
    particles.set_position(0, &center);

    let density = assign_density(&particles, &grid);

    // At the exact cell center, all mass should land in that one cell.
    let expected_density = particles.mass_particle / grid.cell_volume();
    let computed = density.get(3, 3, 3);
    let rel_err = (computed - expected_density).abs() / expected_density;

    assert!(
        rel_err < 1e-12,
        "cell-center particle: expected density {expected_density}, got {computed}"
    );

    // All other cells should be zero.
    let total = density.sum() * grid.cell_volume();
    let rel_err_total = (total - particles.mass_particle).abs() / particles.mass_particle;
    assert!(
        rel_err_total < 1e-12,
        "total mass check failed: {rel_err_total}"
    );
}

#[test]
fn particle_at_cell_corner() {
    let grid = Grid::new(8, 80.0);
    let h = grid.cell_length;
    let mut particles = Particles::zeros(1, 1e10);

    // Place particle at the corner shared by cells (2,2,2), (3,2,2), etc.
    // This is position (3h, 3h, 3h) — the vertex between cells 2 and 3.
    particles.set_position_components(0, [3.0 * h, 3.0 * h, 3.0 * h]);

    let density = assign_density(&particles, &grid);

    // At a corner, mass splits equally among 8 cells (1/8 each).
    let expected_per_cell = particles.mass_particle / (8.0 * grid.cell_volume());

    for dm0 in 0..2_usize {
        for dm1 in 0..2_usize {
            for dm2 in 0..2_usize {
                let density_cell = density.get(2 + dm0, 2 + dm1, 2 + dm2);
                let rel_err = (density_cell - expected_per_cell).abs() / expected_per_cell;
                assert!(
                    rel_err < 1e-12,
                    "corner cell ({}, {}, {}): expected {expected_per_cell}, got {density_cell}",
                    2 + dm0,
                    2 + dm1,
                    2 + dm2
                );
            }
        }
    }
}

#[test]
fn particle_at_cell_face_center() {
    let grid = Grid::new(8, 80.0);
    let h = grid.cell_length;
    let mut particles = Particles::zeros(1, 1e10);

    // Place at the face center between cells (3,3,3) and (4,3,3).
    // Position: x = 4h (face), y = 3.5h (center), z = 3.5h (center).
    particles.set_position_components(0, [4.0 * h, 3.5 * h, 3.5 * h]);

    let density = assign_density(&particles, &grid);

    // Mass splits 50/50 between cells (3,3,3) and (4,3,3).
    let expected = particles.mass_particle / (2.0 * grid.cell_volume());

    let density_left = density.get(3, 3, 3);
    let density_right = density.get(4, 3, 3);

    assert!(
        (density_left - expected).abs() / expected < 1e-12,
        "face left: expected {expected}, got {density_left}"
    );
    assert!(
        (density_right - expected).abs() / expected < 1e-12,
        "face right: expected {expected}, got {density_right}"
    );
}

// ============================================================================
// Periodic boundary wrapping
// ============================================================================

#[test]
fn particle_near_boundary_wraps() {
    let grid = Grid::new(8, 80.0);
    let h = grid.cell_length;
    let mut particles = Particles::zeros(1, 1e10);

    // Place near the upper boundary: just past cell 7, wraps to cell 0.
    particles.set_position_components(0, [7.5 * h + 0.5 * h, 3.5 * h, 3.5 * h]);

    let density = assign_density(&particles, &grid);

    // Should deposit in cells 7 and 0 (periodic wrap).
    let total = density.sum() * grid.cell_volume();
    let rel_err = (total - particles.mass_particle).abs() / particles.mass_particle;
    assert!(
        rel_err < 1e-12,
        "boundary wrap: mass not conserved, rel_err {rel_err}"
    );

    // Cells 7 and 0 should have nonzero density.
    assert!(
        density.get(7, 3, 3) > 0.0 || density.get(0, 3, 3) > 0.0,
        "boundary particle should deposit in cells 7 and/or 0"
    );
}

// ============================================================================
// Force interpolation
// ============================================================================

#[test]
fn uniform_force_returns_same_for_all_particles() {
    let grid = Grid::new(8, 80.0);
    let particles = Particles::on_lattice(4, &grid, 1e-7);

    // Uniform acceleration field: a = (1.0, 2.0, 3.0) everywhere.
    // Force = m_p * acceleration.
    let mut accel_field = VectorField::zeros(&grid);
    accel_field.data[0].fill(1.0);
    accel_field.data[1].fill(2.0);
    accel_field.data[2].fill(3.0);

    let forces = interpolate_force(&accel_field, &particles, &grid);
    let mass_particle = particles.mass_particle;

    for p in 0..particles.count() {
        let force_p = forces.force_on(p);
        assert!(
            (force_p.component(&[0]) - 1.0 * mass_particle).abs() / mass_particle < 1e-12,
            "particle {p}: Fx = {}, expected {}",
            force_p.component(&[0]),
            1.0 * mass_particle
        );
        assert!(
            (force_p.component(&[1]) - 2.0 * mass_particle).abs() / mass_particle < 1e-12,
            "particle {p}: Fy = {}, expected {}",
            force_p.component(&[1]),
            2.0 * mass_particle
        );
        assert!(
            (force_p.component(&[2]) - 3.0 * mass_particle).abs() / mass_particle < 1e-12,
            "particle {p}: Fz = {}, expected {}",
            force_p.component(&[2]),
            3.0 * mass_particle
        );
    }
}

#[test]
fn interpolation_at_cell_center_reads_cell_value() {
    let grid = Grid::new(8, 80.0);
    let mut particles = Particles::zeros(1, 1e10);
    let center = grid.cell_center(3, 3, 3);
    particles.set_position(0, &vector_from_array(center));

    let mut accel_field = VectorField::zeros(&grid);
    accel_field.data[0][[3, 3, 3]] = 7.0;
    accel_field.data[1][[3, 3, 3]] = 8.0;
    accel_field.data[2][[3, 3, 3]] = 9.0;

    let forces = interpolate_force(&accel_field, &particles, &grid);
    let force_0 = forces.force_on(0);
    let mass_particle = particles.mass_particle;

    assert!(
        (force_0.component(&[0]) - 7.0 * mass_particle).abs() < 1e-2,
        "Fx at cell center: expected {}, got {}",
        7.0 * mass_particle,
        force_0.component(&[0])
    );
    assert!(
        (force_0.component(&[1]) - 8.0 * mass_particle).abs() < 1e-2,
        "Fy at cell center: expected {}, got {}",
        8.0 * mass_particle,
        force_0.component(&[1])
    );
    assert!(
        (force_0.component(&[2]) - 9.0 * mass_particle).abs() < 1e-2,
        "Fz at cell center: expected {}, got {}",
        9.0 * mass_particle,
        force_0.component(&[2])
    );
}

// ============================================================================
// Density field properties
// ============================================================================

#[test]
fn density_non_negative() {
    let grid = Grid::new(16, 100.0);
    let particles = Particles::on_lattice(8, &grid, 1e-7);
    let density = assign_density(&particles, &grid);

    for &density_cell in density.data.iter() {
        assert!(
            density_cell >= 0.0,
            "density must be non-negative, got {density_cell}"
        );
    }
}

#[test]
fn lattice_produces_uniform_density() {
    let grid = Grid::new(8, 80.0);
    let density_mean = 1e-7;
    let particles = Particles::on_lattice(8, &grid, density_mean);
    let density = assign_density(&particles, &grid);

    // With particles on a lattice matching the grid, every cell gets
    // exactly one particle's worth of mass.
    for m0 in 0..8 {
        for m1 in 0..8 {
            for m2 in 0..8 {
                let density_cell = density.get(m0, m1, m2);
                let rel_err = (density_cell - density_mean).abs() / density_mean;
                assert!(
                    rel_err < 1e-12,
                    "cell ({m0},{m1},{m2}): expected {density_mean}, got {density_cell}"
                );
            }
        }
    }
}
