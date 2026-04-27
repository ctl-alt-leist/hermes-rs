use hermes_rs::grid::Grid;

// ============================================================================
// Construction and derived quantities
// ============================================================================

#[test]
fn grid_cell_length() {
    let grid = Grid::new(64, 100_000.0);

    let expected = 100_000.0 / 64.0;
    assert!(
        (grid.cell_length - expected).abs() < 1e-10,
        "cell_length: expected {expected}, got {}",
        grid.cell_length
    );
}

#[test]
fn grid_total_cells() {
    let grid = Grid::new(32, 100_000.0);
    assert_eq!(grid.total_cells(), 32 * 32 * 32);
}

#[test]
fn grid_cell_volume() {
    let grid = Grid::new(64, 100_000.0);
    let h = grid.cell_length;
    let expected = h * h * h;

    assert!(
        (grid.cell_volume() - expected).abs() / expected < 1e-12,
        "cell_volume mismatch"
    );
}

#[test]
fn grid_box_volume() {
    let grid = Grid::new(64, 100_000.0);
    let expected = 100_000.0_f64.powi(3);

    assert!(
        (grid.box_volume() - expected).abs() / expected < 1e-12,
        "box_volume mismatch"
    );
}

// ============================================================================
// Periodic index wrapping
// ============================================================================

#[test]
fn wrap_index_identity() {
    let grid = Grid::new(64, 100_000.0);

    for m in 0..64 {
        assert_eq!(grid.wrap_index(m as isize), m);
    }
}

#[test]
fn wrap_index_negative() {
    let grid = Grid::new(64, 100_000.0);
    assert_eq!(grid.wrap_index(-1), 63);
    assert_eq!(grid.wrap_index(-2), 62);
    assert_eq!(grid.wrap_index(-64), 0);
}

#[test]
fn wrap_index_overflow() {
    let grid = Grid::new(64, 100_000.0);
    assert_eq!(grid.wrap_index(64), 0);
    assert_eq!(grid.wrap_index(65), 1);
    assert_eq!(grid.wrap_index(128), 0);
}

// ============================================================================
// Periodic position wrapping
// ============================================================================

#[test]
fn wrap_position_in_range() {
    let grid = Grid::new(64, 100.0);

    for x in [0.0, 25.0, 50.0, 99.99] {
        let wrapped = grid.wrap_position(x);
        assert!(
            (wrapped - x).abs() < 1e-10,
            "position {x} should be unchanged, got {wrapped}"
        );
    }
}

#[test]
fn wrap_position_overflow() {
    let grid = Grid::new(64, 100.0);

    let wrapped = grid.wrap_position(105.0);
    assert!(
        (wrapped - 5.0).abs() < 1e-10,
        "105.0 should wrap to 5.0, got {wrapped}"
    );
}

#[test]
fn wrap_position_negative() {
    let grid = Grid::new(64, 100.0);

    let wrapped = grid.wrap_position(-3.0);
    assert!(
        (wrapped - 97.0).abs() < 1e-10,
        "-3.0 should wrap to 97.0, got {wrapped}"
    );
}

#[test]
fn wrap_position_exact_boundary() {
    let grid = Grid::new(64, 100.0);

    // Exactly at box_length should wrap to 0
    let wrapped = grid.wrap_position(100.0);
    assert!(
        wrapped.abs() < 1e-10,
        "100.0 should wrap to 0.0, got {wrapped}"
    );
}

// ============================================================================
// Cell centers
// ============================================================================

#[test]
fn cell_center_first() {
    let grid = Grid::new(64, 100_000.0);
    let center = grid.cell_center(0, 0, 0);
    let half_h = grid.cell_length / 2.0;

    for d in 0..3 {
        assert!(
            (center[d] - half_h).abs() < 1e-10,
            "cell (0,0,0) center[{d}] should be {half_h}, got {}",
            center[d]
        );
    }
}

#[test]
fn cell_center_last() {
    let grid = Grid::new(4, 100.0);
    let center = grid.cell_center(3, 3, 3);
    let expected = (3.0 + 0.5) * 25.0;

    for d in 0..3 {
        assert!(
            (center[d] - expected).abs() < 1e-10,
            "cell (3,3,3) center[{d}] should be {expected}, got {}",
            center[d]
        );
    }
}

// ============================================================================
// Wavevector components
// ============================================================================

#[test]
fn wavevector_zero_mode() {
    let grid = Grid::new(64, 100_000.0);
    let k0 = grid.wavevector_component(0);

    assert!(k0.abs() < 1e-15, "k(0) should be 0, got {k0}");
}

#[test]
fn wavevector_fundamental() {
    let grid = Grid::new(64, 100_000.0);
    let k1 = grid.wavevector_component(1);
    let expected = 2.0 * std::f64::consts::PI / 100_000.0;

    assert!(
        (k1 - expected).abs() / expected < 1e-12,
        "k(1) should be {expected}, got {k1}"
    );
}

#[test]
fn wavevector_negative_frequency() {
    let grid = Grid::new(64, 100_000.0);

    // Index 63 should give the same magnitude as index 1, but negative
    let k_pos = grid.wavevector_component(1);
    let k_neg = grid.wavevector_component(63);

    assert!(
        (k_neg + k_pos).abs() < 1e-12,
        "k(63) should be -k(1): k(63) = {k_neg}, k(1) = {k_pos}"
    );
}
