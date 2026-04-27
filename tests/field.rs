use hermes_rs::field::{ScalarField, VectorField};
use hermes_rs::grid::Grid;

// ============================================================================
// ScalarField
// ============================================================================

#[test]
fn scalar_field_zeros() {
    let grid = Grid::new(8, 100.0);
    let field = ScalarField::zeros(&grid);

    assert_eq!(field.data.shape(), &[8, 8, 8]);
    assert!(field.sum().abs() < 1e-15);
}

#[test]
fn scalar_field_sum() {
    let grid = Grid::new(4, 100.0);
    let mut field = ScalarField::zeros(&grid);

    *field.get_mut(0, 0, 0) = 1.0;
    *field.get_mut(1, 2, 3) = 2.5;

    assert!((field.sum() - 3.5).abs() < 1e-15);
}

#[test]
fn scalar_field_scale() {
    let grid = Grid::new(4, 100.0);
    let mut field = ScalarField::zeros(&grid);
    *field.get_mut(0, 0, 0) = 3.0;

    field.scale(2.0);

    assert!((field.get(0, 0, 0) - 6.0).abs() < 1e-15);
}

#[test]
fn scalar_field_add() {
    let grid = Grid::new(4, 100.0);
    let mut a = ScalarField::zeros(&grid);
    let mut b = ScalarField::zeros(&grid);

    *a.get_mut(0, 0, 0) = 1.0;
    *b.get_mut(0, 0, 0) = 2.0;

    let c = &a + &b;
    assert!((c.get(0, 0, 0) - 3.0).abs() < 1e-15);
}

#[test]
fn scalar_field_sub() {
    let grid = Grid::new(4, 100.0);
    let mut a = ScalarField::zeros(&grid);
    let mut b = ScalarField::zeros(&grid);

    *a.get_mut(1, 1, 1) = 5.0;
    *b.get_mut(1, 1, 1) = 2.0;

    let c = &a - &b;
    assert!((c.get(1, 1, 1) - 3.0).abs() < 1e-15);
}

#[test]
fn scalar_field_mul_scalar() {
    let grid = Grid::new(4, 100.0);
    let mut a = ScalarField::zeros(&grid);
    *a.get_mut(2, 2, 2) = 4.0;

    let b = &a * 0.5;
    assert!((b.get(2, 2, 2) - 2.0).abs() < 1e-15);
}

// ============================================================================
// VectorField
// ============================================================================

#[test]
fn vector_field_zeros() {
    let grid = Grid::new(8, 100.0);
    let field = VectorField::zeros(&grid);

    for d in 0..3 {
        assert_eq!(field.component(d).shape(), &[8, 8, 8]);
        assert!(field.component(d).sum().abs() < 1e-15);
    }
}

#[test]
fn vector_field_component_access() {
    let grid = Grid::new(4, 100.0);
    let mut field = VectorField::zeros(&grid);

    field.component_mut(0)[[1, 2, 3]] = 7.0;
    field.component_mut(1)[[1, 2, 3]] = 8.0;
    field.component_mut(2)[[1, 2, 3]] = 9.0;

    assert!((field.get(0, 1, 2, 3) - 7.0).abs() < 1e-15);
    assert!((field.get(1, 1, 2, 3) - 8.0).abs() < 1e-15);
    assert!((field.get(2, 1, 2, 3) - 9.0).abs() < 1e-15);
}

#[test]
fn vector_field_scale() {
    let grid = Grid::new(4, 100.0);
    let mut field = VectorField::zeros(&grid);

    field.component_mut(0)[[0, 0, 0]] = 3.0;
    field.component_mut(1)[[0, 0, 0]] = 4.0;

    field.scale(2.0);

    assert!((field.get(0, 0, 0, 0) - 6.0).abs() < 1e-15);
    assert!((field.get(1, 0, 0, 0) - 8.0).abs() < 1e-15);
}

#[test]
fn vector_field_add() {
    let grid = Grid::new(4, 100.0);
    let mut a = VectorField::zeros(&grid);
    let mut b = VectorField::zeros(&grid);

    a.component_mut(0)[[0, 0, 0]] = 1.0;
    b.component_mut(0)[[0, 0, 0]] = 2.0;

    let c = &a + &b;
    assert!((c.get(0, 0, 0, 0) - 3.0).abs() < 1e-15);
}
