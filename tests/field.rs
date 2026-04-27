use hermes_rs::algebra::{euclidean_3, vector_from_components};
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

// ============================================================================
// Morphis geometric algebra integration
// ============================================================================

#[test]
fn scalar_field_carries_euclidean_metric() {
    let grid = Grid::new(4, 100.0);
    let field = ScalarField::zeros(&grid);

    assert_eq!(field.metric, euclidean_3());
}

#[test]
fn scalar_field_scalar_at() {
    let grid = Grid::new(4, 100.0);
    let mut field = ScalarField::zeros(&grid);
    *field.get_mut(1, 2, 3) = 5.0;

    let s = field.scalar_at(1, 2, 3);

    assert_eq!(s.grade(), 0);
    assert!((s.component(&[]) - 5.0).abs() < 1e-15);
}

#[test]
fn vector_field_carries_euclidean_metric() {
    let grid = Grid::new(4, 100.0);
    let field = VectorField::zeros(&grid);

    assert_eq!(field.metric, euclidean_3());
}

#[test]
fn vector_field_vector_at() {
    let grid = Grid::new(4, 100.0);
    let mut field = VectorField::zeros(&grid);
    field.component_mut(0)[[1, 2, 3]] = 7.0;
    field.component_mut(1)[[1, 2, 3]] = 8.0;
    field.component_mut(2)[[1, 2, 3]] = 9.0;

    let v = field.vector_at(1, 2, 3);

    assert_eq!(v.grade(), 1);
    assert!((v.component(&[0]) - 7.0).abs() < 1e-15);
    assert!((v.component(&[1]) - 8.0).abs() < 1e-15);
    assert!((v.component(&[2]) - 9.0).abs() < 1e-15);
}

#[test]
fn vector_field_set_vector_at() {
    let grid = Grid::new(4, 100.0);
    let mut field = VectorField::zeros(&grid);

    let v = vector_from_components(1.0, 2.0, 3.0);
    field.set_vector_at(1, 1, 1, &v);

    let got = field.vector_at(1, 1, 1);
    assert!((got.component(&[0]) - 1.0).abs() < 1e-15);
    assert!((got.component(&[1]) - 2.0).abs() < 1e-15);
    assert!((got.component(&[2]) - 3.0).abs() < 1e-15);
}

#[test]
fn vector_field_vector_at_norm() {
    let grid = Grid::new(4, 100.0);
    let mut field = VectorField::zeros(&grid);
    field.component_mut(0)[[0, 0, 0]] = 3.0;
    field.component_mut(1)[[0, 0, 0]] = 4.0;

    let v = field.vector_at(0, 0, 0);

    assert!((v.norm() - 5.0).abs() < 1e-12);
}

#[test]
fn vector_field_dot_at() {
    let grid = Grid::new(4, 100.0);
    let mut field_1 = VectorField::zeros(&grid);
    let mut field_2 = VectorField::zeros(&grid);

    field_1.component_mut(0)[[0, 0, 0]] = 1.0;
    field_1.component_mut(1)[[0, 0, 0]] = 2.0;
    field_1.component_mut(2)[[0, 0, 0]] = 3.0;

    field_2.component_mut(0)[[0, 0, 0]] = 4.0;
    field_2.component_mut(1)[[0, 0, 0]] = 5.0;
    field_2.component_mut(2)[[0, 0, 0]] = 6.0;

    let dot = field_1.dot_at(0, 0, 0, &field_2);
    let expected = 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0;

    assert!(
        (dot - expected).abs() < 1e-12,
        "dot_at: expected {expected}, got {dot}"
    );
}

#[test]
fn scalar_field_add_preserves_metric() {
    let grid = Grid::new(4, 100.0);
    let a = ScalarField::zeros(&grid);
    let b = ScalarField::zeros(&grid);
    let c = &a + &b;

    assert_eq!(c.metric, euclidean_3());
}
