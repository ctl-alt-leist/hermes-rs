use hermes_rs::algebra::{
    components_from_vector, euclidean_3, scalar_from_f64, vector_from_array, vector_from_components,
};
use morphis::metric::Signature;
use morphis::ops::{geometric, interior_left, wedge};
use morphis::vector::basis;

// ============================================================================
// Metric
// ============================================================================

#[test]
fn euclidean_3_is_euclidean() {
    let g = euclidean_3();
    assert_eq!(g.sig, Signature::Euclidean);
    assert_eq!(g.diag, [1.0, 1.0, 1.0]);
}

// ============================================================================
// Vector construction roundtrips
// ============================================================================

#[test]
fn vector_roundtrip() {
    let v = vector_from_components(1.0, 2.0, 3.0);
    let c = components_from_vector(&v);

    assert!((c[0] - 1.0).abs() < 1e-15);
    assert!((c[1] - 2.0).abs() < 1e-15);
    assert!((c[2] - 3.0).abs() < 1e-15);
}

#[test]
fn vector_from_array_roundtrip() {
    let v = vector_from_array([4.0, 5.0, 6.0]);

    assert_eq!(v.grade(), 1);
    assert!((v.component(&[1]) - 4.0).abs() < 1e-15);
    assert!((v.component(&[2]) - 5.0).abs() < 1e-15);
    assert!((v.component(&[3]) - 6.0).abs() < 1e-15);
}

#[test]
fn scalar_roundtrip() {
    let s = scalar_from_f64(42.0);

    assert_eq!(s.grade(), 0);
    assert!((s.component(&[]) - 42.0).abs() < 1e-15);
}

// ============================================================================
// Morphis operations — norms, arithmetic, products
// ============================================================================

#[test]
fn vector_norm_euclidean() {
    let v = vector_from_components(3.0, 4.0, 0.0);

    assert!((v.norm() - 5.0).abs() < 1e-12);
}

#[test]
fn vector_norm_squared() {
    let v = vector_from_components(1.0, 2.0, 3.0);

    assert!((v.norm_squared() - 14.0).abs() < 1e-12);
}

#[test]
fn vector_addition() {
    let u = vector_from_components(1.0, 0.0, 0.0);
    let v = vector_from_components(0.0, 2.0, 3.0);
    let w = &u + &v;

    assert!((w.component(&[1]) - 1.0).abs() < 1e-15);
    assert!((w.component(&[2]) - 2.0).abs() < 1e-15);
    assert!((w.component(&[3]) - 3.0).abs() < 1e-15);
}

#[test]
fn vector_scalar_multiplication() {
    let v = vector_from_components(1.0, 2.0, 3.0);
    let scaled = &v * 2.0;

    assert!((scaled.component(&[1]) - 2.0).abs() < 1e-15);
    assert!((scaled.component(&[2]) - 4.0).abs() < 1e-15);
    assert!((scaled.component(&[3]) - 6.0).abs() < 1e-15);
}

#[test]
fn dot_product_via_geometric() {
    let u = vector_from_components(1.0, 2.0, 3.0);
    let v = vector_from_components(4.0, 5.0, 6.0);

    // For grade-1 vectors, u * v = u·v + u∧v.
    // The scalar part is the dot product.
    let product = geometric(&u, &v);
    let dot = product.scalar_part();
    let expected = 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0;

    assert!(
        (dot - expected).abs() < 1e-12,
        "dot product: expected {expected}, got {dot}"
    );
}

#[test]
fn wedge_product_produces_bivector() {
    let e = basis(euclidean_3());
    let bivector = wedge(&e[1], &e[2]);

    assert_eq!(bivector.grade(), 2);
    assert!((bivector.component(&[1, 2]) - 1.0).abs() < 1e-15);
    assert!((bivector.component(&[2, 1]) + 1.0).abs() < 1e-15);
}

#[test]
fn wedge_product_antisymmetric() {
    let u = vector_from_components(1.0, 2.0, 0.0);
    let v = vector_from_components(3.0, 4.0, 0.0);

    let uv = wedge(&u, &v);
    let vu = wedge(&v, &u);

    // u ∧ v = -(v ∧ u)
    let sum = &uv + &vu;
    assert!(sum.is_zero(1e-12), "wedge product should be antisymmetric");
}

#[test]
fn interior_product_contracts_grade() {
    let e = basis(euclidean_3());
    let bivector = wedge(&e[1], &e[2]); // grade 2

    // e1 ⌋ (e1 ∧ e2) should give e2 (grade 1)
    let result = interior_left(&e[1], &bivector);

    assert_eq!(result.grade(), 1);
    assert!((result.component(&[2]) - 1.0).abs() < 1e-12);
}
