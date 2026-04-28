# Morphis Field Abstraction — Design Brief

This document specifies a field abstraction for morphis-rs. A field is a spatially-varying geometric algebra object on a periodic grid — a scalar field, a vector field, a bivector field, etc. The abstraction should feel like a natural extension of morphis's existing `Vector<D>`: the same metric awareness, the same grade structure, the same algebraic operations, but distributed over space.

## Guiding Principle

morphis already has the pointwise algebra: `Vector<D>` for k-vectors, `MultiVector<D>` for mixed grades, and the full product suite (wedge, geometric, interior). A field is a collection of these objects indexed by spatial position. The field abstraction should:

- Follow morphis's existing patterns, not invent new ones
- Be general to any field theory on a periodic domain, not cosmology-specific
- Treat the grid as part of the field's identity (like the metric is for vectors)
- Make extraction of pointwise morphis objects natural: `field.at(x) -> Vector<D>`
- Make spatial derivatives first-class operations that respect grade

## Core Type

```rust
/// A field of grade-k geometric algebra elements on a periodic grid.
///
/// Each grid point holds a `Vector<D>` of the specified grade.
/// The field carries the metric and grid geometry, enabling both
/// pointwise algebraic operations and spatial differential operators.
pub struct Field<const D: usize> {
    /// Field data. For grade k on an N^D grid, shape is [N, N, ..., N, D, D, ..., D]
    /// where the first D entries are spatial indices and the remaining k entries
    /// are the tensor indices of each element.
    pub data: ArrayD<f64>,
    /// Grade of each element (0 = scalar, 1 = vector, 2 = bivector, ...).
    grade: usize,
    /// Metric defining the inner product structure.
    pub metric: Metric<D>,
    /// Grid geometry (cell count per side, box length, cell length).
    pub grid: Grid<D>,
}
```

The `Grid<D>` type should also live in morphis (it's just periodic box geometry, not physics-specific):

```rust
/// Periodic grid geometry in D dimensions.
pub struct Grid<const D: usize> {
    /// Number of cells per side.
    pub n_cells: usize,
    /// Box side length.
    pub box_length: f64,
    /// Cell side length (derived: box_length / n_cells).
    pub cell_length: f64,
}
```

## Constructors

Following the `Vector<D>` pattern — free functions and associated methods:

```rust
/// Zero field: every grid point holds a zero k-vector.
pub fn zeros(grade: usize, grid: &Grid<D>, metric: Metric<D>) -> Field<D>

/// Constant field: every grid point holds the same value.
pub fn constant(value: &Vector<D>, grid: &Grid<D>) -> Field<D>

/// From a function: evaluate f(x) at each grid point.
pub fn from_fn(grade: usize, grid: &Grid<D>, metric: Metric<D>,
               f: impl Fn(&[f64; D]) -> Vector<D>) -> Field<D>

/// Scalar field from a scalar function.
pub fn scalar_field(grid: &Grid<D>, metric: Metric<D>,
                    f: impl Fn(&[f64; D]) -> f64) -> Field<D>
```

## Pointwise Access

The primary interface for extracting and inserting morphis objects:

```rust
/// Extract the value at grid point (m0, m1, m2, ...) as a morphis Vector<D>.
pub fn at(&self, indices: &[usize]) -> Vector<D>

/// Set the value at a grid point from a morphis Vector<D>.
pub fn set(&mut self, indices: &[usize], value: &Vector<D>)

/// Grade of elements in this field.
pub fn grade(&self) -> usize

/// Grid geometry.
pub fn grid(&self) -> &Grid<D>

/// Number of grid points (total).
pub fn n_points(&self) -> usize
```

## Pointwise Algebraic Operations

These lift Vector<D> operations to fields. Every operation that makes sense pointwise on vectors should work on fields:

```rust
// Field + Field (same grade, same grid)
impl Add for &Field<D> -> Field<D>
impl Sub for &Field<D> -> Field<D>
impl Neg for &Field<D> -> Field<D>

// Field * scalar
impl Mul<f64> for &Field<D> -> Field<D>

// Pointwise norm squared: returns a scalar field (grade 0)
pub fn norm_squared(&self) -> Field<D>

// Pointwise reverse
pub fn rev(&self) -> Field<D>
```

Pointwise products between fields of different grades:

```rust
/// Pointwise wedge product of two fields.
/// Grade of result = grade(self) + grade(other).
pub fn wedge(f: &Field<D>, g: &Field<D>) -> Field<D>

/// Pointwise left interior product.
pub fn interior_left(f: &Field<D>, g: &Field<D>) -> Field<D>

/// Pointwise scalar product: returns a scalar field.
pub fn scalar_product(f: &Field<D>, g: &Field<D>) -> Field<D>
```

## Spatial Derivatives

These are the field-theoretic operations that distinguish a field from a collection of vectors. All derivatives on a periodic domain are computed spectrally via FFT.

### Gradient (grade-raising)

The gradient of a grade-k field is a grade-(k+1) field:

```text
grad(f) = ∇f = Σ_a e_a ∧ ∂_a f
```

For a scalar field (grade 0), the gradient is a vector field (grade 1).
For a vector field (grade 1), the gradient includes the curl (grade 2).

```rust
/// Gradient: raises grade by 1 via exterior derivative.
/// Scalar field → vector field, vector field → bivector field, etc.
pub fn grad(&self) -> Field<D>
```

### Divergence (grade-lowering)

The divergence of a grade-k field is a grade-(k-1) field:

```text
div(f) = ∇ · f = Σ_a e_a ⌋ ∂_a f
```

For a vector field (grade 1), the divergence is a scalar field (grade 0).

```rust
/// Divergence: lowers grade by 1 via interior derivative.
/// Vector field → scalar field, bivector field → vector field, etc.
pub fn div(&self) -> Field<D>
```

### Curl (grade-preserving, D=3 only)

In 3D, the curl of a vector field is a vector field (via the Hodge dual of the exterior derivative):

```text
curl(f) = ∇ ∧ f  (for grade-1 fields in 3D, returns grade-2)
```

More precisely, this is the exterior derivative applied to a 1-form, yielding a 2-form. In 3D with the Hodge dual, this maps back to a vector. But in geometric algebra, we should keep the natural grade — the curl of a vector field is a bivector field. The user can Hodge-dualize if they want a vector.

```rust
/// Exterior derivative: ∇ ∧ f.
/// For a vector field in 3D, this is the curl (returns bivector field).
pub fn curl(&self) -> Field<D>
```

### Laplacian (grade-preserving)

The scalar Laplacian acts on each component:

```text
∇²f = Σ_a ∂²_a f
```

```rust
/// Laplacian: grade-preserving second derivative.
pub fn laplacian(&self) -> Field<D>
```

### Implementation: Spectral Derivatives

All derivatives are computed in Fourier space:

- ∂_a f → multiply f_hat(k) by i*k_a
- ∇²f → multiply f_hat(k) by -|k|²
- The gradient, divergence, and curl are built from per-component partial derivatives combined with the appropriate algebraic operation (wedge for grad/curl, interior for div)

The field should support forward and inverse FFT:

```rust
/// Forward FFT: real-space field → Fourier-space representation.
pub fn fft(&self) -> FieldFourier<D>

/// Inverse FFT: Fourier-space → real-space field.
pub fn from_fourier(fourier: &FieldFourier<D>, grid: &Grid<D>, metric: Metric<D>) -> Field<D>
```

The `FieldFourier<D>` type holds the complex Fourier coefficients. Derivative operations can work either by transforming to Fourier space, applying the operator, and transforming back, or by providing a Fourier-space API for chaining multiple operations before a single inverse transform.

## Integration (Scalar Fields)

```rust
/// Volume integral of a scalar field: ∫ f dV.
/// Returns f64 for grade-0 fields.
pub fn integrate(&self) -> f64

/// Sum of all values (no volume weighting).
pub fn sum(&self) -> f64
```

## Interplay with Existing morphis Types

The field abstraction should compose naturally with existing morphis operations:

```rust
let g = euclidean::<3>();
let grid = Grid::<3>::new(64, 100.0);

// Scalar field
let rho = Field::scalar_field(&grid, g, |x| {
    (2.0 * PI * x[0] / 100.0).sin()
});

// Extract a morphis scalar at a point
let rho_center: Vector<3> = rho.at(&[32, 32, 32]);
assert_eq!(rho_center.grade(), 0);

// Gradient: scalar → vector field
let grad_rho = rho.grad();
assert_eq!(grad_rho.grade(), 1);

// Extract a morphis vector at a point
let v: Vector<3> = grad_rho.at(&[0, 0, 0]);
assert_eq!(v.grade(), 1);

// Bivector field (e.g., magnetic field)
let b_field = Field::zeros(2, &grid, g);
assert_eq!(b_field.grade(), 2);

// Divergence of bivector → vector field
let div_b = b_field.div();
assert_eq!(div_b.grade(), 1);

// Pointwise norm squared: bivector field → scalar field
let b_energy = b_field.norm_squared();
assert_eq!(b_energy.grade(), 0);

// Laplacian: grade-preserving
let laplacian_rho = rho.laplacian();
assert_eq!(laplacian_rho.grade(), 0);
```

## Even-Subalgebra Fields (Wavefunctions)

In 3D, the even subalgebra G^+ = G^0 + G^3 is isomorphic to the complex numbers via the pseudoscalar I = e_1 ^ e_2 ^ e_3, which satisfies I^2 = -1. A wavefunction is a field valued in this even subalgebra — each grid point holds a scalar part plus a pseudoscalar part:

```text
α(x) = a(x) + b(x) I
```

This is a **mixed-grade field** (grades 0 and 3 simultaneously), not a pure-grade field. The field abstraction must support this. Options:

1. **MultiVectorField<D>** — a field of MultiVector<D> values. General but heavy.
2. **EvenField<D>** — a specialized field restricted to the even subalgebra. Stores two real arrays (scalar part and pseudoscalar coefficient). Lighter, and the even subalgebra closure under multiplication is a useful compile-time guarantee.

Recommended: **EvenField<D>** as a first-class type, since wavefunctions are the primary use case and the even subalgebra has special algebraic properties (closed under multiplication, isomorphic to C).

```rust
/// A field valued in the even subalgebra G^+ = G^0 ⊕ G^D.
///
/// In 3D, this is isomorphic to a complex-valued field: each point
/// holds a + bI where I is the unit pseudoscalar.
pub struct EvenField<const D: usize> {
    /// Scalar (grade-0) part.
    pub scalar: ArrayD<f64>,      // shape [N, N, ..., N]
    /// Pseudoscalar coefficient (grade-D part, without the I factor).
    pub pseudoscalar: ArrayD<f64>, // shape [N, N, ..., N]
    pub metric: Metric<D>,
    pub grid: Grid<D>,
}
```

Essential operations:

```rust
/// Reversal (complex conjugation): (a + bI) → (a - bI)
pub fn rev(&self) -> EvenField<D>

/// Pointwise product: (a + bI)(c + dI) = (ac - bd) + (ad + bc)I
/// Closed in the even subalgebra.
pub fn mul(&self, other: &EvenField<D>) -> EvenField<D>

/// Norm squared: α * α_rev = a² + b² (grade-0 scalar field)
pub fn norm_squared(&self) -> Field<D>

/// Phase rotation: multiply by exp(I θ) = cos(θ) + sin(θ) I
/// This is the kinetic-step operation in split-step integration.
pub fn rotate_phase(&self, angle: &Field<D>) -> EvenField<D>

/// Extract density: ρ = m * α_rev * α (scalar field)
pub fn density(&self, mass: f64) -> Field<D>

/// Extract velocity field: v = (1/m) ∇S where S is the phase
/// Requires gradient operation on the phase, so this depends on
/// the Field derivative infrastructure.
pub fn velocity(&self, mass: f64) -> Field<D>

/// Pointwise extraction as MultiVector<D>
pub fn at(&self, indices: &[usize]) -> MultiVector<D>
```

The `rotate_phase` method is crucial — it's what makes the split-step integrator work. In Fourier space, the kinetic step is a phase rotation by -hbar |k|^2 dt / (2m), applied pointwise in k-space.

## Conjugation and Reversal on Fields

Reversal is a fundamental operation for fields, not just vectors:

```rust
/// For a pure-grade field: multiply by (-1)^{k(k-1)/2} pointwise.
/// For an even-subalgebra field: flip the pseudoscalar part (complex conjugation).
pub fn rev(&self) -> Self
```

This is needed for:
- Computing density from wavefunctions: ρ = m * α_rev * α
- Computing probability currents: j = (I hbar / 2m)(α_rev ∇α - α ∇α_rev)
- Energy functionals: kinetic energy involves (∇α) ⌋ (∇α_rev)

## Laplacian Inverse (Poisson Solve)

The operation "given a source field f, find the field φ such that ∇²φ = f" is a generic spectral operation: multiply by -1/|k|^2 in Fourier space, with the zero mode set to zero. This is not physics-specific — any field theory with an elliptic constraint uses it.

```rust
/// Solve ∇²φ = f for φ on the periodic domain.
///
/// Spectral method: φ_hat(k) = -f_hat(k) / |k|² with φ_hat(0) = 0.
/// The zero mode is projected out (periodic domain has no unique
/// solution for the mean; this returns the zero-mean solution).
///
/// Grade-preserving: operates on each component independently.
pub fn laplacian_inverse(&self) -> Field<D>
```

This belongs in morphis because:
- It's a property of the Laplacian on periodic domains, not of gravity
- Any Helmholtz decomposition uses it
- The Green's function -1/|k|^2 is pure math, not physics
- hermes would then write `phi = (4 * pi * G * rho).laplacian_inverse()` — the physics is in the prefactor, the math is in the solve

## What Stays in hermes

The field abstraction in morphis is the mathematical substrate. The following are physics and remain in hermes:

- Physical constants and prefactors (4πG, hbar_eff, coupling constants)
- Split-step integrator orchestration (the sequence of kinetic/potential steps)
- Cosmological factors (scale factor, Hubble parameter, FLRW expansion)
- Initial condition generation (Zel'dovich, NFW profiles, wavefunction seeding)
- Visualization, I/O, pipeline

morphis provides: fields, grids, derivatives, pointwise algebra, even-subalgebra operations, Laplacian inverse. hermes provides: the physics that combines them into a simulation.

## Testing Strategy

Following morphis's existing test patterns (integration tests in `tests/`):

**Pure-grade fields:**
- **Gradient of sin(kx)** should give k*cos(kx) — spectral exactness
- **Divergence of gradient** should equal Laplacian — Hodge identity
- **Curl of gradient** should be zero — exact identity
- **Divergence of curl** should be zero — exact identity
- **Laplacian inverse of Laplacian** should roundtrip (up to zero mode)
- **Integration of constant** should give value * volume
- **Pointwise extraction roundtrip**: set a value, extract it, compare
- **Grade propagation**: grad raises, div lowers, laplacian preserves
- **Norm squared of basis vector field** should give constant scalar field of 1.0

**Even-subalgebra fields:**
- **Reversal is involution**: rev(rev(α)) = α
- **Norm squared is real**: α_rev * α has zero pseudoscalar part
- **Phase rotation preserves norm**: |exp(Iθ) α|^2 = |α|^2
- **Product closure**: EvenField * EvenField is EvenField (no odd-grade leakage)
- **Density extraction**: for α = sqrt(ρ/m) exp(I S/hbar), density(m) recovers ρ

**Laplacian inverse:**
- **Roundtrip**: laplacian(laplacian_inverse(f)) = f (for zero-mean f)
- **Known solution**: laplacian_inverse of sin(kx) gives -sin(kx)/k^2
- **Zero mode**: laplacian_inverse of a constant returns zero

## Notes

- The brief specifies D-dimensional fields for generality, but 3D (D=3) is the primary use case. 2D fields would be useful for testing and for 2D toy models.
- The `Grid<D>` type is kept simple (uniform, periodic, cubic). Adaptive grids, non-periodic boundaries, and non-cubic domains are future extensions that don't affect the field algebra.
- FFT implementation: morphis should depend on ndrustfft (already in hermes, pure Rust). The FFT infrastructure is general enough to belong in morphis alongside the field type.
- The EvenField<D> type is specific to odd-dimensional spaces where the pseudoscalar squares to -1. In even dimensions, I^2 = +1 and the even subalgebra splits differently. The 3D case is the priority.
