# Particle-Mesh Method

The particle-mesh (PM) method computes gravitational forces by depositing particle masses onto a grid, solving the Poisson equation in Fourier space, and interpolating the resulting force field back to particle positions. It is the simplest method that does real cosmological physics and the natural foundation for the HCD framework's scale-0 dynamics.

## The Force Chain

The gravitational generator $G^{(0)}$ is built from three operators composed as a pipeline:

$$
\{x_n, p_n\} \xrightarrow{\mathcal{M}} ρ \xrightarrow{\mathcal{P}} F(x_g) \xrightarrow{\mathcal{I}} \{F_n\}
$$

### Mass Assignment ($\mathcal{M}$)

Cloud-in-cell (CIC) interpolation deposits each particle's mass across the 8 surrounding grid cells using trilinear weights. For particle $n$ at position $x_n$ with cell spacing $h$:

$$
ρ(x_g) = \frac{1}{h^3} \sum_n m_p \, W\left(\frac{x_n - x_g}{h}\right)
$$

where $W$ is the product of one-dimensional triangle functions along each axis.

**Implementation:** `physics::cic::assign_density` returns a `ScalarField` carrying the Euclidean 3-metric. The type signature is grade-1 particle list to grade-0 density field.

### Poisson Solve ($\mathcal{P}$)

The overdensity $δ = ρ / \bar{ρ} - 1$ is transformed to Fourier space via a 3D real-to-complex FFT (ndrustfft). The gravitational potential is obtained by multiplication with the discrete Green's function:

$$
\hat{ϕ}(k) = \frac{4π G \bar{ρ} a^2}{k^2_{\text{discrete}}} \hat{δ}(k)
$$

where $k^2_{\text{discrete}} = (2/h)^2 (\sin^2 k_x h/2 + \sin^2 k_y h/2 + \sin^2 k_z h/2)$ is the finite-difference Laplacian in Fourier space. The zero mode is set to zero (periodic box has no defined absolute potential). The force field is obtained by multiplication with $-ik_d$ for each component, then inverse FFT.

**Implementation:** `physics::poisson::PoissonSolver` precomputes the Green's function at construction time and reuses FFT handlers across calls. The 3D R2C transform is composed from axis-wise 1D transforms: R2C along axis 2, then C2C along axes 1 and 0.

### Force Interpolation ($\mathcal{I}$)

The grid force field is interpolated back to particle positions using the same CIC kernel as mass assignment. Sharing the kernel between assignment and interpolation is the discrete analog of Newton's third law and ensures momentum conservation to machine precision.

**Implementation:** `physics::cic::interpolate_force` returns `ParticleForces` with morphis-native access via `force_on(n) -> Vector<3>`.

## State Vector

The persistent state at scale 0 is:

$$
Ψ^{(0)} = \big(\{x_n\}, \{p_n\}, a\big)
$$

with comoving positions $x_n \in \mathbb{R}^3 / L\mathbb{Z}^3$ (periodic) and canonical momenta $p_n = m_p a^2 \dot{x}_n$. The grid fields (density, potential, force) are derived per-step and not part of the persistent state.

In the framework's grade vocabulary, $g_{\max}(0) = 1$. Scalars (density, potential) and vectors (positions, momenta, forces) are dynamically active. The angular momentum bivector $L = x \wedge p$ is derivable as a diagnostic via the morphis wedge product but is not separately evolved.

**Implementation:** `physics::particles::Particles` stores positions and momenta as `Array2<f64>` with shape `[3, N_p]` for cache-friendly CIC access. The primary interface returns morphis `Vector<3>` objects: `position_of(n)`, `momentum_of(n)`, `angular_momentum(n)` (grade-2 bivector).

## Time Integration

A symplectic kick-drift-kick leapfrog advances the state with cosmological step factors:

1. **Half-kick:** $p \to p + F \times K(a_n, a_{n+1/2})$
2. **Full drift:** $x \to x + (p / m_p) \times D(a_n, a_{n+1})$
3. **Recompute force** at new positions
4. **Half-kick:** $p \to p + F \times K(a_{n+1/2}, a_{n+1})$

where the kick and drift factors are integrals over the Hubble parameter:

$$
K(a_0, a_1) = \int_{a_0}^{a_1} \frac{da}{a^2 H(a)}, \qquad D(a_0, a_1) = \int_{a_0}^{a_1} \frac{da}{a^3 H(a)}
$$

The closing half-kick force is cached for reuse as the next step's opening half-kick, reducing the per-step cost from two force evaluations to one.

**Implementation:** `physics::integrator` provides `kick`, `drift`, and `step_kdk`. Kicks and drifts are morphis vector operations: `p_new = &p + &(&force * kick_factor)`.

## Initialization

Zel'dovich approximation from a linear power spectrum. Particles start on a uniform Lagrangian lattice and are displaced:

$$
x_n = q_n + D_+(a_{\text{init}}) \, Ψ(q_n), \qquad p_n = m_p a^2 f H \, Ψ(q_n)
$$

where $Ψ(k) = ik / k^2 \, \hat{δ}(k)$ is the displacement field computed from a Gaussian random overdensity. The power spectrum uses the Eisenstein-Hu no-wiggle transfer function normalized by $σ_8$.

**Implementation:** `physics::initial::zeldovich_init` returns `Particles` with morphis-native positions and momenta.

## Conservation

Conservation laws stratify by grade:

| Grade | Quantity | Conservation |
|-------|----------|-------------|
| 0 | Total mass $M = N_p m_p$ | Exact (fixed particle count) |
| 1 | Comoving momentum $P = \sum p_n$ | Machine precision (CIC kernel symmetry) |
| 0 | Energy (Layzer-Irvine) | Approximate (timestep error) |
| 2 | Angular momentum $L = \sum x_n \wedge p_n$ | Not conserved (periodic box) |

**Implementation:** `physics::diagnostics::Diagnostics` computes all quantities through morphis operations: `total_momentum()` via vector sum, `angular_momentum()` via wedge product, `kinetic_energy()` via `norm_squared()`.
