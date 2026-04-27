# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
make lint              # cargo fmt + cargo clippy --fix
make test              # cargo test
make build             # cargo build --release
make clean             # cargo clean + remove .DS_Store
make reset             # clean + build
make publish           # git tag from Cargo.toml version + push
cargo test <test_name> # run a single test
```

## Project Overview

`hermes-rs` implements **Hierarchical Closure Dynamics (HCD)** -- a multi-scale simulation framework where coarse-grained state spaces are connected bidirectionally to fine-grained patches via restriction and prolongation operators. The closure problem (coarse equations needing unresolved fine-scale correlations) is addressed by learned rate-matrix corrections that preserve conservation by construction.

Depends on `morphis` (`../morphis-rs`) for geometric algebra: vectors, multivectors, metrics, products, versors, outermorphisms.

### Core Concepts

- **Grade-stratified state spaces**: Each scale $s$ has state space $\mathcal{V}_s = \oplus_{k=0}^{k_{\max}(s)} \mathcal{G}^k \otimes \mathcal{f}_s$. Grade activation $k_{\max}(s)$ increases with resolution -- cosmological scales carry only scalars; galactic scales activate vectors and bivectors; stellar scales use the full algebra.

- **Master-equation dynamics**: $\dot{\Psi}^{(s)} = [G^{(s)} + \delta Q^{(s)}]\Psi^{(s)} + S^{(s)} + A^{(s)}$ where $G$ is the resolvable generator, $\delta Q$ is the learned closure (rate matrix with zero column sums), $S$ is sources, $A$ is accretion/boundary.

- **Restriction** $R^{(s \leftarrow s+1)}$: spatial averaging + grade projection (fine → coarse). Lossy, deterministic.

- **Prolongation** $P^{(s+1 \leftarrow s)} = P_{\text{refine}} + P_{\text{expand}}$: refinement is deterministic interpolation (grade-preserving); expansion populates newly-activated grades (generative). Compatibility: $R \circ P_{\text{refine}} = \mathbb{1}$, $R \circ P_{\text{expand}} = 0$.

- **Zoom lifecycle**: trigger → prolong coarse to fine → integrate fine dynamics → restrict back → feed residual as closure → discard fine patch.

- **Conservation**: $K_\alpha^{(s)} = K_\alpha^{(s+1)} \circ R^{(s \leftarrow s+1)}$. Rate-matrix closure ensures mass/energy/momentum conservation structurally (zero column sums), not via penalty.

- **Bivector-native physics**: magnetic field as bivector (oriented plane), angular momentum as bivector. Not pseudovectors.

### Boundary Between morphis and hermes

- **morphis knows**: elements, products, linear maps, decompositions
- **hermes knows**: grids, time integration, scale hierarchy, restriction/prolongation operators, closure terms, conservation monitoring, zoom triggers

## Variable Naming

Use **descriptive names** based on the physical quantity, not mathematical symbols. The root should be the **kind of quantity**, with subscripts/superscripts appended as suffixes. This is applied-math/simulation code, not pure-algebra code -- prefer words over single letters (contrast with morphis-rs, which uses single-letter mathematical style).

**Pattern:** `<quantity>_<subscript>_<superscript>`

**Subscript/label attachment rule:** If the base name is a word, use an underscore: `mass_h`, `hubble_z`. If the base is a single letter or symbol, attach directly: `H0`, `Hz`.

| Math Symbol | Variable Name | Description |
|-------------|---------------|-------------|
| $H_0$ | `H0` | Hubble constant |
| $H(a)$ | `hubble_a` | Hubble parameter at scale factor $a$ |
| $\Omega_m$ | `omega_m` | Matter density parameter |
| $\bar{\rho}_m$ | `density_mean` | Mean comoving matter density |
| $m_p$ | `mass_particle` | Particle mass |
| $\mathbf{p}_i$ | `momentum` | Canonical momentum |
| $D_+(a)$ | `growth_factor` | Linear growth factor |
| $f(a)$ | `growth_rate` | Logarithmic growth rate |
| $c_s$ | `speed_sound` | Sound speed |
| $\rho_b$ | `density_baryon` | Baryon density |
| $\mathbf{B}$ | `field_magnetic` | Magnetic bivector field |

For science-coding conventions (config structure, naming patterns, physical modeling style), consult the sibling project `../plexis/` and its CLAUDE.md.

## Units

Internal units: kpc, M_☉, Gyr, eV (k_B = 1). Matches the plexis sibling project. Constants in `src/constants.rs` are ported from `plexis/core/constants.py`.

## Module Layout

- `src/constants.rs` -- physical constants in working units
- `src/cosmology.rs` -- `Cosmology` struct, Friedmann equation, growth factor, kick/drift factors
- `src/config.rs` -- `Configuration`, TOML loading with serde, three-tier deep merge
- `src/error.rs` -- `HermesError` enum

## Configuration

TOML-based, three-tier: embedded `configs/defaults.toml` → optional file override → programmatic overrides. Partial files are deep-merged. Load with `config::load_defaults()` or `config::build_configuration()`.

## Dependencies

- `morphis` -- geometric algebra (Vector, MultiVector, Metric, products, versors, outermorphisms)
- `serde` + `toml` -- configuration
- `thiserror` -- error types

Dev: `proptest`, `approx`

## Testing

All tests are integration tests in `tests/`. No `#[cfg(test)]` unit tests in `src/`. One file per module being tested. Tests verify algebraic laws and conservation properties, using tolerance `1e-12` for floating-point comparisons.

## CI

- Pre-commit hooks: `cargo fmt --check`, `cargo clippy -- -D warnings`, `cargo test`
- GitHub Actions (`.github/workflows/ci.yml`): lint job (fmt + clippy) and test job, on push/PR to main
