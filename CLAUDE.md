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

`hermes-rs` implements **Hierarchical Closure Dynamics (HCD)** -- a multi-scale cosmological simulation framework. Two dynamical modes at scale-0:

- **Particle-mesh (PM):** Dark matter N-body with CIC, FFT-Poisson gravity, symplectic KDK leapfrog
- **Schrodinger-Poisson (SP):** Fuzzy dark matter wavefunction via split-step spectral integrator

The Content abstraction (Particles / Fields / Mixed) with pluggable Dynamics lets both modes coexist. Depends on `morphis` (geometric algebra) for all physical quantities.

### Boundary Between morphis and hermes

- **morphis knows**: elements, products, linear maps, decompositions, fields, spectral operators
- **hermes knows**: grids, time integration, scale hierarchy, closure terms, conservation monitoring, I/O, visualization

## Module Layout

```
src/
  algebra.rs            shared Euclidean 3-metric, morphis conversions
  colormap.rs           density/velocity colormapping
  error.rs              HermesError enum
  config/               TOML configuration (three-tier: ontology/simulation/output)
    ontology.rs         Ontology, Spacetime, ParticleSpecies, FieldSpecies, Lagrangian
    simulation.rs       GridConfig, TimeConfig, InitializationConfig
    output.rs           SnapshotsConfig, DiagnosticsConfig, LoggingConfig, DisplayConfig
    legacy.rs           old Configuration type (migration period)
  engine/               composable physics engine
    state.rs            SimulationState: named maps of particles and fields
    free/               FreeEvolution trait + per-species modules
      schrodinger.rs    SchrodingerEvolution (split-step kinetic propagator)
    coupling/           Coupling trait + cross-species modules
      poisson.rs        PoissonGravity (shared Poisson solver for all species)
  core/                 legacy simulation orchestration
    content.rs          Content enum (Particles / Fields / Mixed)
    dynamics.rs         Dynamics trait (pluggable step function)
    pm_dynamics.rs      ParticleMeshDynamics (delegates to PoissonGravity)
    schrodinger_dynamics.rs  SchrodingerPoissonDynamics (delegates to PoissonGravity)
    simulation.rs       simulation driver (from_scene, from_config, resume)
  physics/              physical models and solvers
    constants.rs        physical constants (kpc, M_sun, Gyr, eV)
    cosmology.rs        FLRW background, growth factor, kick/drift factors
    grid.rs             periodic cubic grid geometry
    field.rs            grade-0 and grade-1 fields with morphis extraction
    particles.rs        SoA particle storage with morphis interface
    cic.rs              cloud-in-cell mass assignment + force interpolation
    poisson.rs          FFT Poisson solver (ndrustfft)
    integrator.rs       symplectic KDK leapfrog
    diagnostics.rs      conservation audits
    spectral.rs         FFT wrappers (fft_3d, ifft_3d)
  io/
    snapshot.rs         Snapshot type, SnapshotContent enum, bincode serialization
    observer.rs         Observer trait (legacy, used by tests)
  run/
    cli.rs              clap-based CLI (--scene, --live, --playback, --resume)
    runner.rs           mode routing (headless, live, playback, record, resume)
    pipeline.rs         threaded pipeline: router, disk writer, precompute, viewer
  physics/initial/      initialization dispatched by config method string
    zeldovich.rs        Zel'dovich particle init from CDM power spectrum
    zeldovich_field.rs  Zel'dovich + random wavefunction init
    nfw.rs              NFW halo particle sampling
    nfw_field.rs        NFW halo plane-wave field init
  scenes/               legacy Scene trait (fallback, being phased out)
  visuals/              (#[cfg(feature = "vis")])
    viewer.rs           static 3D particle viewer (kiss3d)
    plots.rs            density slices, P(k), conservation plots (plotters)
    volumetric_renderer.rs  additive-blended Gaussian point sprites for fields
```

## Variable Naming

Use **descriptive names** based on the physical quantity, not mathematical symbols. This is applied-math/simulation code, not pure-algebra code -- prefer words over single letters (contrast with morphis-rs, which uses single-letter mathematical style).

**Pattern:** `<quantity>_<subscript>_<superscript>`

**Subscript/label attachment rule:** If the base name is a word, use an underscore: `mass_h`, `hubble_z`. If the base is a single letter or symbol, attach directly: `H0`, `Hz`.

| Math Symbol | Variable Name | Description |
|-------------|---------------|-------------|
| $H_0$ | `H0` | Hubble constant |
| $Ω_m$ | `omega_m` | Matter density parameter |
| $Ω_v$ | `omega_v` | Vacuum energy density parameter |
| $\bar{ρ}_m$ | `density_mean` | Mean comoving matter density |
| $m_p$ | `mass_particle` | Particle mass |
| $D_+(a)$ | `growth_factor` | Linear growth factor |
| $f(a)$ | `growth_rate` | Logarithmic growth rate |
| $ℓ/m$ | `length_scale` | Field smoothing length ratio |

In documentation: use α (not ψ) for the dark matter wavefunction field.

## LaTeX in Docs

GitHub's markdown renderer does not support `\,` for thin spaces. Use ` \ ` (backslash-space) instead in all LaTeX within this project's markdown files.

## Units

Internal units: kpc, M_☉, Gyr, eV ($k_B = 1$). Matches the plexis sibling project. Constants in `physics/constants.rs` are ported from `plexis/core/constants.py`.

## File Naming

Non-code files (scenes, snapshots, output) use hyphens: `snapshot-00000.bin`, `cosmic-web-pm`. Rust source files use underscores per Rust convention. Scene naming convention: `-pm` for particle-mesh, `-ft` for field theory.

## Pipeline Architecture

Simulation, disk I/O, and visualization run on independent threads connected by bounded channels. The simulation always runs on a spawned thread; main owns the event loop. `Arc<Snapshot>` enables zero-copy fan-out from the router to multiple consumers. Disk writer is non-droppable (every snapshot saved); viewer is droppable (frames silently dropped under load).

## Configuration

TOML-based with deep merge. Scene configs live in `scenes/*.toml` and define the full simulation. Base defaults are compiled in from `src/config/base.toml`.

**EngineConfig:** Three top-level containers: `[ontology]` (spacetime, particles, fields, lagrangian), `[simulation]` (grid, time, initialization), `[output]` (snapshots, diagnostics, logging, display). The CLI `--scene` flag points to any TOML file; the engine dispatches initialization by `[simulation.initialization.method]`.

**Legacy (Configuration):** Flat sections still used by the simulation driver internals. Built automatically from EngineConfig at runtime.

## Testing

All tests are integration tests in `tests/`. No `#[cfg(test)]` unit tests in `src/`. Tests verify conservation properties and algebraic laws using tolerance `1e-12`.

## CI

- Pre-commit hooks: `cargo fmt --check`, `cargo clippy -- -D warnings`, `cargo test`
- GitHub Actions: lint job (fmt + clippy) and test job, on push/PR to main
