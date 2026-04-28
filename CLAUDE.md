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

`hermes-rs` implements **Hierarchical Closure Dynamics (HCD)** -- a multi-scale cosmological simulation framework. The current implementation is a scale-0 particle-mesh N-body simulator with FFT gravity, Zel'dovich initialization, and a threaded pipeline architecture for I/O and visualization.

Depends on `morphis` (geometric algebra) for all physical quantities: positions, momenta, and forces are morphis grade-1 vectors; angular momentum is a grade-2 bivector via the wedge product.

### Boundary Between morphis and hermes

- **morphis knows**: elements, products, linear maps, decompositions
- **hermes knows**: grids, time integration, scale hierarchy, closure terms, conservation monitoring, I/O, visualization

## Module Layout

```
src/
  algebra.rs          shared Euclidean 3-metric, morphis conversions
  colormap.rs         density/velocity colormapping
  config.rs           TOML configuration with deep merge
  error.rs            HermesError enum
  physics/            simulation core
    constants.rs      physical constants (kpc, M_☉, Gyr, eV)
    cosmology.rs      FLRW background, growth factor, kick/drift factors
    grid.rs           periodic cubic grid
    field.rs          grade-0 and grade-1 fields with morphis extraction
    particles.rs      SoA particle storage with morphis interface
    cic.rs            cloud-in-cell mass assignment + force interpolation
    poisson.rs        FFT Poisson solver (ndrustfft)
    integrator.rs     symplectic KDK leapfrog
    diagnostics.rs    conservation audits
    simulation.rs     simulation driver (from_scene + from_config)
  io/                 data I/O
    snapshot.rs       Snapshot type, bincode serialization
    observer.rs       Observer trait (legacy, used by tests)
  run/                execution
    cli.rs            clap-based CLI
    runner.rs         mode routing (headless, live, playback, record)
    pipeline.rs       threaded pipeline: router, disk writer, precompute, viewer
  scenes/             simulation scenarios (each a subdirectory)
    cosmic_web/       Zel'dovich PM in a 100 Mpc periodic box (default)
    galaxy_group/     constrained Zel'dovich in a 3 Mpc box
  vis/                visualization (#[cfg(feature = "vis")])
    viewer.rs         static 3D particle viewer (kiss3d)
    plots.rs          density slices, P(k), conservation plots (plotters)
```

## Variable Naming

Use **descriptive names** based on the physical quantity, not mathematical symbols. This is applied-math/simulation code, not pure-algebra code -- prefer words over single letters (contrast with morphis-rs, which uses single-letter mathematical style).

**Pattern:** `<quantity>_<subscript>_<superscript>`

**Subscript/label attachment rule:** If the base name is a word, use an underscore: `mass_h`, `hubble_z`. If the base is a single letter or symbol, attach directly: `H0`, `Hz`.

| Math Symbol | Variable Name | Description |
|-------------|---------------|-------------|
| $H_0$ | `H0` | Hubble constant |
| $Ω_m$ | `omega_m` | Matter density parameter |
| $\bar{ρ}_m$ | `density_mean` | Mean comoving matter density |
| $m_p$ | `mass_particle` | Particle mass |
| $D_+(a)$ | `growth_factor` | Linear growth factor |
| $f(a)$ | `growth_rate` | Logarithmic growth rate |
| $ρ_b$ | `density_baryon` | Baryon density |

For science-coding conventions (config structure, naming patterns, physical modeling style), consult the sibling project `../plexis/` and its CLAUDE.md.

## Units

Internal units: kpc, M_☉, Gyr, eV ($k_B = 1$). Matches the plexis sibling project. Constants in `physics/constants.rs` are ported from `plexis/core/constants.py`.

## File Naming

Non-code files (configs, snapshots, output) use hyphens: `snapshot-00000.bin`, `cosmic-web`. Timestamp directories use `<date>_<time>` format. Rust source files use underscores per Rust convention.

## Pipeline Architecture

Simulation, disk I/O, and visualization run on independent threads connected by bounded channels. The simulation always runs on a spawned thread; main owns the event loop. `Arc<Snapshot>` enables zero-copy fan-out from the router to multiple consumers.

## Configuration

TOML-based, three-tier: embedded `configs/defaults.toml` → optional file override → CLI overrides. Partial files are deep-merged.

## Testing

All tests are integration tests in `tests/`. No `#[cfg(test)]` unit tests in `src/`. Tests verify conservation properties and algebraic laws using tolerance `1e-12`.

## CI

- Pre-commit hooks: `cargo fmt --check`, `cargo clippy -- -D warnings`, `cargo test`
- GitHub Actions: lint job (fmt + clippy) and test job, on push/PR to main
