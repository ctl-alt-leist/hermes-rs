# Project Overview

Hermes is the dynamical simulation layer of the Hierarchical Closure Dynamics (HCD) framework. Its algebraic substrate is [morphis-rs](https://github.com/ctl-alt-leist/morphis-rs), a Rust geometric algebra library. morphis handles elements, products, and linear maps; hermes handles grids, time integration, scale hierarchies, and closure terms.

## Scope

The current implementation is a scale-0 particle-mesh (PM) dark-matter-only cosmological N-body simulation. This is the simplest instantiation of the generator-equation form

$$
\dot{Ψ}^{(s)} = G^{(s)} Ψ^{(s)} + S^{(s)} + C^{(s)} + A^{(s)}
$$

with the closure $C^{(0)} = 0$, sources $S^{(0)} = 0$, and the gravitational generator $G^{(0)}$ built from the FFT-Poisson force chain.

### What's Built

| Component | Module | Status |
|-----------|--------|--------|
| FLRW cosmological background | `physics::cosmology` | Complete |
| Periodic cubic grid | `physics::grid` | Complete |
| Scalar and vector fields | `physics::field` | Complete |
| Particle storage (SoA) | `physics::particles` | Complete |
| Cloud-in-cell mass assignment | `physics::cic` | Complete |
| FFT Poisson solver | `physics::poisson` | Complete |
| Zel'dovich initialization | `physics::initial` | Complete |
| Symplectic KDK leapfrog | `physics::integrator` | Complete |
| Conservation diagnostics | `physics::diagnostics` | Complete |
| Simulation driver | `physics::simulation` | Complete |
| Pipeline threading | `run::pipeline` | Complete |
| Snapshot I/O (bincode) | `io::snapshot` | Complete |
| Observer pattern | `io::observer` | Complete |
| Live 3D viewer | `run::pipeline` (vis) | Complete |
| Playback viewer | `run::pipeline` (vis) | Complete |
| GIF recording | `run::runner` | Complete |
| CLI with clap | `run::cli` | Complete |
| Scene system | `scenes` | Extensible |

### What's Next

| Component | Description |
|-----------|-------------|
| Baryonic gas (Euler hydro) | MUSCL-Hancock finite-volume on the existing grid |
| Scale-1 zoom patches | MHD with bivector B-field, restriction/prolongation |
| Learned closures | Rate-matrix-parametrized neural network corrections |
| Retained-mode rendering | GPU point cloud mesh for >100K particles at 60fps |

## Module Layout

```
src/
  algebra.rs          shared Euclidean 3-metric, morphis conversions
  colormap.rs         density/velocity colormapping
  config.rs           TOML configuration with deep merge
  error.rs            HermesError enum
  physics/
    constants.rs      physical constants (kpc, M_sun, Gyr, eV)
    cosmology.rs      FLRW background, growth factor, kick/drift factors
    grid.rs           periodic cubic grid geometry
    field.rs          grade-0 and grade-1 fields with morphis extraction
    particles.rs      SoA particle storage with morphis interface
    cic.rs            cloud-in-cell mass assignment + force interpolation
    poisson.rs        FFT Poisson solver (ndrustfft)
    initial.rs        Zel'dovich ICs from Eisenstein-Hu power spectrum
    integrator.rs     symplectic KDK leapfrog
    diagnostics.rs    conservation audits (mass, momentum, energy, L bivector)
    simulation.rs     simulation driver
  io/
    snapshot.rs       Snapshot type, bincode serialization
    observer.rs       Observer trait, FileObserver, MemoryObserver
  run/
    cli.rs            clap-based CLI
    runner.rs         mode routing (headless, live, playback, record)
    pipeline.rs       threaded pipeline: router, disk writer, precompute, viewer
  scenes/
    cosmic_web.rs     Zel'dovich PM in a periodic box (default scene)
  vis/                (#[cfg(feature = "vis")])
    viewer.rs         static 3D particle viewer (kiss3d)
    plots.rs          density slices, P(k), conservation plots (plotters)
```

## Units

Internal units: kpc, M_sun, Gyr, eV ($k_B = 1$). Constants in `physics::constants` are ported from the [plexis](https://github.com/ctl-alt-leist/plexis) sibling project.

## Documentation

- [Particle-Mesh Method](1_particle-mesh.md) -- the PM force chain and its physics
- [Pipeline Architecture](2_pipeline.md) -- threading model and data flow
- [Configuration](3_configuration.md) -- TOML schema and hierarchy
- [HCD Context](4_hcd-context.md) -- how this connects to the multi-scale framework

## Resources

### Physical Background

- Dodelson & Schmidt, *Modern Cosmology* (2020) -- cosmological perturbation theory
- Hockney & Eastwood, *Computer Simulation Using Particles* (1988) -- PM method
- Springel, *The cosmological simulation code GADGET-2* (2005) -- reference N-body implementation

### Implementation References

- [JaxPM](https://github.com/DifferentiableUniverseInitiative/JaxPM) -- differentiable PM in JAX
- [pmwd](https://github.com/eelregit/pmwd) -- PM with derivatives
- [morphis-rs](https://github.com/ctl-alt-leist/morphis-rs) -- geometric algebra substrate
