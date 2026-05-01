# Project Overview

Hermes is the dynamical simulation layer of the Hierarchical Closure Dynamics (HCD) framework. Its algebraic substrate is [morphis-rs](https://github.com/ctl-alt-leist/morphis-rs), a Rust geometric algebra library. morphis handles elements, products, and linear maps; hermes handles grids, time integration, scale hierarchies, and closure terms.

## Scope

The current implementation supports two dynamical modes at scale-0:

- **Particle-mesh (PM)**: Dark matter N-body with CIC mass assignment, FFT-Poisson gravity, and symplectic KDK leapfrog. Zel'dovich initial conditions from the linear power spectrum.
- **Schrodinger-Poisson (SP)**: Fuzzy dark matter wavefunction evolution via split-step spectral integrator. Even-subalgebra fields from morphis with self-gravity through the Poisson equation.

Both are instantiations of the generator-equation form

$$
\dot{\Psi}^{(s)} = G^{(s)} \Psi^{(s)} + S^{(s)} + C^{(s)} + A^{(s)}
$$

with closure $C^{(0)} = 0$, sources $S^{(0)} = 0$, and the gravitational generator $G^{(0)}$ built from the FFT-Poisson force chain.

The content abstraction (Particles / Fields / Mixed) allows both modes to coexist in the same framework, with pluggable dynamics selected per scene.

### What's Built

| Component | Module | Status |
|-----------|--------|--------|
| FLRW cosmological background | `physics::cosmology` | Complete |
| Periodic cubic grid | `physics::grid` | Complete |
| Scalar and vector fields | `physics::field` | Complete |
| Particle storage (SoA) | `physics::particles` | Complete |
| Cloud-in-cell mass assignment | `physics::cic` | Complete |
| FFT Poisson solver | `physics::poisson` | Complete |
| Symplectic KDK leapfrog | `physics::integrator` | Complete |
| Conservation diagnostics | `physics::diagnostics` | Complete |
| Content abstraction | `core::content` | Complete |
| Dynamics trait | `core::dynamics` | Complete |
| PM dynamics | `core::pm_dynamics` | Complete |
| Schrodinger-Poisson dynamics | `core::schrodinger_dynamics` | Complete |
| Simulation driver | `core::simulation` | Complete |
| Cosmic web scene (PM) | `scenes::cosmic_web` | Complete |
| Galaxy group scene (PM) | `scenes::galaxy_group` | Complete |
| Fuzzy dark matter scene (SP) | `scenes::fuzzy_dm` | Complete |
| Pipeline threading | `run::pipeline` | Complete |
| Snapshot I/O (bincode) | `io::snapshot` | Complete |
| Observer pattern | `io::observer` | Complete |
| Live 3D viewer | `run::pipeline` (vis) | Complete |
| Playback viewer with controls | `run::pipeline` (vis) | Complete |
| Volumetric field renderer | `visuals::volumetric_renderer` (vis) | Complete |
| GIF recording | `run::runner` | Complete |
| CLI with clap | `run::cli` | Complete |
| Scene system | `scenes` | Extensible |
| Resume from snapshot | `run::runner` | Complete |
| Configurable initialization | `config::InitializationConfig` | Complete |

### What's Next

| Component | Description |
|-----------|-------------|
| Unified physics engine | Content-driven module composition (see `docs/unified-engine.md`) |
| Self-describing snapshots | Metadata headers for scene-independent resume |
| Baryonic gas (Euler hydro) | MUSCL-Hancock finite-volume on the existing grid |
| Scale-1 zoom patches | MHD with bivector B-field, restriction/prolongation |
| Learned closures | Rate-matrix-parametrized neural network corrections |
| Volume raycasting | wgpu-based 3D density field rendering |

## Module Layout

```
src/
  algebra.rs          shared Euclidean 3-metric, morphis conversions
  colormap.rs         density/velocity colormapping
  config.rs           TOML configuration with deep merge
  error.rs            HermesError enum
  core/               simulation orchestration
    content.rs        Content enum (Particles / Fields / Mixed)
    dynamics.rs       Dynamics trait (pluggable step function)
    pm_dynamics.rs    ParticleMeshDynamics (KDK + Poisson)
    schrodinger_dynamics.rs  SchrodingerPoissonDynamics (split-step spectral)
    simulation.rs     simulation driver (from_scene, from_config, resume)
  physics/            physical models and solvers
    constants.rs      physical constants (kpc, M_sun, Gyr, eV)
    cosmology.rs      FLRW background, growth factor, kick/drift factors
    grid.rs           periodic cubic grid geometry
    field.rs          grade-0 and grade-1 fields with morphis extraction
    particles.rs      SoA particle storage with morphis interface
    cic.rs            cloud-in-cell mass assignment + force interpolation
    poisson.rs        FFT Poisson solver (ndrustfft)
    integrator.rs     symplectic KDK leapfrog
    diagnostics.rs    conservation audits (mass, momentum, energy, L bivector)
  io/
    snapshot.rs       Snapshot type, SnapshotContent enum, bincode serialization
    observer.rs       Observer trait, FileObserver, MemoryObserver
  run/
    cli.rs            clap-based CLI (--scene, --live, --playback, --resume, etc.)
    runner.rs         mode routing (headless, live, playback, record, resume)
    pipeline.rs       threaded pipeline: router, disk writer, precompute, viewer
  scenes/             each a subdirectory with init.rs + defaults.toml
    cosmic_web/       Zel'dovich PM in a 100 Mpc periodic box (default)
    galaxy_group/     3 colliding NFW halos in an 8 Mpc box
    fuzzy_dm/         Schrodinger-Poisson wavefunction in a 10 Mpc box
  visuals/            (#[cfg(feature = "vis")])
    viewer.rs         static 3D particle viewer (kiss3d)
    plots.rs          density slices, P(k), conservation plots (plotters)
    volumetric_renderer.rs  additive-blended Gaussian point sprites for fields
```

## Units

Internal units: kpc, M_sun, Gyr, eV ($k_B = 1$). Constants in `physics::constants` are ported from the [plexis](https://github.com/ctl-alt-leist/plexis) sibling project.

## Documentation

- [Particle-Mesh Method](1_particle-mesh.md) -- the PM force chain and its physics
- [Pipeline Architecture](2_pipeline.md) -- threading model and data flow
- [Configuration](3_configuration.md) -- TOML schema and hierarchy
- [HCD Context](4_hcd-context.md) -- how this connects to the multi-scale framework
- [Morphis Fields](5_morphis-fields.md) -- field-theoretic formulation and morphis integration
- [Unified Engine Plan](a_plans-for-a-unified-engine.md) -- content-driven physics composition
- [Efficient Snapshots Plan](b_plan-for-efficient-snapshots.md) -- I/O bottleneck analysis

## Resources

### Physical Background

- Dodelson & Schmidt, *Modern Cosmology* (2020) -- cosmological perturbation theory
- Hockney & Eastwood, *Computer Simulation Using Particles* (1988) -- PM method
- Springel, *The cosmological simulation code GADGET-2* (2005) -- reference N-body implementation
- Schive et al., *Understanding the Core-Halo Relation of Quantum Wave Dark Matter* (2014) -- FDM reference

### Implementation References

- [JaxPM](https://github.com/DifferentiableUniverseInitiative/JaxPM) -- differentiable PM in JAX
- [pmwd](https://github.com/eelregit/pmwd) -- PM with derivatives
- [morphis-rs](https://github.com/ctl-alt-leist/morphis-rs) -- geometric algebra substrate
