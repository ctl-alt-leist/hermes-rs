# Unified Physics Engine

## Vision

A simulation is defined by its content (what is in the box) and its physics (what laws act on it). The engine should not need to be told which scene it came from — it should look at what's in the snapshot, look at which physics modules are enabled, and evolve the system forward.

The current architecture couples content to dynamics through the scene: cosmic-web produces particles and selects ParticleMeshDynamics, fuzzy-dm produces fields and selects SchrodingerPoissonDynamics. The scene is the matchmaker. The goal is to remove that coupling so the engine can operate on any content with any combination of physics.

## Content as State

The simulation state is its content — particles, fields, or both. A snapshot fully describes this state at a moment in time. To resume, you load a snapshot and apply physics. No scene required.

Current snapshot carries:
- Particle positions, momenta, mass per particle
- Field density on a grid
- Scale factor, step number

What it needs to also carry:
- Grid parameters (n_cells, box_length) — the Poisson solver needs these
- Cosmology parameters — H(a), density_mean depend on these
- Content type tag — so the engine knows what it's looking at
- Optionally: the full wavefunction (scalar + pseudoscalar components), not just density, for field content that needs to be resumed

With this metadata, a snapshot becomes a self-describing initial condition. Point the engine at any snapshot file, give it a time range, and it runs.

## Physics as Modules

Instead of a single `Dynamics` trait that owns the whole step, physics is decomposed into modules that each know how to act on content:

```
Gravity          — Poisson solve, applies to both particles and fields
QuantumPressure  — kinetic step of SP integrator, applies to fields only
Hydrodynamics    — (future) baryon pressure, cooling, applies to baryon fields
Electromagnetism — (future) Maxwell equations, applies to EM bivector field
```

Each module has a simple interface:
- Can it act on this content? (particles, fields, mixed)
- Apply one step given the current state and dt

The engine composes modules: if gravity is on and quantum pressure is on, the step is a split-step that interleaves them. If gravity is on and quantum pressure is off, it's pure PM. The composition logic lives in the engine, not in the scenes.

Configuration would look like:

```toml
[physics]
gravity          = true
quantum_pressure = true
hydrodynamics    = false
electromagnetic  = false
```

Scenes become convenience presets that set content type, physics modules, and default parameters — but they're not load-bearing. A resumed simulation doesn't need a scene; it needs a snapshot and a physics configuration.

## Implications for Snapshots

The snapshot format needs to be self-describing. One approach:

```
Header:
  format_version: u32
  content_type: "particles" | "fields" | "mixed"
  grid: { n_cells, box_length }
  cosmology: { hubble, omega_m, omega_v, ... }
  scale_factor: f64
  step: usize
  
Body:
  (content-type-dependent data)
```

For field content, the body must store the full wavefunction (both components of the even field), not just the density. Density is a derived quantity — you can't resume a Schrodinger-Poisson simulation from density alone because the phase information is lost.

For particle content, positions and momenta suffice. The KDK leapfrog is self-starting from a single state (no need for two snapshots to infer momentum — momentum is stored directly as canonical momentum p = a² m dx/dt).

## Two-Snapshot Question

Some integrators require information from two consecutive states to restart (e.g., Verlet without stored velocities, or methods that cache forces from the previous step). The current KDK leapfrog caches `forces_prev` for reuse, but this is an optimization, not a requirement — the first step of a resumed run simply recomputes forces from scratch, which is correct (just slightly less efficient for one step).

The Schrodinger-Poisson split-step is fully self-starting from a single wavefunction state. No two-snapshot issue there.

If a future integrator genuinely requires two states, the snapshot format could store both, or the engine could accept a pair of snapshots. This is a bridge to cross when we get there.

## Relationship to the Scale Hierarchy

The full HCD framework envisions multiple scales running simultaneously, with restriction/prolongation operators coupling them. Each scale has its own content and its own physics modules. The unified engine at a single scale is the building block — the multi-scale orchestrator composes engines at different resolutions.

The snapshot format should anticipate this: a multi-scale snapshot is a collection of single-scale snapshots at different resolutions, each self-describing. But single-scale operation is the foundation and should work independently.

## Current State

- Content abstraction exists (Particles / Fields / Mixed enum)
- Dynamics trait exists (step method takes Content + Cosmology)
- Two dynamics implementations: ParticleMeshDynamics, SchrodingerPoissonDynamics
- Resume works for particles via --resume flag, but requires --scene to select dynamics
- Snapshots carry content but not grid/cosmology metadata
- Physics modules are monolithic (each Dynamics impl owns its full step)

## Next Steps (not prioritized)

1. Add grid and cosmology metadata to the snapshot format
2. Store full wavefunction (not just density) for field snapshots
3. Factor gravity out of both PM and SP dynamics into a shared module
4. Make --resume infer dynamics from snapshot content type
5. Add physics module configuration to TOML
