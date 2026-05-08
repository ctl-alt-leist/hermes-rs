# Species Specification: Design Notes

Notes on the eventual design for user-defined species — the layer above the plumbing. The multi-species pipeline (snapshots, visualization, engine composition) is being built now; this document collects first thoughts on the physics-specification layer that will sit on top of it.

The goal is that a user writes a TOML config declaring an arbitrary collection of named species — fields and particles — each with its geometry, dynamics, couplings, and visualization. Nothing is hardcoded to "alpha" or "dark matter." The framework dispatches initialization, evolution, and rendering from the config alone.

## Fields: Geometry and Lagrangian

A field species should be fully specified by:

1. **Algebraic grade** — already in config as `grade = [0, 3]` (even subalgebra), `grade = 0` (scalar), `grade = 2` (bivector). This determines what kind of object lives at each grid point.

2. **Lagrangian form** — first-degree (Schrodinger-type, norm-preserving) or second-degree (wave-equation, propagating). Currently handled by the `free` string ("schrodinger", "wave"). This could become richer: the form determines the sector type, the elementary flows, and the conservation laws.

3. **Self-interaction** — currently `self_interaction` for Gross-Pitaevskii. More general: any scalar function of the field's norm, specified as a potential `V(|psi|^2)`. The simplest cases (quadratic mass, quartic GP) cover most near-term needs. Whether to support arbitrary user-defined potentials or a menu of named options is an open question.

4. **Name and symbol** — the TOML key is the name. A display symbol (Unicode Greek letter) could be added for visualization labels and diagnostics output. E.g., `symbol = "α"` in the config, used by the pipeline for legends and by diagnostics for conservation readouts.

The science notes (section 1.13, "From Geometric Specification to Lagrangian") lay out a seven-step recipe for constructing a sector's Lagrangian from its grade and desired physics. The question is whether this recipe should be encoded in the config (declarative: "I want a first-degree even-subalgebra field with quartic self-interaction") or remain a paper-to-code workflow where the user writes a Sector implementation. The declarative approach is cleaner for standard cases; custom Sector implementations handle anything exotic.

## Particles: Properties and Forces

A particle species should be specified by:

1. **Count and mass** — already in config.
2. **Charge** — already in config (for electromagnetic coupling, unused so far).
3. **Softening** — gravitational softening length, needed for particle-particle interactions. Currently implicit in the Poisson solver's Green's function; should be per-species and explicit.
4. **Deposition kernel** — CIC, TSC, or higher. Currently hardcoded to CIC. Should be per-species or at least per-simulation.

Particles don't have Lagrangians in the same sense as fields — their dynamics are Hamilton's equations with forces interpolated from the grid. The "evolution rule" for particles is the force law, which comes from the coupling structure rather than from a per-species Lagrangian. The interesting physics questions are about what forces act on each particle species and how those forces are computed.

## Couplings: Which Species Talk to Which

The current config has `gravity = true` (universal, all massive species) and `electromagnetic = ["beta", "gamma"]` (selective). This is the right structure but needs to grow.

Near-term additions:

- **Gravitational mass fraction** — when multiple field species share a Poisson potential, how much of the total density does each contribute? Currently computed from `mass * |psi|^2`; for cosmological splits (dark matter vs baryons), the mass parameters handle this naturally. But if two species have the same mass parameter, we might want explicit density weights.

- **Coupling constants** — the gravitational constant G is universal, but other couplings (electromagnetic, baryon-magnetic) have per-pair coupling strengths. The config should specify these on the coupling entry.

- **Coupling topology** — which species source which potentials. The current architecture (all sectors deposit into one gravity solver) works for universal gravity. For selective couplings (only baryons source the magnetic field), we need a way to say "species X sources coupling Y but species Z does not." The `electromagnetic = ["beta", "gamma"]` pattern is a start.

Longer-term, the coupling specification merges with the Lagrangian specification: the joint Lagrangian has free terms (per-species) and interaction terms (per-pair), and the config declares both. The engine reads the joint Lagrangian and constructs the appropriate sectors and solvers. This is the full vision from the science notes; the near-term steps are incremental moves toward it.

## Visualization: Per-Species Rendering

Already implemented: per-species colormaps via `[output.display.species.<name>]`. Future additions:

- **Display symbol** — a Unicode character shown in legends and diagnostics. `symbol = "α"` in the species config.
- **Render mode per species** — some species might render as volumetric blobs, others as points, others as streamlines (for vector/bivector fields). Currently all fields are volumetric and all particles are points; the render mode could be per-species.
- **Opacity and size per species** — override `blob_alpha`, `blob_size` per species for visual emphasis (e.g., dim the dark matter, brighten baryons).
- **Slice rendering** — 2D density slices through the box, one per species, as a diagnostic output.

## What We Build Now vs. Later

**Now (plumbing):** Unified snapshot format with named species, `Snapshot::capture_from_state`, per-species particle visualization, ParticleSector, retire hardcoded Content types from the engine path. This makes arbitrary named species flow end-to-end without caring what physics they carry.

**Next conversation:** Lagrangian-form dispatch from config (declarative sector construction for standard cases), coupling topology in config, per-species softening and deposition kernels for particles.

**Later:** Full Lagrangian specification in config or code, gauge couplings, non-Abelian structure, sub-grid closures from the multi-scale machinery.
