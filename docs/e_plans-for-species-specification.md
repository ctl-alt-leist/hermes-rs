# Species Specification: Design Notes

Notes on the design for user-defined species — the layer above the plumbing.

The goal: a user writes a TOML config declaring an arbitrary collection of named species — fields
and particles — each with its geometry, dynamics, couplings, and visualization. Nothing is hardcoded.
The framework dispatches initialization, evolution, and rendering from the config alone.

## Naming

Two levels of identity per species:

- **Name** — the TOML key. Human-readable, used as the BTreeMap key everywhere in the code, in
  snapshot species labels, in coupling references, in visualization config lookups.
  Examples: `"dark matter"`, `"baryonic matter"`, `"magnetic field"`.

- **Symbol** — a Unicode display character for diagnostics and visualization labels.
  Examples: `"α"`, `"β"`, `"γ"`, `"φ"`.

```toml
[ontology.fields."dark matter"]
symbol = "α"
grade  = [0, 3]
n      = 64
# ...
```

## Fields

A field species is specified by:

```toml
[ontology.fields."dark matter"]
symbol           = "α"
grade            = [0, 3]          # even subalgebra
n                = 64              # cells per side (n³ grid points)
mass             = 1e10            # M_sun
length_scale     = 4000.0          # l/m: diffusivity (kpc² / Gyr)
free             = "schrodinger"   # Lagrangian form
self_interaction = 1e6             # Gross-Pitaevskii coupling (optional)
```

- **`grade`** — algebraic grade in the geometric algebra. Determines the target subspace.
- **`free`** — Lagrangian form: `"schrodinger"` (first-degree, norm-preserving) or
  `"wave"` (second-degree, propagating). Determines which Sector implementation the engine
  constructs.
- **`self_interaction`** — optional quartic (Gross-Pitaevskii) coupling constant. When present,
  the engine uses a GrossPitaevskiiSector instead of a SchrodingerSector.
- **`n`** — grid cells per side (total grid points = n³). Same parameter name as particles.

The engine reads `free` and `self_interaction` to dispatch the correct sector type automatically.
Custom Sector implementations handle anything beyond the standard menu.

## Particles

A particle species is specified by:

```toml
[ontology.particles."dark matter"]
symbol    = "α"
n         = 64           # particles per side (n³ total)
mass      = 1e10         # M_sun per particle
softening = 50.0         # gravitational softening length (kpc)
kernel    = "cic"        # deposition kernel: "cic", "tsc", "pcs"
```

- **`n`** — particles per side (total count = n³). Consistent with field `n`.
- **`softening`** — spatial extent of the particle. Applied in the Poisson solver's Green's
  function. Per-species because different populations can have different effective sizes.
- **`kernel`** — deposition/interpolation kernel. Per-species for flexibility, though most
  simulations will use the same kernel for all species. Default: `"cic"`.

## Couplings

Couplings are declared as a list of interaction terms. Each entry names its kind, lists
participating species, and carries coupling-specific parameters.

```toml
[[ontology.coupling]]
kind    = "gravity"
species = ["dark matter", "baryonic matter"]

[[ontology.coupling]]
kind    = "electromagnetic"
species = [
    { name = "baryonic matter", charge = 1.0 },
    { name = "magnetic field" },
]
```

- **`gravity`** — universal Poisson coupling. All listed species deposit density and feel the
  resulting potential. The gravitational constant G is universal; no per-pair parameters needed.
- **`electromagnetic`** — selective coupling with per-species charge. Charge lives on the coupling
  entry, not on the species, because it describes how a species participates in a specific
  interaction.

Future coupling kinds: baryon pressure, magnetic sourcing, custom closure terms.

## Visualization

Per-species display config via `[output.display.species."<name>"]`:

```toml
[output.display.species."dark matter"]
colormap       = "cool"
colormap_range = [0.3, 3.0]
blob_size      = 28.0
blob_alpha     = 0.08
render_mode    = "volumetric"

[output.display.species."baryonic matter"]
colormap    = "ember"
blob_alpha  = 0.15
render_mode = "volumetric"
```

- **`colormap`** — named colormap: `"hot"`, `"cool"`, `"ember"`, `"verdant"`.
- **`colormap_range`** — density floor/ceiling in units of rho_mean. Falls back to global default.
- **`blob_size`**, **`blob_alpha`** — per-species rendering overrides.
- **`render_mode`** — `"volumetric"` (blobs), `"points"`, or `"streamlines"` (future).
  Default: `"volumetric"` for fields, `"points"` for particles.

## Sector Dispatch

The engine constructs sectors automatically from the config:

| `free` value     | `self_interaction` | Sector type              |
| ---------------- | ------------------ | ------------------------ |
| `"schrodinger"`  | absent             | `SchrodingerSector`      |
| `"schrodinger"`  | present            | `GrossPitaevskiiSector`  |
| `"wave"`         | —                  | `WaveSector` (future)    |
| (particles)      | —                  | `ParticleSector`         |

This dispatch table is the bridge between declarative config and the engine's trait-based
architecture. Adding a new sector type means adding a row to the table and implementing the
Sector trait.

## What to Build Now

1. Add `name` (from TOML key) and `symbol` to `FieldSpecies` and `ParticleSpecies`.
2. Change particle `count` to `n` (per side, total = n³). Match field `n`.
3. Remove `charge` from species configs (move to coupling entries when EM is built).
4. Add `blob_size`, `blob_alpha`, `render_mode` to `SpeciesDisplayConfig`.
5. Wire symbol through to diagnostics and print statements.
6. Build sector dispatch from config (auto-construct sectors from `free` + `self_interaction`).
7. Restructure `[ontology.lagrangian]` into `[[ontology.coupling]]` list.
8. Update all scene TOMLs.
