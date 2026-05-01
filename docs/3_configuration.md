# Configuration

Hermes uses a four-tier TOML configuration system with deep merge.

## Hierarchy

1. **Embedded defaults** (`configs/defaults.toml`) — always loaded, compiled into the binary
2. **Scene defaults** (`src/scenes/*/defaults.toml`) — scene-specific overrides
3. **Config file** (optional) — user-provided partial TOML, passed as a CLI argument
4. **CLI overrides** (`--steps`, `--particles`, `--grid`) — highest priority

Each tier is deep-merged into the previous: only the fields you want to change need to appear. Unspecified fields inherit from the tier below.

User config files in `configs/` with `.local.toml` extension are gitignored.

## Schema

```toml
[cosmology]
hubble         = 0.674       # dimensionless h (H0 = 100h km/s/Mpc)
omega_m        = 0.315       # total matter density parameter
omega_b        = 0.0493      # baryon density parameter
omega_r        = 9.15e-5     # radiation density parameter
omega_k        = 0.0         # spatial curvature parameter
omega_v        = 0.6849085   # vacuum energy density parameter
sigma_8        = 0.811       # RMS fluctuation amplitude at 8 h^-1 Mpc
spectral_index = 0.965       # primordial power spectrum slope n_s

[simulation]
n_grid         = 32          # grid cells per side (N_g^3 total)
n_particles    = 32          # particles per side (N_p^3 total)
box_length     = 100000.0    # comoving box side in kpc (100 Mpc)

[time]
scale_factor_range    = [0.02, 1.0]  # [initial, final]: z ~ 49 to z = 0
scale_factor_stepping = "log"        # "log" or "linear"
n_steps               = 300

[output]
write_interval      = 1      # save a snapshot every N steps
diagnostic_interval = 10     # compute full diagnostics every N steps

[initialization]
spectrum               = "power"  # "power" (CDM P(k)) or "random" (synthetic)
perturbation_amplitude = 0.1      # RMS density perturbation at initialization
band_pass              = [1.5, 0.5]  # [k_min / k_fundamental, k_max / k_nyquist]

[field]
length_scale = 2000.0    # ell/m: smoothing length ratio (kpc^2 / Gyr)
mass         = 1e10      # field mass parameter (M_sun)

[visualization]
point_size        = 5.0              # screen-space point size (particles)
blob_size         = 18.0             # screen-space blob size (volumetric fields)
blob_alpha        = 0.12             # per-blob opacity for additive blending
blob_falloff      = 10.0             # Gaussian falloff rate for volumetric blobs
camera_distance   = 1.9              # distance from origin (box is [-0.5, 0.5])
camera_angle      = [0.56, 0.42, 0.69]  # direction vector (multiplied by distance)
colormap_range    = [0.3, 3.0]       # [floor, ceiling] as density / rho_mean
jitter            = 0.3              # grid-point jitter as fraction of cell size
gif_resolution    = 512              # pixel size for GIF recording
gif_point_radius  = 1                # point radius in GIF frames
```

## Validation

The cosmology section is validated at load time:

- `hubble`, `omega_m`, `omega_b`, `sigma_8` must be positive
- `omega_v`, `omega_r` must be non-negative
- `omega_b <= omega_m` (baryons can't exceed total matter)
- $Ω_v + Ω_k + Ω_m + Ω_r = 1$ within $10^{-6}$

## Scene Defaults

Each scene provides a `defaults.toml` that overrides only what differs from the global defaults. The cosmic-web scene sets 64^3 resolution and the z=49 starting point. The fuzzy-dm scene sets "random" spectrum and linear stepping. Scene defaults are merged between the global defaults and the user config file.

## CLI Flags

| Flag | Overrides |
|------|-----------|
| `--steps N` | `time.n_steps` |
| `--particles N` | `simulation.n_particles` |
| `--grid N` | `simulation.n_grid` |
| `--scene NAME` | selects scene defaults |
| `--seed N` | RNG seed (not in config) |

## Custom Config File

Create a TOML file with only the fields you want to override:

```toml
# configs/my-run.local.toml
[time]
scale_factor_range = [0.5, 1.5]
n_steps = 200
```

Run with: `hermes --scene cosmic-web configs/my-run.local.toml`

## Implementation

- `config::load_defaults()` — parse embedded defaults
- `config::build_configuration(file, overrides)` — four-tier deep merge
- Deep merge via recursive `toml::Value::Table` traversal: tables merge recursively, scalars overwrite
- `VisualizationConfig`, `InitializationConfig`, `FieldConfig` have `Default` impls so they can be omitted from partial config files
