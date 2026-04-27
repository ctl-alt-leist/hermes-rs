# Configuration

Hermes uses a three-tier TOML configuration system with deep merge, following the same pattern as [plexis](https://github.com/ctl-alt-leist/plexis).

## Hierarchy

1. **Embedded defaults** (`configs/defaults.toml`) -- always loaded, compiled into the binary
2. **Config file** (optional) -- partial TOML overriding any defaults
3. **CLI overrides** (`--steps`, `--particles`) -- highest priority

Each tier is deep-merged into the previous: only the fields you want to change need to appear. Unspecified fields inherit from the tier below.

## Schema

```toml
[cosmology]
hubble         = 0.674       # dimensionless h (H0 = 100h km/s/Mpc)
omega_m        = 0.315       # total matter density parameter
omega_b        = 0.0493      # baryon density parameter
omega_r        = 9.15e-5     # radiation density parameter
omega_k        = 0.0         # spatial curvature parameter
omega_lambda   = 0.6849085   # vacuum energy (cosmological constant)
sigma_8        = 0.811       # RMS fluctuation amplitude at 8 h^-1 Mpc
spectral_index = 0.965       # primordial power spectrum slope n_s

[simulation]
n_cells      = 64            # grid cells per side (total = n_cells^3)
n_particles  = 64            # particles per side (total = n_particles^3)
box_length   = 100000.0      # comoving box side in kpc

[time]
scale_factor_initial = 0.02  # start at z ~ 49
scale_factor_final   = 1.0   # end at z = 0
n_steps              = 200   # number of time steps
stepping             = "log_a"  # "log_a" or "linear_a"

[output]
directory         = "output"
snapshot_interval = 10       # diagnostics every N steps
```

## Validation

The cosmology section is validated at load time:

- `hubble`, `omega_m`, `omega_b`, `sigma_8` must be positive
- `omega_lambda`, `omega_r` must be non-negative
- `omega_b <= omega_m` (baryons can't exceed total matter)
- $Ω_Λ + Ω_k + Ω_m + Ω_r = 1$ within $10^{-6}$

## Custom Config File

Create a TOML file with only the fields you want to override:

```toml
# my-run.toml -- a smaller, faster run
[simulation]
n_cells    = 16
n_particles = 16

[time]
n_steps = 50
```

Run with: `cargo run --release -- my-run.toml`

## Presets

The defaults ship Planck 2018 best-fit cosmology. Different cosmologies (WMAP, Einstein-de Sitter) can be specified as complete override files:

```toml
# einstein-de-sitter.toml
[cosmology]
hubble       = 0.7
omega_m      = 1.0
omega_b      = 0.05
omega_r      = 0.0
omega_k      = 0.0
omega_lambda = 0.0
sigma_8      = 0.811
spectral_index = 1.0
```

## Implementation

- `config::load_defaults()` -- parse embedded defaults
- `config::build_configuration(file, overrides)` -- three-tier deep merge
- Deep merge via recursive `toml::Value::Table` traversal: tables merge recursively, scalars overwrite
