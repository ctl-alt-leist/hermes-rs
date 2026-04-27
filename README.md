# Hermes

Hermes is a cosmological particle-mesh simulator built on [morphis](https://github.com/ctl-alt-leist/morphis-rs), a geometric algebra library. Named for the messenger god who moves between scales, Hermes implements the dynamical simulation layer of the Hierarchical Closure Dynamics framework -- a multi-scale approach to cosmological structure formation where coarse-grained state spaces are connected bidirectionally to fine-grained patches via restriction and prolongation operators.

## Features

- **Particle-mesh N-body**: FFT-based Poisson gravity with cloud-in-cell mass assignment
- **Geometric algebra substrate**: all physical quantities are morphis objects -- positions, momenta, forces are grade-1 vectors; angular momentum is a grade-2 bivector via the wedge product
- **Symplectic integration**: kick-drift-kick leapfrog with cosmological step factors
- **Zel'dovich initialization**: Eisenstein-Hu transfer function, sigma-8-normalized power spectrum
- **Pipeline architecture**: simulation, disk I/O, and visualization run on independent threads connected by bounded channels
- **Live 3D viewer**: real-time particle visualization during simulation via kiss3d
- **TOML configuration**: three-tier hierarchy with deep merge (defaults, presets, overrides)
- **Snapshot I/O**: bincode serialization with morphis vector roundtrip fidelity

## Quick Start

```bash
# Run a simulation (saves snapshots to data/<timestamp>/)
cargo run --release -- --particles 32 --steps 200

# Run with live 3D viewer
cargo run --release --features vis -- --live --particles 32 --steps 300

# Play back saved snapshots
cargo run --release --features vis -- --playback data/<dir>

# Record to GIF
cargo run --release -- --playback data/<dir> --record output/cosmic-web.gif --fps 20

# See all options
cargo run --release -- --help
```

## CLI

```
hermes [OPTIONS] [CONFIG_FILE]

Arguments:
  [CONFIG_FILE]        TOML config file (overrides defaults)

Options:
  --scene <NAME>       Simulation scene [default: cosmic-web]
  --live               Open live 3D viewer (requires --features vis)
  --save [DIR]         Save snapshots (default: data/<timestamp>/)
  --no-save            Don't save snapshots
  --playback DIR       Play back saved snapshots
  --record FILE        Record playback as GIF
  --fps N              Playback/recording framerate [default: 15]
  --seed N             RNG seed [default: 42]
  --steps N            Override time steps
  --particles N        Override particles per side
  -q, --quiet          Suppress output
```

## Documentation

- [Project Overview](docs/0_project-overview.md) -- scope, architecture, and reading guide
- [Particle-Mesh Method](docs/1_particle-mesh.md) -- the PM force chain and its implementation
- [Pipeline Architecture](docs/2_pipeline.md) -- threading, channels, and data flow
- [Configuration](docs/3_configuration.md) -- TOML schema and config hierarchy
- [HCD Context](docs/4_hcd-context.md) -- how hermes connects to the multi-scale framework

## Development

| Command | Description |
|---------|-------------|
| `make lint` | `cargo fmt` + `cargo clippy --fix` |
| `make test` | `cargo test` |
| `make build` | `cargo build --release` |
| `make clean` | `cargo clean` + remove `.DS_Store` |

Pre-commit hooks enforce `cargo fmt --check`, `cargo clippy -- -D warnings`, and `cargo test` on every commit.

### Code Style

- Rust 2024 edition, `rustfmt` defaults, `clippy` clean
- All tests are integration tests in `tests/`
- `thiserror` for error types, `serde` + `toml` for config, `bincode` for snapshots
- morphis geometric algebra for all physical quantities
- Descriptive variable names (`density_mean`, `mass_particle`), not mathematical symbols

## License

MIT

---

Built with [morphis](https://github.com/ctl-alt-leist/morphis-rs) geometric algebra. Developed with [Claude Code](https://claude.ai/code).
