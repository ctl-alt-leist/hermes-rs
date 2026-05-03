//! Command-line interface for hermes.

use clap::Parser;

/// Hermes — Cosmological particle-mesh simulation.
#[derive(Parser, Debug)]
#[command(name = "hermes", about = "Hierarchical Closure Dynamics simulator")]
pub struct Cli {
    /// TOML config file (overrides defaults).
    pub config_file: Option<String>,

    /// Simulation scene.
    #[arg(long, default_value = "cosmic-web-pm")]
    pub scene: String,

    /// Open live 3D viewer during simulation (requires --features vis).
    #[arg(long)]
    pub live: bool,

    /// Save snapshots to directory (default: data/<timestamp>/).
    #[arg(long)]
    pub save: Option<Option<String>>,

    /// Don't save snapshots.
    #[arg(long, conflicts_with = "save")]
    pub no_save: bool,

    /// Play back saved snapshots from directory (no simulation).
    #[arg(long)]
    pub playback: Option<String>,

    /// Save playback as GIF (e.g. --record output.gif). Use with --playback.
    #[arg(long)]
    pub record: Option<String>,

    /// Playback framerate in fps.
    #[arg(long, default_value = "30")]
    pub fps: u64,

    /// RNG seed.
    #[arg(long, default_value = "42")]
    pub seed: u64,

    /// Override number of time steps.
    #[arg(long)]
    pub steps: Option<usize>,

    /// Override particles per side (N_p).
    #[arg(long)]
    pub particles: Option<usize>,

    /// Override grid cells per side (N_g).
    #[arg(long)]
    pub grid: Option<usize>,

    /// Resume simulation from the last snapshot in a directory.
    #[arg(long)]
    pub resume: Option<String>,

    /// Suppress terminal output.
    #[arg(short, long)]
    pub quiet: bool,
}

impl Cli {
    /// Resolve the save directory.
    ///
    /// - `--save dir` → use that dir
    /// - `--save` (no arg) → `data/<scene>/`
    /// - no flag and no `--no-save` → `data/<scene>/` (save by default)
    /// - `--no-save` → None
    pub fn save_directory(&self) -> Option<String> {
        if self.no_save {
            return None;
        }

        match &self.save {
            Some(Some(dir)) => Some(dir.clone()),
            Some(None) | None => Some(format!("data/{}", self.scene)),
        }
    }
}
