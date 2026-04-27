//! Command-line interface for hermes.

use clap::Parser;

/// Hermes — Cosmological particle-mesh simulation.
#[derive(Parser, Debug)]
#[command(name = "hermes", about = "Hierarchical Closure Dynamics simulator")]
pub struct Cli {
    /// TOML config file (overrides defaults).
    pub config_file: Option<String>,

    /// Simulation scene.
    #[arg(long, default_value = "cosmic_web")]
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

    /// RNG seed.
    #[arg(long, default_value = "42")]
    pub seed: u64,

    /// Override number of time steps.
    #[arg(long)]
    pub steps: Option<usize>,

    /// Override particles per side.
    #[arg(long)]
    pub particles: Option<usize>,

    /// Suppress terminal output.
    #[arg(short, long)]
    pub quiet: bool,
}

impl Cli {
    /// Resolve the save directory.
    ///
    /// - `--save dir` → use that dir
    /// - `--save` (no arg) → `data/<timestamp>/`
    /// - no flag and no `--no-save` → `data/<timestamp>/` (save by default)
    /// - `--no-save` → None
    pub fn save_directory(&self) -> Option<String> {
        if self.no_save {
            return None;
        }

        match &self.save {
            Some(Some(dir)) => Some(dir.clone()),
            Some(None) | None => Some(timestamped_dir()),
        }
    }
}

fn timestamped_dir() -> String {
    let now = chrono::Local::now();

    format!("data/{}", now.format("%Y-%m-%d_%H%M%S"))
}
