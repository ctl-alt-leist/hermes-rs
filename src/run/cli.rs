//! Command-line interface for hermes.

use clap::Parser;

/// Hermes — Hierarchical Closure Dynamics simulator.
#[derive(Parser, Debug)]
#[command(name = "hermes")]
pub struct Cli {
    /// Scene config: path to a .toml file (with or without extension).
    #[arg(long, default_value = "scenes/cosmic-web-pm")]
    pub scene: String,

    /// Base config file merged under the scene config (optional overrides).
    #[arg(long)]
    pub config: Option<String>,

    /// Save snapshots (default dir: next to scene TOML; or specify a path).
    #[arg(long)]
    pub save: Option<Option<String>>,

    /// Open live 3D viewer during simulation (requires --features vis).
    #[arg(long)]
    pub live: bool,

    /// Play back saved snapshots from directory (requires --features vis).
    #[arg(long)]
    pub playback: Option<String>,

    /// Save playback as GIF (use with --playback).
    #[arg(long)]
    pub record: Option<String>,

    /// Playback/recording framerate.
    #[arg(long, default_value = "30")]
    pub fps: u64,

    /// RNG seed.
    #[arg(long, default_value = "42")]
    pub seed: u64,

    /// Override number of time steps.
    #[arg(long)]
    pub steps: Option<usize>,

    /// Override particles per side.
    #[arg(long)]
    pub particles: Option<usize>,

    /// Override grid cells per side.
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
    /// Resolve the save directory. None if --save was not passed.
    ///
    /// - `--save dir` → use that dir
    /// - `--save` (no arg) → directory next to the scene TOML (strip .toml)
    /// - no --save flag → None (don't save)
    pub fn save_directory(&self) -> Option<String> {
        match &self.save {
            Some(Some(dir)) => Some(dir.clone()),
            Some(None) => {
                let base = self.scene.strip_suffix(".toml").unwrap_or(&self.scene);
                Some(base.to_string())
            }
            None => None,
        }
    }
}
