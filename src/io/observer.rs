//! Observer trait and implementations for simulation monitoring.
//!
//! Observers are notified at snapshot intervals during the simulation run.
//! They never modify the simulation — they only receive read-only snapshots.
//!
//! Implementations:
//! - `FileObserver` — writes snapshots to `./data/` as bincode files
//! - `NullObserver` — does nothing (headless batch runs)

use std::path::PathBuf;

use crate::io::snapshot::{Snapshot, save_snapshot};

/// Trait for observing simulation progress.
///
/// Observers receive snapshots at configured intervals and are notified
/// when the simulation completes. Observers must be `Send` to support
/// future threaded viewer implementations.
pub trait Observer: Send {
    /// Called when a snapshot is available.
    fn on_snapshot(&mut self, snapshot: &Snapshot);

    /// Called when the simulation run completes.
    fn on_complete(&mut self);
}

/// Observer that writes snapshots to disk as bincode files.
///
/// Files are written to `{directory}/snapshot_{step:05}.bin`.
/// The directory is created automatically if it doesn't exist.
pub struct FileObserver {
    directory: PathBuf,
    n_saved: usize,
}

impl FileObserver {
    /// Create a file observer that writes to the given directory.
    pub fn new(directory: impl Into<PathBuf>) -> Self {
        Self {
            directory: directory.into(),
            n_saved: 0,
        }
    }

    /// Number of snapshots saved so far.
    pub fn n_saved(&self) -> usize {
        self.n_saved
    }
}

impl Observer for FileObserver {
    fn on_snapshot(&mut self, snapshot: &Snapshot) {
        let filename = format!("snapshot_{:05}.bin", snapshot.step);
        let path = self.directory.join(filename);

        match save_snapshot(snapshot, &path) {
            Ok(()) => self.n_saved += 1,
            Err(e) => eprintln!("warning: failed to save snapshot: {e}"),
        }
    }

    fn on_complete(&mut self) {
        if self.n_saved > 0 {
            println!(
                "FileObserver: saved {n} snapshots to {dir}",
                n = self.n_saved,
                dir = self.directory.display()
            );
        }
    }
}

/// Observer that does nothing. For headless batch runs.
pub struct NullObserver;

impl Observer for NullObserver {
    fn on_snapshot(&mut self, _snapshot: &Snapshot) {}
    fn on_complete(&mut self) {}
}

/// Observer that collects snapshots in memory for post-run analysis.
pub struct MemoryObserver {
    snapshots: Vec<Snapshot>,
}

impl MemoryObserver {
    pub fn new() -> Self {
        Self {
            snapshots: Vec::new(),
        }
    }

    /// Access all recorded snapshots.
    pub fn snapshots(&self) -> &[Snapshot] {
        &self.snapshots
    }

    /// Take ownership of the recorded snapshots.
    pub fn into_snapshots(self) -> Vec<Snapshot> {
        self.snapshots
    }
}

impl Default for MemoryObserver {
    fn default() -> Self {
        Self::new()
    }
}

impl Observer for MemoryObserver {
    fn on_snapshot(&mut self, snapshot: &Snapshot) {
        self.snapshots.push(snapshot.clone());
    }

    fn on_complete(&mut self) {}
}
