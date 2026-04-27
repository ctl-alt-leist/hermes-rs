//! Simulation driver tying initialization, time stepping, and diagnostics.
//!
//! The `Simulation` struct owns all state needed for a complete PM run:
//! particles, solver, grid, cosmology, and the diagnostics history.
//! Construction from a `Configuration` initializes everything; the `run`
//! method advances through the time-step schedule and notifies observers
//! at snapshot intervals.

use crate::config::Configuration;
use crate::error::HermesError;
use crate::io::observer::Observer;
use crate::io::snapshot::Snapshot;
use crate::physics::cosmology::Cosmology;
use crate::physics::diagnostics::Diagnostics;
use crate::physics::grid::Grid;
use crate::physics::initial::zeldovich_init;
use crate::physics::integrator::{midpoint, scale_factor_schedule, step_kdk};
use crate::physics::particles::Particles;
use crate::physics::poisson::PoissonSolver;

/// Complete state of a particle-mesh cosmological simulation.
pub struct Simulation {
    pub config: Configuration,
    pub cosmology: Cosmology,
    pub grid: Grid,
    pub particles: Particles,
    pub solver: PoissonSolver,
    /// Diagnostics recorded at each snapshot interval.
    pub diagnostics_history: Vec<Diagnostics>,
    /// Current step index.
    pub step: usize,
    /// Current scale factor.
    pub scale_factor: f64,
}

impl Simulation {
    /// Construct a simulation from a configuration.
    ///
    /// Initializes the grid, Poisson solver, and particles via Zel'dovich
    /// approximation. Records initial diagnostics.
    pub fn from_config(config: Configuration, seed: u64) -> Result<Self, HermesError> {
        config.cosmology.validate()?;

        let grid = Grid::new(config.simulation.n_cells, config.simulation.box_length);
        let mut solver = PoissonSolver::new(&grid);
        let cosmology = config.cosmology.clone();
        let scale_factor = config.time.scale_factor_initial;

        let particles = zeldovich_init(
            config.simulation.n_particles,
            &grid,
            &cosmology,
            scale_factor,
            seed,
        )?;

        let initial_diagnostics =
            Diagnostics::compute(&particles, &grid, &cosmology, &mut solver, scale_factor);

        Ok(Self {
            config,
            cosmology,
            grid,
            particles,
            solver,
            diagnostics_history: vec![initial_diagnostics],
            step: 0,
            scale_factor,
        })
    }

    /// Run the full simulation, notifying observers at snapshot intervals.
    ///
    /// Observers receive a `Snapshot` (morphis-native positions and momenta)
    /// at each snapshot interval. Multiple observers can run simultaneously —
    /// e.g. a `FileObserver` writing to disk and a `MemoryObserver` collecting
    /// in memory.
    /// Run the simulation without a progress callback.
    pub fn run(&mut self, observers: &mut [Box<dyn Observer>]) -> Result<usize, HermesError> {
        self.run_with_progress(observers, |_, _| {})
    }

    /// Run the full simulation with a per-step progress callback.
    ///
    /// The callback receives `(step, scale_factor)` after each step.
    pub fn run_with_progress(
        &mut self,
        observers: &mut [Box<dyn Observer>],
        on_step: impl Fn(usize, f64),
    ) -> Result<usize, HermesError> {
        let schedule = scale_factor_schedule(
            self.config.time.scale_factor_initial,
            self.config.time.scale_factor_final,
            self.config.time.n_steps,
            &self.config.time.stepping,
        );

        // Notify observers with initial snapshot (lightweight — no Poisson).
        let initial_snapshot = Snapshot::capture_lightweight(&self.particles, 0, self.scale_factor);
        for observer in observers.iter_mut() {
            observer.on_snapshot(&initial_snapshot);
        }

        let mut forces_prev = None;
        let has_observers = !observers.is_empty();

        for n in 0..self.config.time.n_steps {
            let scale_factor_prev = schedule[n];
            let scale_factor_next = schedule[n + 1];
            let scale_factor_mid = midpoint(scale_factor_prev, scale_factor_next);

            let forces = step_kdk(
                &mut self.particles,
                &mut self.solver,
                &self.grid,
                &self.cosmology,
                scale_factor_prev,
                scale_factor_mid,
                scale_factor_next,
                forces_prev.as_ref(),
            )?;

            self.step = n + 1;
            self.scale_factor = scale_factor_next;
            forces_prev = Some(forces);

            on_step(self.step, self.scale_factor);

            // Lightweight snapshot for observers every step.
            if has_observers {
                let snapshot =
                    Snapshot::capture_lightweight(&self.particles, self.step, self.scale_factor);
                for observer in observers.iter_mut() {
                    observer.on_snapshot(&snapshot);
                }
            }

            // Full diagnostics (expensive Poisson solve) only at the wider interval.
            let is_diagnostic_step = self
                .step
                .is_multiple_of(self.config.output.snapshot_interval)
                || self.step == self.config.time.n_steps;

            if is_diagnostic_step {
                let diagnostics = Diagnostics::compute(
                    &self.particles,
                    &self.grid,
                    &self.cosmology,
                    &mut self.solver,
                    self.scale_factor,
                );
                self.diagnostics_history.push(diagnostics);
            }
        }

        for observer in observers.iter_mut() {
            observer.on_complete();
        }

        Ok(self.step)
    }

    /// Run the simulation, sending snapshots into a pipeline channel.
    ///
    /// The simulation sends `Arc<Snapshot>` through the `SnapshotSender`;
    /// downstream consumers (disk writer, viewer) receive them via the
    /// router. The simulation never blocks on consumers.
    pub fn run_into_pipeline(
        &mut self,
        sender: &crate::run::pipeline::SnapshotSender,
        on_step: impl Fn(usize, f64),
    ) -> Result<usize, HermesError> {
        let schedule = scale_factor_schedule(
            self.config.time.scale_factor_initial,
            self.config.time.scale_factor_final,
            self.config.time.n_steps,
            &self.config.time.stepping,
        );

        // Send initial snapshot.
        let initial = std::sync::Arc::new(Snapshot::capture_lightweight(
            &self.particles,
            0,
            self.scale_factor,
        ));
        sender.send(initial);

        let mut forces_prev = None;

        for n in 0..self.config.time.n_steps {
            let scale_factor_prev = schedule[n];
            let scale_factor_next = schedule[n + 1];
            let scale_factor_mid = midpoint(scale_factor_prev, scale_factor_next);

            let forces = step_kdk(
                &mut self.particles,
                &mut self.solver,
                &self.grid,
                &self.cosmology,
                scale_factor_prev,
                scale_factor_mid,
                scale_factor_next,
                forces_prev.as_ref(),
            )?;

            self.step = n + 1;
            self.scale_factor = scale_factor_next;
            forces_prev = Some(forces);

            on_step(self.step, self.scale_factor);

            // Lightweight snapshot wrapped in Arc for zero-copy fan-out.
            let snapshot = std::sync::Arc::new(Snapshot::capture_lightweight(
                &self.particles,
                self.step,
                self.scale_factor,
            ));
            sender.send(snapshot);

            // Full diagnostics at wider interval (stays on sim thread).
            let is_diagnostic_step = self
                .step
                .is_multiple_of(self.config.output.snapshot_interval)
                || self.step == self.config.time.n_steps;

            if is_diagnostic_step {
                let diagnostics = Diagnostics::compute(
                    &self.particles,
                    &self.grid,
                    &self.cosmology,
                    &mut self.solver,
                    self.scale_factor,
                );
                self.diagnostics_history.push(diagnostics);
            }
        }

        sender.done();

        Ok(self.step)
    }

    /// Latest recorded diagnostics.
    pub fn latest_diagnostics(&self) -> Option<&Diagnostics> {
        self.diagnostics_history.last()
    }
}
