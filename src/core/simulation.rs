//! Simulation driver: content-agnostic orchestration.
//!
//! The `Simulation` struct holds content (particles, fields, or both)
//! and a dynamics module that knows how to evolve it. The driver
//! handles the time-stepping schedule, snapshot capture, diagnostics,
//! and pipeline integration without knowing what kind of content
//! it's running.

use crate::config::Configuration;
use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::error::HermesError;
use crate::io::observer::Observer;
use crate::io::snapshot::Snapshot;
use crate::physics::cosmology::Cosmology;
use crate::physics::diagnostics::Diagnostics;
use crate::physics::integrator::scale_factor_schedule;

/// Complete state of a simulation.
pub struct Simulation {
    pub config: Configuration,
    pub cosmology: Cosmology,
    pub content: Content,
    pub dynamics: Box<dyn Dynamics>,
    /// Diagnostics recorded at each snapshot interval.
    pub diagnostics_history: Vec<Diagnostics>,
    /// Current step index.
    pub step: usize,
    /// Current scale factor.
    pub scale_factor: f64,
}

impl Simulation {
    /// Construct a simulation from a scene and configuration.
    ///
    /// The scene provides the initial content and the dynamics module.
    pub fn from_scene(
        scene: &dyn crate::scenes::Scene,
        config: Configuration,
        seed: u64,
    ) -> Result<Self, HermesError> {
        config.cosmology.validate()?;
        scene.validate(&config)?;

        let cosmology = config.cosmology.clone();
        let scale_factor = config.time.scale_factor_initial();

        let (content, dynamics) = scene.initialize(&config, &cosmology, seed)?;

        // Compute initial diagnostics (particle-only for now).
        let diagnostics_history = if let Some(particles) = content.particles() {
            // Need a temporary solver for initial diagnostics.
            let grid = crate::physics::grid::Grid::new(
                config.simulation.n_cells(),
                config.simulation.box_length,
            );
            let mut solver = crate::physics::poisson::PoissonSolver::new(&grid);
            let diag =
                Diagnostics::compute(particles, &grid, &cosmology, &mut solver, scale_factor);
            vec![diag]
        } else {
            Vec::new()
        };

        Ok(Self {
            config,
            cosmology,
            content,
            dynamics,
            diagnostics_history,
            step: 0,
            scale_factor,
        })
    }

    /// Construct from config using the default cosmic-web scene.
    ///
    /// Convenience method for tests and backward compatibility.
    pub fn from_config(config: Configuration, seed: u64) -> Result<Self, HermesError> {
        let scene = crate::scenes::cosmic_web::CosmicWeb;

        Self::from_scene(&scene, config, seed)
    }

    /// Resume a particle simulation from an existing state.
    ///
    /// Takes pre-existing particles and a starting scale factor (from a
    /// loaded snapshot) and builds the simulation for the remaining time
    /// range defined in the config.
    pub fn resume_particles(
        particles: crate::physics::particles::Particles,
        scale_factor: f64,
        start_step: usize,
        config: Configuration,
    ) -> Result<Self, HermesError> {
        config.cosmology.validate()?;

        let cosmology = config.cosmology.clone();
        let grid = crate::physics::grid::Grid::new(
            config.simulation.n_cells(),
            config.simulation.box_length,
        );

        let content = Content::Particles(particles);
        let dynamics: Box<dyn Dynamics> =
            Box::new(crate::core::pm_dynamics::ParticleMeshDynamics::new(grid));

        Ok(Self {
            config,
            cosmology,
            content,
            dynamics,
            diagnostics_history: Vec::new(),
            step: start_step,
            scale_factor,
        })
    }

    /// Run the simulation without a progress callback.
    pub fn run(&mut self, observers: &mut [Box<dyn Observer>]) -> Result<usize, HermesError> {
        self.run_with_progress(observers, |_, _| {})
    }

    /// Run the full simulation with a per-step progress callback.
    pub fn run_with_progress(
        &mut self,
        observers: &mut [Box<dyn Observer>],
        on_step: impl Fn(usize, f64),
    ) -> Result<usize, HermesError> {
        let schedule = scale_factor_schedule(
            self.config.time.scale_factor_initial(),
            self.config.time.scale_factor_final(),
            self.config.time.n_steps,
            &self.config.time.scale_factor_stepping,
        );

        // Notify observers with initial snapshot.
        {
            let initial_snapshot =
                Snapshot::capture_from_content(&self.content, 0, self.scale_factor);
            for observer in observers.iter_mut() {
                observer.on_snapshot(&initial_snapshot);
            }
        }

        let has_observers = !observers.is_empty();

        for n in 0..self.config.time.n_steps {
            let scale_factor_prev = schedule[n];
            let scale_factor_next = schedule[n + 1];

            self.dynamics.step(
                &mut self.content,
                &self.cosmology,
                scale_factor_prev,
                scale_factor_next,
            )?;

            self.step = n + 1;
            self.scale_factor = scale_factor_next;

            on_step(self.step, self.scale_factor);

            let is_final = n + 1 == self.config.time.n_steps;

            // Lightweight snapshot for observers every step.
            if has_observers {
                let snapshot =
                    Snapshot::capture_from_content(&self.content, self.step, self.scale_factor);
                for observer in observers.iter_mut() {
                    observer.on_snapshot(&snapshot);
                }
            }

            // Full diagnostics at wider interval.
            let is_diagnostic_step = self
                .step
                .is_multiple_of(self.config.output.diagnostic_interval)
                || is_final;

            if is_diagnostic_step && let Some(particles) = self.content.particles() {
                // Need solver access for potential energy.
                // The PM dynamics owns the solver, but diagnostics need it too.
                // For now, create a temporary solver. This is wasteful but correct.
                // TODO: expose solver through dynamics trait or compute diagnostics
                // differently for field content.
                let grid = crate::physics::grid::Grid::new(
                    self.config.simulation.n_cells(),
                    self.config.simulation.box_length,
                );
                let mut solver = crate::physics::poisson::PoissonSolver::new(&grid);
                let diagnostics = Diagnostics::compute(
                    particles,
                    &grid,
                    &self.cosmology,
                    &mut solver,
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
    pub fn run_into_pipeline(
        &mut self,
        sender: &crate::run::pipeline::SnapshotSender,
        on_step: impl Fn(usize, f64),
    ) -> Result<usize, HermesError> {
        let schedule = scale_factor_schedule(
            self.config.time.scale_factor_initial(),
            self.config.time.scale_factor_final(),
            self.config.time.n_steps,
            &self.config.time.scale_factor_stepping,
        );

        // Send initial snapshot (skip for resumed runs to avoid duplicating the last frame).
        let step_offset = self.step;
        if step_offset == 0 {
            let initial = std::sync::Arc::new(Snapshot::capture_from_content(
                &self.content,
                0,
                self.scale_factor,
            ));
            sender.send(initial);
        }

        for n in 0..self.config.time.n_steps {
            let scale_factor_prev = schedule[n];
            let scale_factor_next = schedule[n + 1];

            self.dynamics.step(
                &mut self.content,
                &self.cosmology,
                scale_factor_prev,
                scale_factor_next,
            )?;

            self.step = step_offset + n + 1;
            self.scale_factor = scale_factor_next;

            on_step(self.step, self.scale_factor);

            // Snapshot at write_interval or final step.
            let is_final = n + 1 == self.config.time.n_steps;
            let is_write_step =
                self.step.is_multiple_of(self.config.output.write_interval) || is_final;

            if is_write_step {
                let snapshot = std::sync::Arc::new(Snapshot::capture_from_content(
                    &self.content,
                    self.step,
                    self.scale_factor,
                ));
                sender.send(snapshot);
            }

            // Full diagnostics at wider interval.
            let is_diagnostic_step = self
                .step
                .is_multiple_of(self.config.output.diagnostic_interval)
                || is_final;

            if is_diagnostic_step && let Some(particles) = self.content.particles() {
                let grid = crate::physics::grid::Grid::new(
                    self.config.simulation.n_cells(),
                    self.config.simulation.box_length,
                );
                let mut solver = crate::physics::poisson::PoissonSolver::new(&grid);
                let diagnostics = Diagnostics::compute(
                    particles,
                    &grid,
                    &self.cosmology,
                    &mut solver,
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
