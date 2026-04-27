//! Simulation driver tying initialization, time stepping, and diagnostics.
//!
//! The `Simulation` struct owns all state needed for a complete PM run:
//! particles, solver, grid, cosmology, and the diagnostics history.
//! Construction from a `Configuration` initializes everything; the `run`
//! method advances through the time-step schedule.

use crate::config::Configuration;
use crate::cosmology::Cosmology;
use crate::diagnostics::Diagnostics;
use crate::error::HermesError;
use crate::grid::Grid;
use crate::initial::zeldovich_init;
use crate::integrator::{midpoint, scale_factor_schedule, step_kdk};
use crate::particles::Particles;
use crate::poisson::PoissonSolver;

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

    /// Run the full simulation from initial to final scale factor.
    ///
    /// Returns the number of steps completed.
    pub fn run(&mut self) -> Result<usize, HermesError> {
        let schedule = scale_factor_schedule(
            self.config.time.scale_factor_initial,
            self.config.time.scale_factor_final,
            self.config.time.n_steps,
            &self.config.time.stepping,
        );

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

            // Record diagnostics at snapshot intervals.
            if self
                .step
                .is_multiple_of(self.config.output.snapshot_interval)
                || self.step == self.config.time.n_steps
            {
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

        Ok(self.step)
    }

    /// Latest recorded diagnostics.
    pub fn latest_diagnostics(&self) -> Option<&Diagnostics> {
        self.diagnostics_history.last()
    }
}
