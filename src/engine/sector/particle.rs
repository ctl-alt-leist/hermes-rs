/// Particle sector: discrete phase-space degrees of freedom.
///
/// Advances a named particle species under:
///   T-flow: position drift (x → x + p/m × dt)
///   V-flow: momentum kick from interpolated gravitational force
///
/// The force field is produced by the GravitySolver and delivered
/// through the Potential struct. CIC interpolation maps the grid
/// force to individual particle positions.
use crate::engine::sector::{Potential, Sector};
use crate::engine::state::SimulationState;
use crate::error::HermesError;
use crate::physics::cic::interpolate_force;

/// Particle sector for a named particle species.
pub struct ParticleSector {
    name: String,
}

impl ParticleSector {
    /// Create a particle sector for the named species.
    pub fn new(name: String) -> Self {
        Self { name }
    }
}

impl Sector for ParticleSector {
    fn name(&self) -> &str {
        &self.name
    }

    fn kinetic_flow(
        &mut self,
        state: &mut SimulationState,
        _scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError> {
        let particles = state.particles.get_mut(&self.name).ok_or_else(|| {
            HermesError::Config(format!("particles '{}' not found in state", self.name))
        })?;

        let mass_inv = 1.0 / particles.mass_particle;

        for p in 0..particles.count() {
            let position = particles.position_of(p);
            let momentum = particles.momentum_of(p);
            let displacement = &momentum * (mass_inv * dt);
            let position_new = &position + &displacement;
            particles.set_position(p, &position_new);
        }

        // Wrap positions into the periodic box.
        let box_length = state.grid.box_length;
        for p in 0..particles.count() {
            let pos = particles.position_of(p);
            let wrapped = crate::algebra::vector_from_components(
                ((pos.component(&[1]) % box_length) + box_length) % box_length,
                ((pos.component(&[2]) % box_length) + box_length) % box_length,
                ((pos.component(&[3]) % box_length) + box_length) % box_length,
            );
            particles.set_position(p, &wrapped);
        }

        Ok(())
    }

    fn potential_flow(
        &mut self,
        state: &mut SimulationState,
        potential: &Potential,
        _scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError> {
        let force = potential.force.as_ref().ok_or_else(|| {
            HermesError::Config("particle sector requires force field in potential".to_string())
        })?;

        let particles = state.particles.get_mut(&self.name).ok_or_else(|| {
            HermesError::Config(format!("particles '{}' not found in state", self.name))
        })?;

        let forces = interpolate_force(force, particles, &state.grid);

        for p in 0..particles.count() {
            let momentum = particles.momentum_of(p);
            let f = forces.force_on(p);
            let momentum_new = &momentum + &(&f * dt);
            particles.set_momentum(p, &momentum_new);
        }

        Ok(())
    }

    fn deposit_density(
        &self,
        state: &SimulationState,
    ) -> Result<morphis::field::Field<3>, HermesError> {
        let particles = state.particles.get(&self.name).ok_or_else(|| {
            HermesError::Config(format!("particles '{}' not found in state", self.name))
        })?;

        // CIC deposit onto hermes grid, then convert to morphis Field.
        let hermes_density = crate::physics::cic::assign_density(particles, &state.grid);

        // Convert hermes ScalarField → morphis Field<3>.
        // The callback receives physical positions; convert to grid indices.
        let cell_length = state.grid.cell_length;
        let n = state.grid.n_cells;
        let morphis_density = morphis::field::Field::scalar_field(
            &state.morphis_grid,
            morphis::metric::euclidean::<3>(),
            |pos| {
                let m0 = ((pos[0] / cell_length) as usize).min(n - 1);
                let m1 = ((pos[1] / cell_length) as usize).min(n - 1);
                let m2 = ((pos[2] / cell_length) as usize).min(n - 1);
                hermes_density.data[[m0, m1, m2]]
            },
        );

        Ok(morphis_density)
    }
}
