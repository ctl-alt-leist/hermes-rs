//! Symplectic kick-drift-kick leapfrog with cosmological step factors.
//!
//! The integrator advances the particle state through the gravitational
//! force chain using a time-symmetric KDK scheme. Cosmological expansion
//! enters through the kick and drift factors, which are integrals of
//! powers of the scale factor over the Hubble parameter.
//!
//! A single step from scale factor a_n to a_{n+1} with midpoint a_{n+1/2}:
//!
//! ```text
//! 1. Half-kick:  p → p + F × kick_factor(a_n, a_{n+1/2})
//! 2. Full drift: x → x + (p / m) × drift_factor(a_n, a_{n+1})
//! 3. Recompute force at new positions
//! 4. Half-kick:  p → p + F × kick_factor(a_{n+1/2}, a_{n+1})
//! ```

use crate::error::HermesError;
use crate::physics::cic::{ParticleForces, assign_density, interpolate_force};
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;
use crate::physics::poisson::PoissonSolver;

/// Apply a momentum kick to all particles: p → p + F × kick_factor.
///
/// The kick is a morphis vector operation: each particle's momentum
/// gains a scaled copy of its force vector.
pub fn kick(particles: &mut Particles, forces: &ParticleForces, kick_factor: f64) {
    for p in 0..particles.count() {
        let momentum = particles.momentum_of(p);
        let force = forces.force_on(p);
        let momentum_new = &momentum + &(&force * kick_factor);
        particles.set_momentum(p, &momentum_new);
    }
}

/// Apply a position drift to all particles: x → x + (p / m) × drift_factor.
///
/// Positions are wrapped into the periodic box after the drift.
pub fn drift(particles: &mut Particles, drift_factor: f64, grid: &Grid) {
    let mass_inv = 1.0 / particles.mass_particle;

    for p in 0..particles.count() {
        let position = particles.position_of(p);
        let momentum = particles.momentum_of(p);
        let displacement = &momentum * (mass_inv * drift_factor);
        let position_new = &position + &displacement;
        particles.set_position(p, &position_new);
    }

    particles.wrap_positions(grid);
}

/// Compute the gravitational force on all particles via the CIC → Poisson → CIC chain.
pub fn compute_forces(
    particles: &Particles,
    solver: &mut PoissonSolver,
    grid: &Grid,
    cosmology: &Cosmology,
    scale_factor: f64,
) -> ParticleForces {
    let density = assign_density(particles, grid);
    let density_mean = cosmology.density_matter();

    // Compute overdensity δ = ρ/ρ̄ - 1.
    let mut overdensity = density;
    overdensity.data /= density_mean;
    overdensity.data -= 1.0;

    let force_field = solver.solve(&overdensity, density_mean, scale_factor);

    interpolate_force(&force_field, particles, grid)
}

#[allow(clippy::too_many_arguments)]
/// Execute one full KDK leapfrog step.
///
/// Advances particles from scale factor `scale_factor_prev` to `scale_factor_next`
/// through the midpoint `scale_factor_mid`. If `forces_prev` is provided, it
/// is used for the opening half-kick (avoiding a redundant force evaluation
/// from the previous step's closing half-kick). Returns the force at the
/// new positions for reuse in the next step.
pub fn step_kdk(
    particles: &mut Particles,
    solver: &mut PoissonSolver,
    grid: &Grid,
    cosmology: &Cosmology,
    scale_factor_prev: f64,
    scale_factor_mid: f64,
    scale_factor_next: f64,
    forces_prev: Option<&ParticleForces>,
) -> Result<ParticleForces, HermesError> {
    // Opening half-kick: p^n → p^{n+1/2}
    let kick_factor_open = cosmology.kick_factor(scale_factor_prev, scale_factor_mid);

    let forces_open = match forces_prev {
        Some(f) => f.clone(),
        None => compute_forces(particles, solver, grid, cosmology, scale_factor_prev),
    };

    kick(particles, &forces_open, kick_factor_open);

    // Full drift: x^n → x^{n+1}
    let drift_fac = cosmology.drift_factor(scale_factor_prev, scale_factor_next);
    drift(particles, drift_fac, grid);

    // Recompute force at new positions
    let forces_new = compute_forces(particles, solver, grid, cosmology, scale_factor_next);

    // Closing half-kick: p^{n+1/2} → p^{n+1}
    let kick_factor_close = cosmology.kick_factor(scale_factor_mid, scale_factor_next);
    kick(particles, &forces_new, kick_factor_close);

    Ok(forces_new)
}

/// Generate a schedule of scale factors for time stepping.
///
/// Returns `n_steps + 1` scale factors from `scale_factor_initial` to `scale_factor_final`,
/// spaced logarithmically in a (default) or linearly.
pub fn scale_factor_schedule(
    scale_factor_initial: f64,
    scale_factor_final: f64,
    n_steps: usize,
    stepping: &str,
) -> Vec<f64> {
    let mut schedule = Vec::with_capacity(n_steps + 1);

    match stepping {
        "log_a" => {
            let log_start = scale_factor_initial.ln();
            let log_end = scale_factor_final.ln();
            let d_log = (log_end - log_start) / n_steps as f64;

            for n in 0..=n_steps {
                schedule.push((log_start + n as f64 * d_log).exp());
            }
        }
        "linear_a" => {
            let da = (scale_factor_final - scale_factor_initial) / n_steps as f64;

            for n in 0..=n_steps {
                schedule.push(scale_factor_initial + n as f64 * da);
            }
        }
        _ => {
            // Default to log_a
            let log_start = scale_factor_initial.ln();
            let log_end = scale_factor_final.ln();
            let d_log = (log_end - log_start) / n_steps as f64;

            for n in 0..=n_steps {
                schedule.push((log_start + n as f64 * d_log).exp());
            }
        }
    }

    schedule
}

/// Midpoint scale factor between two schedule entries.
pub fn midpoint(scale_factor_prev: f64, scale_factor_next: f64) -> f64 {
    // Geometric mean for logarithmic stepping.
    (scale_factor_prev * scale_factor_next).sqrt()
}
