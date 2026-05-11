/// Mixed dynamics: particles and fields under shared gravity.
///
/// Composes KDK particle stepping with Schrodinger field split-step
/// in a single Strang-symmetric timestep. Both representations source
/// and feel the gravitational potential independently (each through
/// its own Poisson solve path).
///
/// Strang composition: T/2 → V → T/2
///   T = field kinetic step (Fourier-space phase rotation)
///   V = particle KDK step + field potential step (Poisson gravity)
use crate::core::content::Content;
use crate::core::dynamics::Dynamics;
use crate::core::schrodinger_dynamics::kinetic_step;
use crate::engine::coupling::poisson::PoissonGravity;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::integrator::midpoint;

/// Mixed particle + field dynamics.
pub struct MixedDynamics {
    gravity: PoissonGravity,
}

impl MixedDynamics {
    /// Create mixed dynamics with a gravity module.
    pub fn new(gravity: PoissonGravity) -> Self {
        Self { gravity }
    }
}

impl Dynamics for MixedDynamics {
    fn step(
        &mut self,
        content: &mut Content,
        cosmology: &Cosmology,
        scale_factor_prev: f64,
        scale_factor_next: f64,
    ) -> Result<(), HermesError> {
        let fields = content.fields().ok_or_else(|| {
            HermesError::Config("mixed dynamics requires field content".to_string())
        })?;

        let alpha = fields
            .alpha
            .as_ref()
            .ok_or_else(|| HermesError::Config("mixed dynamics requires α field".to_string()))?;

        let scale_factor = (scale_factor_prev + scale_factor_next) / 2.0;
        let dt = (scale_factor_next - scale_factor_prev)
            / (scale_factor * cosmology.hubble_parameter(scale_factor));

        let ell = fields.params.smoothing_length;
        let mass = fields.params.mass_alpha;
        let grid = alpha.grid;
        let scale_factor_mid = midpoint(scale_factor_prev, scale_factor_next);

        // 1. Field kinetic half-step T(dt/2).
        kinetic_step(
            content.fields_mut().unwrap().alpha.as_mut().unwrap(),
            &grid,
            ell,
            mass,
            scale_factor,
            dt / 2.0,
        );

        // 2. Particle full KDK step (kick → drift → kick).
        self.gravity.step_particles_kdk(
            content,
            cosmology,
            scale_factor_prev,
            scale_factor_mid,
            scale_factor_next,
        )?;

        // 3. Field potential full step V(dt).
        self.gravity
            .potential_step_field(content, cosmology, scale_factor, dt)?;

        // 4. Field kinetic half-step T(dt/2).
        kinetic_step(
            content.fields_mut().unwrap().alpha.as_mut().unwrap(),
            &grid,
            ell,
            mass,
            scale_factor,
            dt / 2.0,
        );

        Ok(())
    }
}
