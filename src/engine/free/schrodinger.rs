/// Free Schrodinger evolution: kinetic step of the split-step integrator.
///
/// Advances an even-subalgebra field under the free kinetic term
/// -l/(2m a²) ∇²α. This is the dispersive part of the Schrodinger
/// equation — wave packets spread, plane waves rotate in phase space,
/// and the norm is exactly preserved.
use crate::engine::free::FreeEvolution;
use crate::engine::state::FieldEntry;
use crate::error::HermesError;

/// Free Schrodinger kinetic evolution.
///
/// Stateless — the phase rotation coefficients are computed from
/// the field parameters and timestep at each call.
pub struct SchrodingerEvolution;

impl FreeEvolution for SchrodingerEvolution {
    fn step(
        &mut self,
        field: &mut FieldEntry,
        grid: &morphis::grid::Grid<3>,
        scale_factor: f64,
        dt: f64,
    ) -> Result<(), HermesError> {
        crate::core::schrodinger_dynamics::kinetic_step(
            &mut field.data,
            grid,
            field.smoothing_length,
            field.mass,
            scale_factor,
            dt,
        );

        Ok(())
    }
}
