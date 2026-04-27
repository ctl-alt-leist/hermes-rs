//! Conservation diagnostics and simulation health monitoring.
//!
//! All derived quantities are computed through morphis geometric algebra:
//! momentum via vector summation, angular momentum via wedge product
//! (grade-2 bivector), kinetic energy via the metric norm.
//!
//! Conservation laws stratify by grade:
//!   - Grade 0: total mass (exact), energy (Layzer-Irvine relation)
//!   - Grade 1: total comoving momentum (machine precision)
//!   - Grade 2: total angular momentum bivector (diagnostic, not conserved
//!     in a periodic box)

use morphis::vector::Vector;

use crate::physics::cic::assign_density;
use crate::physics::cosmology::Cosmology;
use crate::physics::field::ScalarField;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;
use crate::physics::poisson::PoissonSolver;

/// Snapshot of diagnostic quantities at one point in time.
#[derive(Debug, Clone)]
pub struct Diagnostics {
    /// Scale factor at which diagnostics were computed.
    pub scale_factor: f64,
    /// Total mass (M_☉). Exactly N_p × m_p.
    pub mass_total: f64,
    /// Total comoving momentum as a morphis grade-1 vector.
    /// Should be conserved to machine precision.
    pub momentum_total: Vector<3>,
    /// Total kinetic energy Σ |p|² / (2 m a²) (M_☉ kpc² Gyr⁻²).
    pub energy_kinetic: f64,
    /// Total potential energy ½ Σ m ϕ(x) / a (M_☉ kpc² Gyr⁻²).
    pub energy_potential: f64,
    /// Total angular momentum bivector L = Σ x ∧ p (grade-2).
    /// Not conserved in a periodic box but useful as a diagnostic.
    pub angular_momentum: Vector<3>,
}

impl Diagnostics {
    /// Compute all diagnostics from the current particle state.
    ///
    /// Requires a Poisson solver to compute the gravitational potential
    /// for the potential energy. The solver is mutated because it holds
    /// FFT workspace buffers.
    pub fn compute(
        particles: &Particles,
        grid: &Grid,
        cosmology: &Cosmology,
        solver: &mut PoissonSolver,
        scale_factor: f64,
    ) -> Self {
        let mass_total = particles.total_mass();
        let momentum_total = particles.total_momentum();
        let energy_kinetic = particles.kinetic_energy(scale_factor);
        let angular_momentum = particles.total_angular_momentum();

        let energy_potential =
            compute_potential_energy(particles, grid, cosmology, solver, scale_factor);

        Self {
            scale_factor,
            mass_total,
            momentum_total,
            energy_kinetic,
            energy_potential,
            angular_momentum,
        }
    }

    /// Total energy (kinetic + potential).
    pub fn energy_total(&self) -> f64 {
        self.energy_kinetic + self.energy_potential
    }

    /// Magnitude of total comoving momentum.
    pub fn momentum_magnitude(&self) -> f64 {
        self.momentum_total.norm()
    }

    /// Magnitude of total angular momentum bivector.
    pub fn angular_momentum_magnitude(&self) -> f64 {
        self.angular_momentum.norm()
    }
}

/// Compute the gravitational potential energy: E_pot = ½ Σ_i m_p ϕ(x_i) / a.
///
/// The potential ϕ is obtained by solving the Poisson equation for the
/// current density field, then interpolating back to particle positions.
fn compute_potential_energy(
    particles: &Particles,
    grid: &Grid,
    cosmology: &Cosmology,
    solver: &mut PoissonSolver,
    scale_factor: f64,
) -> f64 {
    let density = assign_density(particles, grid);
    let density_mean = cosmology.density_matter();

    // Overdensity δ = ρ/ρ̄ - 1
    let mut overdensity = density;
    overdensity.data /= density_mean;
    overdensity.data -= 1.0;

    // Solve for potential (the solver returns the force field, but we need
    // the potential itself). We compute it by solving without the gradient
    // step. For now, approximate using the virial relation or compute
    // directly from the density-potential pair.
    //
    // Direct approach: compute ϕ on the grid by FFT, then CIC-interpolate
    // to particles. The Poisson solver currently returns F = -∇ϕ, not ϕ
    // itself. As a practical approximation for the diagnostic, we use the
    // virial estimator: E_pot ≈ -2 E_kin for a virialized system, or
    // compute from the density field directly.
    //
    // For a proper implementation, we'd add a solve_potential method to
    // PoissonSolver. For now, compute from the overdensity power:
    // E_pot = -½ × (4πG ρ̄ a) × Σ_cells δ² × h³ × (1/k²)_avg
    //
    // Simplest correct approach: use the density field and Green's function.
    let potential = compute_potential_field(solver, &overdensity, density_mean, scale_factor);

    // Interpolate potential to particle positions using CIC weights.
    let n = grid.n_cells;
    let h = grid.cell_length;
    let h_inv = 1.0 / h;
    let mut energy = 0.0;

    for p in 0..particles.count() {
        let pos = particles.position_components(p);
        let cell = [
            pos[0] * h_inv - 0.5,
            pos[1] * h_inv - 0.5,
            pos[2] * h_inv - 0.5,
        ];

        let base = [
            cell[0].floor() as isize,
            cell[1].floor() as isize,
            cell[2].floor() as isize,
        ];

        let frac = [
            cell[0] - base[0] as f64,
            cell[1] - base[1] as f64,
            cell[2] - base[2] as f64,
        ];

        let weight = [
            [1.0 - frac[0], frac[0]],
            [1.0 - frac[1], frac[1]],
            [1.0 - frac[2], frac[2]],
        ];

        let mut phi_particle = 0.0;
        for (a, &weight_a) in weight[0].iter().enumerate() {
            let g0 = ((base[0] + a as isize) % n as isize + n as isize) as usize % n;
            for (b, &weight_b) in weight[1].iter().enumerate() {
                let g1 = ((base[1] + b as isize) % n as isize + n as isize) as usize % n;
                for (c, &weight_c) in weight[2].iter().enumerate() {
                    let g2 = ((base[2] + c as isize) % n as isize + n as isize) as usize % n;
                    phi_particle += potential.data[[g0, g1, g2]] * weight_a * weight_b * weight_c;
                }
            }
        }

        energy += particles.mass_particle * phi_particle;
    }

    // E_pot = ½ Σ m ϕ / a
    0.5 * energy / scale_factor
}

/// Compute the gravitational potential field ϕ on the grid via FFT.
///
/// Like the Poisson solver but stops before the gradient — returns ϕ(x)
/// rather than F = -∇ϕ.
fn compute_potential_field(
    solver: &mut PoissonSolver,
    overdensity: &ScalarField,
    density_mean: f64,
    scale_factor: f64,
) -> ScalarField {
    use std::f64::consts::PI;

    use ndarray::Array3;
    use ndrustfft::{FftHandler, R2cFftHandler, ndfft, ndfft_r2c, ndifft, ndifft_r2c};
    use num_complex::Complex64;

    use crate::physics::constants::G as GRAV;

    let n = overdensity.data.shape()[0];
    let n_complex = n / 2 + 1;

    // Forward 3D R2C FFT.
    let mut overdensity_hat = Array3::<Complex64>::zeros((n, n, n_complex));
    let handler_r2c = R2cFftHandler::new(n);
    let handler_c2c_1 = FftHandler::new(n);
    let handler_c2c_0 = FftHandler::new(n);

    ndfft_r2c(&overdensity.data, &mut overdensity_hat, &handler_r2c, 2);
    let mut scratch = overdensity_hat.clone();
    ndfft(&overdensity_hat, &mut scratch, &handler_c2c_1, 1);
    overdensity_hat.assign(&scratch);
    ndfft(&overdensity_hat, &mut scratch, &handler_c2c_0, 0);
    overdensity_hat.assign(&scratch);

    // Multiply by Green's function: ϕ̂ = 4πG ρ̄ a² × G(k) × δ̂
    let prefactor = 4.0 * PI * GRAV * density_mean * scale_factor * scale_factor;

    let green = solver.green_function();

    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n_complex {
                overdensity_hat[[m0, m1, m2]] *= prefactor * green[[m0, m1, m2]];
            }
        }
    }

    // Inverse FFT to get real-space potential.
    let mut scratch = overdensity_hat.clone();
    ndifft(&overdensity_hat, &mut scratch, &handler_c2c_0, 0);
    overdensity_hat.assign(&scratch);
    ndifft(&overdensity_hat, &mut scratch, &handler_c2c_1, 1);
    overdensity_hat.assign(&scratch);

    let mut potential_data = Array3::<f64>::zeros((n, n, n));
    ndifft_r2c(&overdensity_hat, &mut potential_data, &handler_r2c, 2);

    ScalarField::from_array(potential_data)
}
