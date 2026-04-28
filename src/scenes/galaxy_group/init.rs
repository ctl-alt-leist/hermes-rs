//! Galaxy group initial conditions: colliding dark matter halos.
//!
//! Places two spherical halos with NFW-like radial density profiles,
//! separated in space with a relative bulk velocity. The particles
//! are distributed within each halo using rejection sampling from
//! the density profile.

use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::algebra::vector_from_components;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// Generate initial conditions for two colliding halos.
///
/// Each halo has an NFW-like density profile truncated at the virial
/// radius. The halos are placed symmetrically about the box center
/// with a relative approach velocity.
pub fn colliding_halos_init(
    n_per_side: usize,
    grid: &Grid,
    cosmology: &Cosmology,
    _scale_factor_initial: f64,
    seed: u64,
) -> Result<Particles, HermesError> {
    let n_total = n_per_side * n_per_side * n_per_side;
    let box_length = grid.box_length;
    let box_center = box_length / 2.0;

    // Total mass in the box from cosmology.
    let density_mean = cosmology.density_matter();
    let mass_total = density_mean * grid.box_volume();
    let mass_particle = mass_total / n_total as f64;

    // Two equal-mass halos, each getting half the particles.
    let n_per_halo = n_total / 2;
    let mass_halo = mass_particle * n_per_halo as f64;

    // Halo parameters.
    let radius_virial = box_length * 0.12; // virial radius ~ 12% of box
    let concentration = 10.0; // typical NFW concentration
    let scale_radius = radius_virial / concentration;

    // Separation: halos placed 40% of box apart, centered in the box.
    let separation = box_length * 0.4;

    // Approach velocity: roughly the circular velocity at the virial radius.
    // v_circ = sqrt(G M / r_vir)
    let g = crate::physics::constants::G;
    let velocity_approach = (g * mass_halo / radius_virial).sqrt() * 0.5;

    let mut particles = Particles::zeros(n_total, mass_particle);
    let mut rng = ChaCha20Rng::seed_from_u64(seed);

    // Halo 1: left of center, moving right.
    let center_1 = [box_center - separation / 2.0, box_center, box_center];
    let velocity_1 = [velocity_approach, 0.0, 0.0];

    sample_nfw_halo(
        &mut particles,
        0,
        n_per_halo,
        &center_1,
        &velocity_1,
        radius_virial,
        scale_radius,
        mass_particle,
        &mut rng,
    );

    // Halo 2: right of center, moving left.
    let center_2 = [box_center + separation / 2.0, box_center, box_center];
    let velocity_2 = [-velocity_approach, 0.0, 0.0];

    sample_nfw_halo(
        &mut particles,
        n_per_halo,
        n_per_halo,
        &center_2,
        &velocity_2,
        radius_virial,
        scale_radius,
        mass_particle,
        &mut rng,
    );

    particles.wrap_positions(grid);

    Ok(particles)
}

#[allow(clippy::too_many_arguments)]
/// Sample particles from an NFW density profile using rejection sampling.
///
/// NFW profile: rho(r) = rho_0 / [(r/r_s)(1 + r/r_s)^2]
///
/// Particles are placed with the given bulk velocity plus a small
/// isotropic velocity dispersion for stability.
fn sample_nfw_halo(
    particles: &mut Particles,
    start_index: usize,
    n_particles: usize,
    center: &[f64; 3],
    bulk_velocity: &[f64; 3],
    radius_virial: f64,
    scale_radius: f64,
    mass_particle: f64,
    rng: &mut ChaCha20Rng,
) {
    let g = crate::physics::constants::G;

    // Velocity dispersion: approximate as fraction of circular velocity.
    let mass_enclosed = mass_particle * n_particles as f64;
    let velocity_dispersion = (g * mass_enclosed / radius_virial).sqrt() * 0.3;

    // NFW profile peaks at r = 0. The maximum density is at the smallest
    // radius we sample. We use rejection sampling in spherical shells.
    let mut placed = 0;

    while placed < n_particles {
        // Sample uniformly in the sphere of radius r_vir.
        let x: f64 = rng.random_range(-1.0..1.0);
        let y: f64 = rng.random_range(-1.0..1.0);
        let z: f64 = rng.random_range(-1.0..1.0);
        let r2 = x * x + y * y + z * z;

        if !(1e-6..=1.0).contains(&r2) {
            continue; // outside unit sphere or too close to center
        }

        let r = r2.sqrt() * radius_virial;
        let s = r / scale_radius;

        // NFW density (unnormalized): 1 / [s * (1 + s)^2]
        let density = 1.0 / (s * (1.0 + s) * (1.0 + s));

        // Rejection: accept with probability proportional to density.
        // Max density is at r → 0 (s → 0), which diverges. We cap at
        // s_min corresponding to the inner resolution limit.
        let s_min = 0.01;
        let density_max = 1.0 / (s_min * (1.0 + s_min) * (1.0 + s_min));
        let accept_probability = density / density_max;

        let u: f64 = rng.random_range(0.0..1.0);
        if u > accept_probability {
            continue;
        }

        let p = start_index + placed;

        // Position: center + r * (x, y, z) / |xyz|
        let position = vector_from_components(
            center[0] + x * radius_virial,
            center[1] + y * radius_virial,
            center[2] + z * radius_virial,
        );
        particles.set_position(p, &position);

        // Momentum: bulk velocity + isotropic dispersion.
        let vx: f64 = rng.random_range(-1.0..1.0);
        let vy: f64 = rng.random_range(-1.0..1.0);
        let vz: f64 = rng.random_range(-1.0..1.0);

        // p = m * a^2 * dx/dt. At a = 1 (or whatever a_init), p = m * v.
        let momentum = vector_from_components(
            mass_particle * (bulk_velocity[0] + velocity_dispersion * vx),
            mass_particle * (bulk_velocity[1] + velocity_dispersion * vy),
            mass_particle * (bulk_velocity[2] + velocity_dispersion * vz),
        );
        particles.set_momentum(p, &momentum);

        placed += 1;
    }
}
