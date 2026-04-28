//! Galaxy group initial conditions: multiple colliding dark matter halos.
//!
//! Places N halos with NFW density profiles at positions within
//! the central region of the box, with infall velocities directed
//! toward the group center of mass.

use std::f64::consts::PI;

use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::algebra::vector_from_components;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// Configuration for a single halo in the group.
pub struct HaloConfig {
    /// Fraction of total mass in this halo.
    pub mass_fraction: f64,
    /// NFW concentration parameter.
    pub concentration: f64,
}

/// Default group: 3 halos with mass ratio 1.0 : 0.6 : 0.3.
pub fn default_halo_configs() -> Vec<HaloConfig> {
    vec![
        HaloConfig {
            mass_fraction: 1.0 / 1.9,
            concentration: 8.0,
        },
        HaloConfig {
            mass_fraction: 0.6 / 1.9,
            concentration: 10.0,
        },
        HaloConfig {
            mass_fraction: 0.3 / 1.9,
            concentration: 14.0,
        },
    ]
}

/// Generate initial conditions for multiple colliding halos.
pub fn colliding_halos_init(
    n_per_side: usize,
    grid: &Grid,
    cosmology: &Cosmology,
    _scale_factor_initial: f64,
    seed: u64,
) -> Result<Particles, HermesError> {
    let halo_configs = default_halo_configs();
    let n_halos = halo_configs.len();

    let n_total = n_per_side * n_per_side * n_per_side;
    let box_length = grid.box_length;
    let box_center = box_length / 2.0;

    let density_mean = cosmology.density_matter();
    let mass_total = density_mean * grid.box_volume();
    let mass_particle = mass_total / n_total as f64;

    let g = crate::physics::constants::G;

    // Distribute particles among halos proportional to mass fractions.
    let mut n_per_halo: Vec<usize> = halo_configs
        .iter()
        .map(|h| (h.mass_fraction * n_total as f64) as usize)
        .collect();

    // Assign remaining particles to the largest halo.
    let assigned: usize = n_per_halo.iter().sum();
    n_per_halo[0] += n_total - assigned;

    // Halo masses.
    let halo_masses: Vec<f64> = n_per_halo
        .iter()
        .map(|&n| n as f64 * mass_particle)
        .collect();

    // Place halos on a circle in the central region.
    let spread_radius = box_length * 0.15;

    let mut rng = ChaCha20Rng::seed_from_u64(seed);

    let mut centers = Vec::with_capacity(n_halos);
    for k in 0..n_halos {
        let angle = 2.0 * PI * k as f64 / n_halos as f64;
        let perturbation_x: f64 = rng.random_range(-0.1..0.1) * spread_radius;
        let perturbation_y: f64 = rng.random_range(-0.1..0.1) * spread_radius;
        let perturbation_z: f64 = rng.random_range(-0.2..0.2) * spread_radius;

        centers.push([
            box_center + spread_radius * angle.cos() + perturbation_x,
            box_center + spread_radius * angle.sin() + perturbation_y,
            box_center + perturbation_z,
        ]);
    }

    // Center of mass.
    let total_mass: f64 = halo_masses.iter().sum();
    let com = [
        centers
            .iter()
            .zip(halo_masses.iter())
            .map(|(c, &m)| c[0] * m)
            .sum::<f64>()
            / total_mass,
        centers
            .iter()
            .zip(halo_masses.iter())
            .map(|(c, &m)| c[1] * m)
            .sum::<f64>()
            / total_mass,
        centers
            .iter()
            .zip(halo_masses.iter())
            .map(|(c, &m)| c[2] * m)
            .sum::<f64>()
            / total_mass,
    ];

    // Virial radii: r_vir = (3M / (4pi * 200 * rho_c))^(1/3)
    let density_critical = cosmology.density_critical(1.0);
    let virial_radii: Vec<f64> = halo_masses
        .iter()
        .map(|&m| (3.0 * m / (4.0 * PI * 200.0 * density_critical)).powf(1.0 / 3.0))
        .collect();

    let mut particles = Particles::zeros(n_total, mass_particle);
    let mut particle_index = 0;

    for k in 0..n_halos {
        let radius_virial = virial_radii[k];
        let scale_radius = radius_virial / halo_configs[k].concentration;

        // Infall velocity toward center of mass.
        let dx = com[0] - centers[k][0];
        let dy = com[1] - centers[k][1];
        let dz = com[2] - centers[k][2];
        let distance = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0);
        let infall_speed = (g * total_mass / distance).sqrt() * 0.4;

        let bulk_velocity = [
            infall_speed * dx / distance,
            infall_speed * dy / distance,
            infall_speed * dz / distance,
        ];

        let velocity_dispersion = (g * halo_masses[k] / radius_virial).sqrt() * 0.3;

        sample_nfw_halo(
            &mut particles,
            particle_index,
            n_per_halo[k],
            &centers[k],
            &bulk_velocity,
            radius_virial,
            scale_radius,
            mass_particle,
            velocity_dispersion,
            &mut rng,
        );

        particle_index += n_per_halo[k];
    }

    particles.wrap_positions(grid);

    Ok(particles)
}

/// Sample particles from an NFW density profile using rejection sampling.
#[allow(clippy::too_many_arguments)]
fn sample_nfw_halo(
    particles: &mut Particles,
    start_index: usize,
    n_particles: usize,
    center: &[f64; 3],
    bulk_velocity: &[f64; 3],
    radius_virial: f64,
    scale_radius: f64,
    mass_particle: f64,
    velocity_dispersion: f64,
    rng: &mut ChaCha20Rng,
) {
    let mut placed = 0;

    while placed < n_particles {
        let x: f64 = rng.random_range(-1.0..1.0);
        let y: f64 = rng.random_range(-1.0..1.0);
        let z: f64 = rng.random_range(-1.0..1.0);
        let r2 = x * x + y * y + z * z;

        if !(1e-6..=1.0).contains(&r2) {
            continue;
        }

        let r = r2.sqrt() * radius_virial;
        let s = r / scale_radius;

        let density = 1.0 / (s * (1.0 + s) * (1.0 + s));

        let s_min = 0.01;
        let density_max = 1.0 / (s_min * (1.0 + s_min) * (1.0 + s_min));
        let accept_probability = density / density_max;

        let u: f64 = rng.random_range(0.0..1.0);
        if u > accept_probability {
            continue;
        }

        let p = start_index + placed;

        let position = vector_from_components(
            center[0] + x * radius_virial,
            center[1] + y * radius_virial,
            center[2] + z * radius_virial,
        );
        particles.set_position(p, &position);

        let vx: f64 = rng.random_range(-1.0..1.0);
        let vy: f64 = rng.random_range(-1.0..1.0);
        let vz: f64 = rng.random_range(-1.0..1.0);

        let momentum = vector_from_components(
            mass_particle * (bulk_velocity[0] + velocity_dispersion * vx),
            mass_particle * (bulk_velocity[1] + velocity_dispersion * vy),
            mass_particle * (bulk_velocity[2] + velocity_dispersion * vz),
        );
        particles.set_momentum(p, &momentum);

        placed += 1;
    }
}
