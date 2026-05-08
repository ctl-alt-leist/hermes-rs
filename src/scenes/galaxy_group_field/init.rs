//! Galaxy group field initial conditions.
//!
//! Places NFW density profiles on the grid with bulk infall velocities
//! encoded as linear phase ramps in the even subalgebra. The density
//! and halo geometry match the particle-mesh galaxy group scene so the
//! two representations start from the same physical state.
//!
//! Each halo contributes a term
//!
//!   α_k(x) = sqrt(ρ_k(x) / m) * exp(I * m * v_k · (x - x_k) / l)
//!
//! and the total field is the coherent sum α = Σ_k α_k. The gradient
//! of the linear phase gives ∇(v_k · (x - x_k)) = v_k exactly, so
//! each halo carries the correct bulk velocity. This is structurally
//! distinct from the Zel'dovich case, which requires a velocity
//! potential because v(x) varies spatially.

use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::metric;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::core::content::FieldParams;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid as HermesGrid;
use crate::scenes::galaxy_group_pm::init::default_halo_configs;

/// Generate colliding NFW halos as a wavefunction field.
///
/// The density is a superposition of NFW profiles placed at the same
/// positions as the particle scene. The velocity is encoded as a
/// linear phase ramp per halo: exp(I * m * v_k · (x - x_k) / l).
pub fn colliding_halos_field(
    grid: &HermesGrid,
    cosmology: &Cosmology,
    params: &FieldParams,
    _scale_factor_initial: f64,
    seed: u64,
) -> EvenField<3> {
    let halo_configs = default_halo_configs();
    let n_halos = halo_configs.len();

    let n = grid.n_cells;
    let box_length = grid.box_length;
    let box_center = box_length / 2.0;
    let cell_length = grid.cell_length;

    let density_mean = cosmology.density_matter();
    let mass_total = density_mean * grid.box_volume();
    let g_const = crate::physics::constants::G;

    let ell = params.smoothing_length;
    let mass = params.mass_alpha;
    let nu = ell / mass;

    // Halo mass fractions -> absolute masses.
    let halo_masses: Vec<f64> = halo_configs
        .iter()
        .map(|h| h.mass_fraction * mass_total)
        .collect();

    // Place halos on a circle (same geometry as PM scene).
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

    // Virial radii.
    let density_critical = cosmology.density_critical(1.0);
    let virial_radii: Vec<f64> = halo_masses
        .iter()
        .map(|&m| (3.0 * m / (4.0 * PI * 200.0 * density_critical)).powf(1.0 / 3.0))
        .collect();

    let scale_radii: Vec<f64> = virial_radii
        .iter()
        .zip(halo_configs.iter())
        .map(|(&r_vir, config)| r_vir / config.concentration)
        .collect();

    // Infall velocities toward center of mass (same formula as PM scene).
    let velocities: Vec<[f64; 3]> = centers
        .iter()
        .map(|center| {
            let dx = com[0] - center[0];
            let dy = com[1] - center[1];
            let dz = com[2] - center[2];
            let distance = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0);
            let infall_speed = (g_const * total_mass / distance).sqrt() * 0.4;

            [
                infall_speed * dx / distance,
                infall_speed * dy / distance,
                infall_speed * dz / distance,
            ]
        })
        .collect();

    // Log regime diagnostics.
    for (k, (vel, &r_vir)) in velocities.iter().zip(virial_radii.iter()).enumerate() {
        let v_mag = (vel[0].powi(2) + vel[1].powi(2) + vel[2].powi(2)).sqrt();
        let nyquist_ratio = mass * v_mag * cell_length / ell;
        let coherence_ratio = r_vir * v_mag / nu;

        eprintln!(
            "halo {k}:  |v| = {v_mag:.1} kpc/Gyr,  m|v|Δx / ℓ = {nyquist_ratio:.3} (< π = {pi:.3}),  r_h v / ν = {coherence_ratio:.2}",
            pi = PI,
        );
    }

    // Pre-compute NFW normalization constants.
    let rho_0: Vec<f64> = (0..n_halos)
        .map(|k| {
            let c = halo_configs[k].concentration;
            let nfw_integral = (1.0 + c).ln() - c / (1.0 + c);
            halo_masses[k] / (4.0 * PI * scale_radii[k].powi(3) * nfw_integral)
        })
        .collect();

    // Build α as a coherent sum of per-halo plane-wave contributions.
    //
    // α_k(x) = sqrt(ρ_k(x) / m) * exp(I * m * v_k · (x - x_k) / l)
    //
    // The phase is referenced from each halo center so absolute phase
    // values stay small across each lump, keeping the discrete
    // representation clean.
    let morphis_grid = morphis::grid::Grid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();

    EvenField::from_fn(&morphis_grid, g, |x| {
        let mut scalar_sum = 0.0;
        let mut pseudo_sum = 0.0;

        for k in 0..n_halos {
            let dx = wrap_distance(x[0] - centers[k][0], box_length);
            let dy = wrap_distance(x[1] - centers[k][1], box_length);
            let dz = wrap_distance(x[2] - centers[k][2], box_length);
            let r = (dx * dx + dy * dy + dz * dz).sqrt().max(0.5 * cell_length);

            // NFW density profile with taper beyond virial radius.
            let s = r / scale_radii[k];
            let nfw_profile = 1.0 / (s * (1.0 + s) * (1.0 + s));
            let taper = if r < virial_radii[k] {
                1.0
            } else {
                let excess = (r - virial_radii[k]) / scale_radii[k];
                (-excess * excess).exp()
            };
            let rho_k = (rho_0[k] * nfw_profile * taper).max(density_mean * 0.01 / n_halos as f64);

            let amplitude = (rho_k / mass).sqrt();

            // Linear phase ramp: v_k · (x - x_k), using wrapped displacement.
            let phase =
                (velocities[k][0] * dx + velocities[k][1] * dy + velocities[k][2] * dz) / nu;

            scalar_sum += amplitude * phase.cos();
            pseudo_sum += amplitude * phase.sin();
        }

        (scalar_sum, pseudo_sum)
    })
}

/// Wrap a distance into [-L/2, L/2] for periodic boundary.
fn wrap_distance(d: f64, box_length: f64) -> f64 {
    let mut wrapped = d % box_length;
    if wrapped > box_length / 2.0 {
        wrapped -= box_length;
    }
    if wrapped < -box_length / 2.0 {
        wrapped += box_length;
    }

    wrapped
}
