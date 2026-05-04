//! Galaxy group field initial conditions.
//!
//! Places NFW density profiles on the grid with bulk infall velocities,
//! then converts to a wavefunction via the inverse Madelung transform.
//! The density and halo geometry match the particle-mesh galaxy group
//! scene so the two representations start from the same physical state.

use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::field::Field;
use morphis::metric;
use ndarray::Array3;
use num_complex::Complex64;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::core::content::FieldParams;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid as HermesGrid;
use crate::physics::spectral::ifft_3d;
use crate::scenes::galaxy_group_pm::init::default_halo_configs;

/// Generate colliding NFW halos as a wavefunction field.
///
/// The density is a superposition of NFW profiles placed at the same
/// positions as the particle scene. The velocity is encoded as a
/// velocity potential (smooth, periodic) via the inverse Madelung map.
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

    // Halo mass fractions → absolute masses.
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

    // Build density on the grid: superposition of NFW profiles.
    let mut density = Array3::<f64>::zeros((n, n, n));

    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * cell_length;
        for m1 in 0..n {
            let y = (m1 as f64 + 0.5) * cell_length;
            for m2 in 0..n {
                let z = (m2 as f64 + 0.5) * cell_length;

                let mut rho = 0.0;
                for k in 0..n_halos {
                    let dx = wrap_distance(x - centers[k][0], box_length);
                    let dy = wrap_distance(y - centers[k][1], box_length);
                    let dz = wrap_distance(z - centers[k][2], box_length);
                    let r = (dx * dx + dy * dy + dz * dz).sqrt().max(0.5 * cell_length);

                    let s = r / scale_radii[k];
                    let nfw_profile = 1.0 / (s * (1.0 + s) * (1.0 + s));

                    // Normalize so the halo integrates to its mass within r_vir.
                    let c = halo_configs[k].concentration;
                    let nfw_integral = (1.0 + c).ln() - c / (1.0 + c);
                    let rho_0 = halo_masses[k] / (4.0 * PI * scale_radii[k].powi(3) * nfw_integral);

                    // Taper beyond virial radius.
                    let taper = if r < virial_radii[k] {
                        1.0
                    } else {
                        let excess = (r - virial_radii[k]) / scale_radii[k];
                        (-excess * excess).exp()
                    };

                    rho += rho_0 * nfw_profile * taper;
                }

                // Floor at the mean density to avoid near-zero regions.
                density[[m0, m1, m2]] = rho.max(density_mean * 0.01);
            }
        }
    }

    // Normalize total mass to match the box.
    let cell_volume = cell_length * cell_length * cell_length;
    let mass_on_grid: f64 = density.iter().sum::<f64>() * cell_volume;
    density *= mass_total / mass_on_grid;

    // Build the velocity potential on the grid.
    //
    // Each halo has an infall velocity toward the center of mass.
    // The velocity field is v = grad(phi_v). We build phi_v by
    // solving for it spectrally from the velocity divergence:
    //   div(v) = laplacian(phi_v)  →  phi_v = laplacian_inverse(div(v))
    let mut div_v = Array3::<f64>::zeros((n, n, n));

    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * cell_length;
        for m1 in 0..n {
            let y = (m1 as f64 + 0.5) * cell_length;
            for m2 in 0..n {
                let z = (m2 as f64 + 0.5) * cell_length;

                // Compute the velocity divergence as a weighted sum
                // of halo contributions, each flowing toward the COM.
                let mut local_div = 0.0;
                for k in 0..n_halos {
                    let dx = wrap_distance(x - centers[k][0], box_length);
                    let dy = wrap_distance(y - centers[k][1], box_length);
                    let dz = wrap_distance(z - centers[k][2], box_length);
                    let r = (dx * dx + dy * dy + dz * dz).sqrt().max(0.5 * cell_length);

                    // Weight by the NFW density profile (velocity follows mass).
                    let s = r / scale_radii[k];
                    let weight = 1.0 / (s * (1.0 + s) * (1.0 + s));

                    // Infall speed toward COM.
                    let d_com_x = com[0] - centers[k][0];
                    let d_com_y = com[1] - centers[k][1];
                    let d_com_z = com[2] - centers[k][2];
                    let d_com = (d_com_x * d_com_x + d_com_y * d_com_y + d_com_z * d_com_z)
                        .sqrt()
                        .max(1.0);
                    let infall_speed = (g_const * total_mass / d_com).sqrt() * 0.4;

                    // The divergence of a radial infall is negative (converging).
                    // div(v_infall) ~ -3 * v_infall / r for uniform infall.
                    let taper = if r < virial_radii[k] {
                        1.0
                    } else {
                        let excess = (r - virial_radii[k]) / scale_radii[k];
                        (-excess * excess).exp()
                    };

                    local_div += -3.0 * infall_speed / d_com * weight * taper;
                }

                div_v[[m0, m1, m2]] = local_div;
            }
        }
    }

    // Solve for phi_v spectrally: phi_v = laplacian_inverse(div_v).
    let n_complex = n / 2 + 1;
    let div_v_hat = crate::physics::spectral::fft_3d(&div_v, n);

    let mut phi_v_hat = Array3::<Complex64>::zeros((n, n, n_complex));
    for m0 in 0..n {
        let kx = grid.wavevector_component(m0);
        for m1 in 0..n {
            let ky = grid.wavevector_component(m1);
            for m2 in 0..n_complex {
                let kz = grid.wavevector_component(m2);
                let k2 = kx * kx + ky * ky + kz * kz;
                if k2 < 1e-30 {
                    continue;
                }
                phi_v_hat[[m0, m1, m2]] = div_v_hat[[m0, m1, m2]] / (-k2);
            }
        }
    }
    let phi_v = ifft_3d(&phi_v_hat, n);

    // Inverse Madelung transform: α = sqrt(ρ/m) * exp(I * m * phi_v / l).
    let morphis_grid = morphis::grid::Grid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let nu = ell / mass;

    let rho_field = Field::scalar_field(&morphis_grid, g, |x| {
        let m0 = ((x[0] / cell_length) as usize).min(n - 1);
        let m1 = ((x[1] / cell_length) as usize).min(n - 1);
        let m2 = ((x[2] / cell_length) as usize).min(n - 1);
        density[[m0, m1, m2]].max(1e-30)
    });

    let phi_v_field = Field::scalar_field(&morphis_grid, g, |x| {
        let m0 = ((x[0] / cell_length) as usize).min(n - 1);
        let m1 = ((x[1] / cell_length) as usize).min(n - 1);
        let m2 = ((x[2] / cell_length) as usize).min(n - 1);
        phi_v[[m0, m1, m2]]
    });

    EvenField::madelung_inverse(&rho_field, &phi_v_field, mass, nu)
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
