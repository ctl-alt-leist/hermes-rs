//! Density colormapping for cosmological visualization.

use crate::physics::cic::assign_density;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// Map a normalized value in [0, 1] to an RGB color on a hot colormap.
///
/// The colormap runs: black → deep blue → cyan → white, designed
/// for cosmological density fields on a dark background.
pub fn colormap_hot(value: f64) -> [f32; 3] {
    let t = value.clamp(0.0, 1.0) as f32;

    let r = (3.0 * t - 1.0).clamp(0.0, 1.0);
    let g = (3.0 * t - 2.0).clamp(0.0, 1.0);
    let b = (2.0 * t)
        .clamp(0.0, 1.0)
        .min(1.0 - (3.0 * t - 2.5).clamp(0.0, 1.0));

    [r, g, b]
}

/// Compute per-particle density estimates by CIC-depositing onto the grid
/// and interpolating back. Returns a normalized [0, 1] density for each
/// particle (log-scaled relative to the mean).
pub fn particle_density_colors(particles: &Particles, grid: &Grid) -> Vec<[f32; 3]> {
    let density = assign_density(particles, grid);
    let density_mean = density.sum() / grid.total_cells() as f64;

    let h_inv = 1.0 / grid.cell_length;
    let n = grid.n_cells;

    let mut density_per_particle = vec![0.0_f64; particles.count()];

    for (p, density_p) in density_per_particle.iter_mut().enumerate() {
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

        let mut rho = 0.0;
        for (a, &weight_a) in weight[0].iter().enumerate() {
            let g0 = ((base[0] + a as isize) % n as isize + n as isize) as usize % n;
            for (b, &weight_b) in weight[1].iter().enumerate() {
                let g1 = ((base[1] + b as isize) % n as isize + n as isize) as usize % n;
                for (c, &weight_c) in weight[2].iter().enumerate() {
                    let g2 = ((base[2] + c as isize) % n as isize + n as isize) as usize % n;
                    rho += density.data[[g0, g1, g2]] * weight_a * weight_b * weight_c;
                }
            }
        }

        *density_p = rho;
    }

    let log_min = (density_mean * 0.1).ln();
    let log_max = density_per_particle
        .iter()
        .copied()
        .fold(0.0_f64, f64::max)
        .max(density_mean)
        .ln();
    let log_range = (log_max - log_min).max(1e-10);

    density_per_particle
        .iter()
        .map(|&rho| {
            let log_rho = rho.max(1e-30).ln();
            let normalized = ((log_rho - log_min) / log_range).clamp(0.0, 1.0);

            colormap_hot(normalized)
        })
        .collect()
}
