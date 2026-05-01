//! 2D diagnostic plots via plotters.

use std::path::Path;

use crate::colormap::colormap_hot;
use crate::physics::cic::assign_density;
use crate::physics::diagnostics::Diagnostics;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// Render a 2D projected density slice as a PNG image.
///
/// Projects a thin slab of thickness `slab_thickness` (in kpc) centered
/// at `slab_center` along the z-axis. The projection sums particle masses
/// in the slab and deposits onto a 2D histogram.
pub fn render_density_slice(
    particles: &Particles,
    grid: &Grid,
    slab_center: f64,
    slab_thickness: f64,
    output_path: &Path,
    resolution: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let half_slab = slab_thickness / 2.0;
    let bin_size = grid.box_length / resolution as f64;

    let mut histogram = vec![vec![0.0_f64; resolution]; resolution];

    for p in 0..particles.count() {
        let pos = particles.position_components(p);
        let z = pos[2];

        let dz = (z - slab_center)
            .abs()
            .min(grid.box_length - (z - slab_center).abs());
        if dz > half_slab {
            continue;
        }

        let bin_x = ((pos[0] / bin_size) as usize).min(resolution - 1);
        let bin_y = ((pos[1] / bin_size) as usize).min(resolution - 1);
        histogram[bin_x][bin_y] += particles.mass_particle;
    }

    let max_density = histogram
        .iter()
        .flat_map(|row| row.iter())
        .copied()
        .fold(0.0_f64, f64::max);
    let min_nonzero = histogram
        .iter()
        .flat_map(|row| row.iter())
        .copied()
        .filter(|&v| v > 0.0)
        .fold(f64::MAX, f64::min);

    let log_min = min_nonzero.max(1e-30).ln();
    let log_max = max_density.max(1e-30).ln();
    let log_range = (log_max - log_min).max(1e-10);

    let root =
        BitMapBackend::new(output_path, (resolution as u32, resolution as u32)).into_drawing_area();
    root.fill(&BLACK)?;

    for (bx, row) in histogram.iter().enumerate() {
        for (by, &mass) in row.iter().enumerate() {
            let normalized = if mass > 0.0 {
                ((mass.ln() - log_min) / log_range).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let color = colormap_hot(normalized);
            let r = (color[0] * 255.0) as u8;
            let g = (color[1] * 255.0) as u8;
            let b = (color[2] * 255.0) as u8;

            root.draw_pixel((bx as i32, by as i32), &RGBColor(r, g, b))?;
        }
    }

    root.present()?;

    Ok(())
}

/// Compute and plot the matter power spectrum P(k).
pub fn plot_power_spectrum(
    particles: &Particles,
    grid: &Grid,
    output_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::f64::consts::PI;

    use ndrustfft::{FftHandler, R2cFftHandler, ndfft, ndfft_r2c};
    use num_complex::Complex64;
    use plotters::prelude::*;

    let density = assign_density(particles, grid);
    let density_mean = density.sum() / grid.total_cells() as f64;
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;

    let mut overdensity = density;
    overdensity.data /= density_mean;
    overdensity.data -= 1.0;

    let mut overdensity_hat = ndarray::Array3::<Complex64>::zeros((n, n, n_complex));
    let handler_r2c = R2cFftHandler::new(n);
    let handler_c2c_1 = FftHandler::new(n);
    let handler_c2c_0 = FftHandler::new(n);

    ndfft_r2c(&overdensity.data, &mut overdensity_hat, &handler_r2c, 2);
    let mut scratch = overdensity_hat.clone();
    ndfft(&overdensity_hat, &mut scratch, &handler_c2c_1, 1);
    overdensity_hat.assign(&scratch);
    ndfft(&overdensity_hat, &mut scratch, &handler_c2c_0, 0);
    overdensity_hat.assign(&scratch);

    let k_nyquist = PI * n as f64 / grid.box_length;
    let n_bins = n / 2;
    let dk = k_nyquist / n_bins as f64;

    let mut power_sum = vec![0.0_f64; n_bins];
    let mut mode_count = vec![0_usize; n_bins];

    for m0 in 0..n {
        let kx = grid.wavevector_component(m0);
        for m1 in 0..n {
            let ky = grid.wavevector_component(m1);
            for m2 in 0..n_complex {
                let kz = grid.wavevector_component(m2);
                let k_mag = (kx * kx + ky * ky + kz * kz).sqrt();

                if k_mag < 1e-30 {
                    continue;
                }

                let bin = (k_mag / dk) as usize;
                if bin < n_bins {
                    let norm = overdensity_hat[[m0, m1, m2]].norm_sqr();
                    power_sum[bin] += norm;
                    mode_count[bin] += 1;
                }
            }
        }
    }

    let volume_box = grid.box_volume();
    let k_values: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * dk).collect();
    let power_values: Vec<f64> = (0..n_bins)
        .map(|b| {
            if mode_count[b] > 0 {
                volume_box * power_sum[b] / mode_count[b] as f64
            } else {
                0.0
            }
        })
        .collect();

    let points: Vec<(f64, f64)> = k_values
        .iter()
        .zip(power_values.iter())
        .filter(|(_, power)| **power > 0.0)
        .map(|(k, power)| (*k, *power))
        .collect();

    if points.is_empty() {
        return Ok(());
    }

    let k_min = points.first().unwrap().0;
    let k_max = points.last().unwrap().0;
    let power_min = points.iter().map(|(_, p)| *p).fold(f64::MAX, f64::min);
    let power_max = points.iter().map(|(_, p)| *p).fold(0.0_f64, f64::max);

    let root = BitMapBackend::new(output_path, (800, 600)).into_drawing_area();
    root.fill(&BLACK)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Matter Power Spectrum",
            ("sans-serif", 24).into_font().color(&WHITE),
        )
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (k_min..k_max).log_scale(),
            (power_min * 0.5..power_max * 2.0).log_scale(),
        )?;

    chart
        .configure_mesh()
        .x_desc("k (1/kpc)")
        .y_desc("P(k) (kpc³)")
        .label_style(("sans-serif", 14).into_font().color(&WHITE))
        .axis_desc_style(("sans-serif", 16).into_font().color(&WHITE))
        .light_line_style(RGBColor(40, 40, 40))
        .bold_line_style(RGBColor(80, 80, 80))
        .draw()?;

    chart.draw_series(LineSeries::new(points, &CYAN))?;

    root.present()?;

    Ok(())
}

/// Plot conservation diagnostics over the simulation history.
pub fn plot_conservation(
    diagnostics: &[Diagnostics],
    output_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    if diagnostics.is_empty() {
        return Ok(());
    }

    let scale_factors: Vec<f64> = diagnostics.iter().map(|d| d.scale_factor).collect();
    let momentum_magnitudes: Vec<f64> =
        diagnostics.iter().map(|d| d.momentum_magnitude()).collect();
    let kinetic_energies: Vec<f64> = diagnostics.iter().map(|d| d.energy_kinetic).collect();
    let potential_energies: Vec<f64> = diagnostics.iter().map(|d| d.energy_potential).collect();
    let total_energies: Vec<f64> = diagnostics.iter().map(|d| d.energy_total()).collect();

    let a_min = scale_factors.first().copied().unwrap_or(0.0);
    let a_max = scale_factors.last().copied().unwrap_or(1.0);

    let root = BitMapBackend::new(output_path, (800, 900)).into_drawing_area();
    root.fill(&BLACK)?;

    let panels = root.split_evenly((3, 1));

    // Panel 1: Momentum magnitude
    {
        let max_momentum = momentum_magnitudes
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
            .max(1e-30);
        let mut chart = ChartBuilder::on(&panels[0])
            .caption(
                "Total Momentum",
                ("sans-serif", 18).into_font().color(&WHITE),
            )
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(60)
            .build_cartesian_2d(a_min..a_max, 0.0..max_momentum * 1.2)?;

        chart
            .configure_mesh()
            .label_style(("sans-serif", 12).into_font().color(&WHITE))
            .light_line_style(RGBColor(40, 40, 40))
            .draw()?;

        let points: Vec<(f64, f64)> = scale_factors
            .iter()
            .zip(momentum_magnitudes.iter())
            .map(|(&a, &p)| (a, p))
            .collect();
        chart.draw_series(LineSeries::new(points, &CYAN))?;
    }

    // Panel 2: Energies
    {
        let all_energies: Vec<f64> = kinetic_energies
            .iter()
            .chain(potential_energies.iter())
            .chain(total_energies.iter())
            .copied()
            .collect();
        let energy_min = all_energies.iter().copied().fold(f64::MAX, f64::min);
        let energy_max = all_energies.iter().copied().fold(f64::MIN, f64::max);
        let margin = (energy_max - energy_min).abs() * 0.1;

        let mut chart = ChartBuilder::on(&panels[1])
            .caption("Energy", ("sans-serif", 18).into_font().color(&WHITE))
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(60)
            .build_cartesian_2d(a_min..a_max, (energy_min - margin)..(energy_max + margin))?;

        chart
            .configure_mesh()
            .label_style(("sans-serif", 12).into_font().color(&WHITE))
            .light_line_style(RGBColor(40, 40, 40))
            .draw()?;

        let kinetic_points: Vec<(f64, f64)> = scale_factors
            .iter()
            .zip(kinetic_energies.iter())
            .map(|(&a, &e)| (a, e))
            .collect();
        let potential_points: Vec<(f64, f64)> = scale_factors
            .iter()
            .zip(potential_energies.iter())
            .map(|(&a, &e)| (a, e))
            .collect();
        let total_points: Vec<(f64, f64)> = scale_factors
            .iter()
            .zip(total_energies.iter())
            .map(|(&a, &e)| (a, e))
            .collect();

        chart
            .draw_series(LineSeries::new(kinetic_points, &RED))?
            .label("Kinetic")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
        chart
            .draw_series(LineSeries::new(potential_points, &BLUE))?
            .label("Potential")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));
        chart
            .draw_series(LineSeries::new(total_points, &WHITE))?
            .label("Total")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], WHITE));

        chart
            .configure_series_labels()
            .background_style(RGBColor(20, 20, 20).mix(0.8))
            .label_font(("sans-serif", 12).into_font().color(&WHITE))
            .border_style(WHITE)
            .draw()?;
    }

    // Panel 3: Angular momentum magnitude
    {
        let angular_momentum_magnitudes: Vec<f64> = diagnostics
            .iter()
            .map(|d| d.angular_momentum_magnitude())
            .collect();
        let max_angular_momentum = angular_momentum_magnitudes
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
            .max(1e-30);

        let mut chart = ChartBuilder::on(&panels[2])
            .caption(
                "Angular Momentum |L|",
                ("sans-serif", 18).into_font().color(&WHITE),
            )
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(60)
            .build_cartesian_2d(a_min..a_max, 0.0..max_angular_momentum * 1.2)?;

        chart
            .configure_mesh()
            .x_desc("Scale factor a")
            .label_style(("sans-serif", 12).into_font().color(&WHITE))
            .axis_desc_style(("sans-serif", 14).into_font().color(&WHITE))
            .light_line_style(RGBColor(40, 40, 40))
            .draw()?;

        let points: Vec<(f64, f64)> = scale_factors
            .iter()
            .zip(angular_momentum_magnitudes.iter())
            .map(|(&a, &l)| (a, l))
            .collect();
        chart.draw_series(LineSeries::new(points, &YELLOW))?;
    }

    root.present()?;

    Ok(())
}
