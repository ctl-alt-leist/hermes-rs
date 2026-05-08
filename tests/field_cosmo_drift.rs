/// Diagnostic test: does a field lump translate under cosmological SP dynamics?
///
/// This test isolates whether the full SchrodingerPoissonDynamics
/// (kinetic + potential with cosmological scale factors) produces
/// bulk translation of a density lump with initial velocity.
use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric::euclidean;

use hermes_rs::core::content::{Content, FieldParams, FieldState};
use hermes_rs::core::dynamics::Dynamics;
use hermes_rs::core::schrodinger_dynamics::SchrodingerPoissonDynamics;
use hermes_rs::engine::coupling::poisson::PoissonGravity;
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;

/// Density-weighted centroid (circular mean for periodic box).
fn field_centroid(alpha: &EvenField<3>, box_length: f64) -> [f64; 3] {
    let n = alpha.grid.n_cells;
    let cell_length = box_length / n as f64;
    let mut total_weight = 0.0;
    let mut sin_sum = [0.0; 3];
    let mut cos_sum = [0.0; 3];

    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * cell_length;
        for m1 in 0..n {
            let y = (m1 as f64 + 0.5) * cell_length;
            for m2 in 0..n {
                let z = (m2 as f64 + 0.5) * cell_length;
                let a = alpha.scalar[[m0, m1, m2]];
                let b = alpha.pseudoscalar[[m0, m1, m2]];
                let weight = a * a + b * b;
                total_weight += weight;
                for (d, pos) in [x, y, z].iter().enumerate() {
                    let theta = 2.0 * PI * pos / box_length;
                    sin_sum[d] += weight * theta.sin();
                    cos_sum[d] += weight * theta.cos();
                }
            }
        }
    }

    let mut centroid = [0.0; 3];
    for d in 0..3 {
        let angle = (sin_sum[d] / total_weight).atan2(cos_sum[d] / total_weight);
        centroid[d] = angle * box_length / (2.0 * PI);
        if centroid[d] < 0.0 {
            centroid[d] += box_length;
        }
    }

    centroid
}

fn wrap_distance(d: f64, l: f64) -> f64 {
    let mut w = d % l;
    if w > l / 2.0 {
        w -= l;
    }
    if w < -l / 2.0 {
        w += l;
    }

    w
}

#[test]
fn lump_translates_under_cosmological_sp_dynamics() {
    let n = 16;
    let box_length = 4000.0;
    let cell_length = box_length / n as f64;

    let cosmology = planck_2018();
    let hermes_grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let mass = 1e10;
    let nu = 2000.0;
    let ell = nu * mass;

    let sigma = 400.0;
    let density_peak = 50.0;
    let density_floor = cosmology.density_matter() * 0.01;
    let velocity = [20.0, 0.0, 0.0];
    let center = [2000.0, 2000.0, 2000.0];

    // Regime diagnostics.
    let v_mag = velocity[0];
    eprintln!("m|v|Δx / ℓ = {:.3} (limit: π)", v_mag * cell_length / nu);
    eprintln!("r v / ν = {:.2}", sigma * v_mag / nu);

    let alpha = EvenField::from_fn(&morphis_grid, g, |x| {
        let mut r2 = 0.0;
        let mut v_dot_dx = 0.0;
        for d in 0..3 {
            let dx = wrap_distance(x[d] - center[d], box_length);
            r2 += dx * dx;
            v_dot_dx += velocity[d] * dx;
        }
        let rho = density_peak * (-r2 / (2.0 * sigma * sigma)).exp() + density_floor;
        let amplitude = (rho / mass).sqrt();
        let phase = v_dot_dx / nu;

        (amplitude * phase.cos(), amplitude * phase.sin())
    });

    let centroid_0 = field_centroid(&alpha, box_length);
    eprintln!("t=0 (a=1.000): centroid_x = {:.1}", centroid_0[0]);

    let field_state = FieldState {
        grid: morphis_grid,
        alpha: Some(alpha),
        beta: None,
        gamma: None,
        params: FieldParams {
            smoothing_length: ell,
            mass_alpha: mass,
        },
    };

    let mut content = Content::Fields(field_state);
    let gravity = PoissonGravity::new(hermes_grid);
    let mut dynamics = SchrodingerPoissonDynamics::new(gravity);

    // Evolve from a=1.0 to a=1.5 in 200 steps (moderate expansion).
    let a_start = 1.0;
    let a_end = 1.5;
    let n_steps = 50;
    let da = (a_end - a_start) / n_steps as f64;

    for step in 0..n_steps {
        let a_prev = a_start + step as f64 * da;
        let a_next = a_prev + da;

        dynamics
            .step(&mut content, &cosmology, a_prev, a_next)
            .unwrap();

        if (step + 1) % 50 == 0 {
            let alpha = content.fields().unwrap().alpha.as_ref().unwrap();
            let centroid = field_centroid(alpha, box_length);
            let displacement = wrap_distance(centroid[0] - centroid_0[0], box_length);
            eprintln!(
                "step {} (a={:.3}): centroid_x = {:.1}, displacement = {:.1} kpc",
                step + 1,
                a_next,
                centroid[0],
                displacement
            );
        }
    }

    let alpha = content.fields().unwrap().alpha.as_ref().unwrap();
    let centroid_final = field_centroid(alpha, box_length);
    let total_displacement = wrap_distance(centroid_final[0] - centroid_0[0], box_length);

    eprintln!("total displacement: {total_displacement:.1} kpc");

    // The lump should have moved substantially.
    assert!(
        total_displacement.abs() > 50.0,
        "lump did not translate under cosmological SP: displacement = {total_displacement:.1} kpc"
    );
}
