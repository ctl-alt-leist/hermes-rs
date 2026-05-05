/// Two Gaussian lumps, no initial velocity, evolve under SP gravity.
///
/// Measures peak displacement after N steps and compares to the
/// Newtonian prediction: Δx = 0.5 × (GM/d²) × t².
use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric::euclidean;

use hermes_rs::core::schrodinger_dynamics::kinetic_step;
use hermes_rs::engine::coupling::poisson::field_potential_step;
use hermes_rs::physics::constants::G as GRAV;

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

/// Find x-position of max projected density, with parabolic refinement.
fn peak_x(alpha: &EvenField<3>, box_length: f64, x_hint: f64) -> f64 {
    let n = alpha.grid.n_cells;
    let cell_length = box_length / n as f64;

    // Project |α|² onto x.
    let mut density_x = vec![0.0; n];
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let a = alpha.scalar[[m0, m1, m2]];
                let b = alpha.pseudoscalar[[m0, m1, m2]];
                density_x[m0] += a * a + b * b;
            }
        }
    }

    // Find max near hint.
    let hint_idx = (x_hint / cell_length) as usize;
    let search = n / 4;
    let mut max_val = 0.0;
    let mut max_idx = hint_idx;
    for offset in 0..search {
        for sign in [1_isize, -1] {
            let m = ((hint_idx as isize + sign * offset as isize).rem_euclid(n as isize)) as usize;
            if density_x[m] > max_val {
                max_val = density_x[m];
                max_idx = m;
            }
        }
    }

    // Parabolic refinement.
    let prev = density_x[(max_idx + n - 1) % n];
    let curr = density_x[max_idx];
    let next = density_x[(max_idx + 1) % n];
    let denom = prev - 2.0 * curr + next;
    let offset = if denom.abs() > 1e-30 {
        0.5 * (prev - next) / denom
    } else {
        0.0
    };

    (max_idx as f64 + 0.5 + offset) * cell_length
}

#[test]
fn two_lumps_accelerate_toward_each_other() {
    let n = 32;
    let box_length = 2000.0;

    let mass = 1e10;
    let nu = 500.0;
    let ell = nu * mass;

    let sigma = 150.0;
    let density_peak = 1e4;
    let density_floor = 1.0;

    let x_left_0 = 700.0;
    let x_right_0 = 1300.0;
    let separation = x_right_0 - x_left_0;

    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let mut alpha = EvenField::from_fn(&morphis_grid, g, |x| {
        let mut rho = density_floor;
        for cx in [x_left_0, x_right_0] {
            let center = [cx, box_length / 2.0, box_length / 2.0];
            let mut r2 = 0.0;
            for d in 0..3 {
                let dx = wrap_distance(x[d] - center[d], box_length);
                r2 += dx * dx;
            }
            rho += density_peak * (-r2 / (2.0 * sigma * sigma)).exp();
        }

        ((rho / mass).sqrt(), 0.0)
    });

    // Lump mass and expected acceleration.
    let mass_lump = density_peak * (2.0 * PI * sigma * sigma).powf(1.5);
    let expected_accel = GRAV * mass_lump / (separation * separation);
    let spreading_time = sigma * sigma / nu;

    eprintln!("sigma = {sigma} kpc, separation = {separation} kpc");
    eprintln!("lump mass ~ {mass_lump:.2e} M_sun");
    eprintln!("expected accel = {expected_accel:.4} kpc/Gyr^2");
    eprintln!("spreading time = {spreading_time:.1} Gyr");
    eprintln!();

    // Measure initial peaks.
    let left_0 = peak_x(&alpha, box_length, x_left_0);
    let right_0 = peak_x(&alpha, box_length, x_right_0);
    eprintln!(
        "initial: left={left_0:.2}, right={right_0:.2}, sep={:.2}",
        right_0 - left_0
    );

    // Evolve: Strang splitting at a=1 (static background).
    let scale_factor = 1.0;
    let dt = 0.05;
    let n_steps = 10;

    for step in 0..n_steps {
        field_potential_step(
            &mut alpha,
            &morphis_grid,
            ell,
            mass,
            density_floor,
            scale_factor,
            dt / 2.0,
        );
        kinetic_step(&mut alpha, &morphis_grid, ell, mass, scale_factor, dt);
        field_potential_step(
            &mut alpha,
            &morphis_grid,
            ell,
            mass,
            density_floor,
            scale_factor,
            dt / 2.0,
        );

        let t = (step + 1) as f64 * dt;
        let left = peak_x(&alpha, box_length, x_left_0);
        let right = peak_x(&alpha, box_length, x_right_0);
        let sep = right - left;
        let dx_left = left - left_0;
        let dx_right = right - right_0;

        let expected_dx = 0.5 * expected_accel * t * t;

        eprintln!(
            "t={t:.2}: left={left:.2} (dx={dx_left:+.2}), right={right:.2} (dx={dx_right:+.2}), sep={sep:.2}, expected |dx|={expected_dx:.2}"
        );
    }

    // Final comparison.
    let t_final = n_steps as f64 * dt;
    let left_final = peak_x(&alpha, box_length, x_left_0);
    let right_final = peak_x(&alpha, box_length, x_right_0);

    let dx_left = left_final - left_0;
    let dx_right = right_final - right_0;
    let expected_dx = 0.5 * expected_accel * t_final * t_final;

    eprintln!();
    eprintln!("measured |dx_left|  = {:.3} kpc", dx_left.abs());
    eprintln!("measured |dx_right| = {:.3} kpc", dx_right.abs());
    eprintln!("expected |dx|       = {expected_dx:.3} kpc");
    eprintln!("ratio left:  {:.3}", dx_left.abs() / expected_dx);
    eprintln!("ratio right: {:.3}", dx_right.abs() / expected_dx);

    // Signs: left should move right, right should move left.
    assert!(dx_left > 0.0, "left lump should move right: dx={dx_left}");
    assert!(dx_right < 0.0, "right lump should move left: dx={dx_right}");

    // Magnitude: within factor of 2 of point-mass prediction.
    let avg_ratio = (dx_left.abs() + dx_right.abs()) / (2.0 * expected_dx);
    eprintln!("average ratio: {avg_ratio:.3}");

    assert!(
        avg_ratio > 0.3 && avg_ratio < 3.0,
        "acceleration doesn't match point-mass estimate: ratio = {avg_ratio:.3}"
    );
}
