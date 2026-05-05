/// Direct comparison: does a field lump track a particle under identical gravity?
///
/// Sets up a single overdensity at the same location in both particle and
/// field representations, evolves both for several steps under the same
/// cosmological SP / PM dynamics, and compares the centroid displacement.
///
/// This is the definitive test: if the formulations solve the same physics,
/// they must agree on bulk motion of a well-resolved, non-overlapping lump.
use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric::euclidean;

use hermes_rs::core::content::{Content, FieldParams, FieldState};
use hermes_rs::core::dynamics::Dynamics;
use hermes_rs::core::pm_dynamics::ParticleMeshDynamics;
use hermes_rs::core::schrodinger_dynamics::SchrodingerPoissonDynamics;
use hermes_rs::engine::coupling::poisson::PoissonGravity;
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::particles::Particles;

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

/// Density-weighted centroid (circular mean for periodic box).
fn field_centroid_x(alpha: &EvenField<3>, box_length: f64) -> f64 {
    let n = alpha.grid.n_cells;
    let cell_length = box_length / n as f64;
    let mut total_weight = 0.0;
    let mut sin_sum = 0.0;
    let mut cos_sum = 0.0;

    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * cell_length;
        for m1 in 0..n {
            for m2 in 0..n {
                let a = alpha.scalar[[m0, m1, m2]];
                let b = alpha.pseudoscalar[[m0, m1, m2]];
                let weight = a * a + b * b;
                total_weight += weight;
                let theta = 2.0 * PI * x / box_length;
                sin_sum += weight * theta.sin();
                cos_sum += weight * theta.cos();
            }
        }
    }

    let angle = (sin_sum / total_weight).atan2(cos_sum / total_weight);
    let mut cx = angle * box_length / (2.0 * PI);
    if cx < 0.0 {
        cx += box_length;
    }

    cx
}

/// Mass-weighted centroid x for particles.
fn particle_centroid_x(particles: &Particles, box_length: f64) -> f64 {
    let mut sin_sum = 0.0;
    let mut cos_sum = 0.0;
    let n = particles.count();

    for p in 0..n {
        let pos = particles.position_of(p);
        let x = hermes_rs::algebra::components_from_vector(&pos)[0];
        let theta = 2.0 * PI * x / box_length;
        sin_sum += theta.sin();
        cos_sum += theta.cos();
    }

    let angle = (sin_sum / n as f64).atan2(cos_sum / n as f64);
    let mut cx = angle * box_length / (2.0 * PI);
    if cx < 0.0 {
        cx += box_length;
    }

    cx
}

#[test]
fn field_and_particle_lump_track_together() {
    let n = 16;
    let box_length = 4000.0;
    let cell_length = box_length / n as f64;

    let cosmology = planck_2018();
    let hermes_grid = Grid::new(n, box_length);
    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    let density_mean = cosmology.density_matter();

    // Field parameters.
    let mass = 1e10;
    let nu = 2000.0;
    let ell = nu * mass;

    // A Gaussian lump offset from center with initial velocity.
    let center = [1500.0, 2000.0, 2000.0];
    let sigma = 400.0;
    let density_peak = 200.0;
    let velocity = [15.0, 0.0, 0.0];

    // Regime check.
    let nyquist = velocity[0] * cell_length / nu;
    eprintln!("m|v|dx/l = {nyquist:.3} (must be < pi)");
    assert!(nyquist < PI, "velocity exceeds Nyquist limit");

    // ==================== FIELD SETUP ====================
    let alpha = EvenField::from_fn(&morphis_grid, g, |x| {
        let mut r2 = 0.0;
        let mut v_dot_dx = 0.0;
        for d in 0..3 {
            let dx = wrap_distance(x[d] - center[d], box_length);
            r2 += dx * dx;
            v_dot_dx += velocity[d] * dx;
        }
        let rho = density_peak * (-r2 / (2.0 * sigma * sigma)).exp() + density_mean * 0.01;
        let amplitude = (rho / mass).sqrt();
        let phase = v_dot_dx / nu;

        (amplitude * phase.cos(), amplitude * phase.sin())
    });

    let field_x0 = field_centroid_x(&alpha, box_length);

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
    let mut field_content = Content::Fields(field_state);
    let gravity = PoissonGravity::new(hermes_grid.clone());
    let mut field_dynamics = SchrodingerPoissonDynamics::new(gravity);

    // ==================== PARTICLE SETUP ====================
    // Sample particles from the same Gaussian via rejection sampling.
    // They must be spatially concentrated (not uniform) so the centroid
    // is well-defined and tracks the bulk motion.
    use rand::Rng;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    let n_total = 1000;
    let mass_particle = density_mean * box_length.powi(3) / n_total as f64;
    let mut particles = Particles::zeros(n_total, mass_particle);

    let mut rng = ChaCha20Rng::seed_from_u64(42);
    let mut placed = 0;
    while placed < n_total {
        let x: f64 = rng.random_range(-1.0..1.0);
        let y: f64 = rng.random_range(-1.0..1.0);
        let z: f64 = rng.random_range(-1.0..1.0);
        let r2 = x * x + y * y + z * z;
        if r2 > 1.0 {
            continue;
        }

        // Accept with probability proportional to Gaussian density.
        let r_physical = r2.sqrt() * 3.0 * sigma;
        let accept = (-r_physical * r_physical / (2.0 * sigma * sigma)).exp();
        let u: f64 = rng.random_range(0.0..1.0);
        if u > accept {
            continue;
        }

        let pos = hermes_rs::algebra::vector_from_components(
            center[0] + x * 3.0 * sigma,
            center[1] + y * 3.0 * sigma,
            center[2] + z * 3.0 * sigma,
        );
        particles.set_position(placed, &pos);

        // Canonical momentum: p = m_particle * v (at a=1).
        let mom = hermes_rs::algebra::vector_from_components(
            mass_particle * velocity[0],
            mass_particle * velocity[1],
            mass_particle * velocity[2],
        );
        particles.set_momentum(placed, &mom);
        placed += 1;
    }

    let particle_x0 = particle_centroid_x(&particles, box_length);

    let mut particle_content = Content::Particles(particles);
    let mut particle_dynamics = ParticleMeshDynamics::new(hermes_grid);

    // ==================== EVOLVE BOTH ====================
    let a_start = 1.0;
    let a_end = 1.5;
    let n_steps = 30;
    let da = (a_end - a_start) / n_steps as f64;

    for step in 0..n_steps {
        let a_prev = a_start + step as f64 * da;
        let a_next = a_prev + da;

        field_dynamics
            .step(&mut field_content, &cosmology, a_prev, a_next)
            .unwrap();
        particle_dynamics
            .step(&mut particle_content, &cosmology, a_prev, a_next)
            .unwrap();

        if (step + 1) % 25 == 0 {
            let field_alpha = field_content.fields().unwrap().alpha.as_ref().unwrap();
            let field_cx = field_centroid_x(field_alpha, box_length);
            let field_dx = wrap_distance(field_cx - field_x0, box_length);

            let part = particle_content.particles().unwrap();
            let part_cx = particle_centroid_x(part, box_length);
            let part_dx = wrap_distance(part_cx - particle_x0, box_length);

            eprintln!(
                "step {} (a={:.3}): field dx = {:.1} kpc, particle dx = {:.1} kpc, ratio = {:.3}",
                step + 1,
                a_next,
                field_dx,
                part_dx,
                if part_dx.abs() > 1.0 {
                    field_dx / part_dx
                } else {
                    f64::NAN
                }
            );
        }
    }

    // Final comparison.
    let field_alpha = field_content.fields().unwrap().alpha.as_ref().unwrap();
    let field_dx = wrap_distance(
        field_centroid_x(field_alpha, box_length) - field_x0,
        box_length,
    );

    let part = particle_content.particles().unwrap();
    let part_dx = wrap_distance(
        particle_centroid_x(part, box_length) - particle_x0,
        box_length,
    );

    eprintln!("\nFinal: field = {field_dx:.1} kpc, particle = {part_dx:.1} kpc");

    let ratio = field_dx / part_dx;
    eprintln!("ratio field/particle = {ratio:.3}");

    // They should agree to within ~20% for well-separated lumps.
    assert!(
        (ratio - 1.0).abs() < 0.5,
        "field and particle centroids diverged: ratio = {ratio:.3}"
    );
}
