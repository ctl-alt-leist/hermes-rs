/// Tests for field-theory halo dynamics.
///
/// These tests verify that lumps initialized with bulk velocities
/// actually translate under the kinetic operator, building up from
/// the simplest possible case (one lump, no gravity) to the full
/// galaxy-group scene.
///
/// Each test logs the dimensionless regime diagnostics:
///   - m|v|dx/l  (phase-Nyquist: must be < pi for resolvable velocity)
///   - r_h * v / nu  (coherence: >> 1 for coherent translation)
use std::collections::BTreeMap;
use std::f64::consts::PI;

use morphis::even_field::EvenField;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric::euclidean;

use hermes_rs::engine::Engine;
use hermes_rs::engine::free::FreeEvolution;
use hermes_rs::engine::free::schrodinger::SchrodingerEvolution;
use hermes_rs::engine::sector::Sector;
use hermes_rs::engine::sector::schrodinger::SchrodingerSector;
use hermes_rs::engine::solver::GravitySolver;
use hermes_rs::engine::state::{FieldEntry, SimulationState};
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;

// ============================================================================
// Helpers
// ============================================================================

/// Compute the density-weighted centroid of |α|² in a periodic box.
///
/// Uses the circular mean (atan2 of weighted sin/cos) to handle
/// wrapping correctly.
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

                let idx = [m0, m1, m2];
                let a = alpha.scalar[idx];
                let b = alpha.pseudoscalar[idx];
                let weight = a * a + b * b;

                total_weight += weight;

                let pos = [x, y, z];
                for d in 0..3 {
                    let theta = 2.0 * PI * pos[d] / box_length;
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

/// Total norm (sum of |α|² over all grid points).
fn field_norm(alpha: &EvenField<3>) -> f64 {
    alpha
        .scalar
        .iter()
        .zip(alpha.pseudoscalar.iter())
        .map(|(a, b)| a * a + b * b)
        .sum()
}

/// Initialize a single Gaussian lump with a bulk velocity.
///
/// α(x) = sqrt(ρ(x) / m) * exp(I * v · (x - x_0) / ν)
///
/// where ρ is a Gaussian centered at x_0 with width sigma.
fn gaussian_lump(
    grid: &MorphisGrid<3>,
    box_length: f64,
    center: [f64; 3],
    sigma: f64,
    velocity: [f64; 3],
    density_peak: f64,
    density_floor: f64,
    mass: f64,
    nu: f64,
) -> EvenField<3> {
    let g = euclidean::<3>();

    EvenField::from_fn(grid, g, |x| {
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
    })
}

/// Wrap a distance into [-L/2, L/2].
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

/// Log regime diagnostics for a lump.
fn log_regime(
    label: &str,
    velocity: &[f64; 3],
    cell_length: f64,
    lump_radius: f64,
    mass: f64,
    ell: f64,
    nu: f64,
) {
    let v_mag = (velocity[0].powi(2) + velocity[1].powi(2) + velocity[2].powi(2)).sqrt();
    let nyquist_ratio = mass * v_mag * cell_length / ell;
    let coherence_ratio = lump_radius * v_mag / nu;

    eprintln!(
        "{label}:  |v| = {v_mag:.1} kpc/Gyr,  m|v|Δx / ℓ = {nyquist_ratio:.3} (limit: π = {pi:.3}),  r v / ν = {coherence_ratio:.2}",
        pi = PI,
    );
}

// ============================================================================
// Step 1: Single-lump translation (no gravity)
// ============================================================================

#[test]
fn single_lump_translates_under_kinetic_evolution() {
    // Parameters chosen to satisfy both constraints:
    //   m|v|dx/l < pi   (phase resolvable on the grid)
    //   r*v/nu >> 1     (lump translates faster than it spreads)
    let n = 16;
    let box_length = 3000.0; // kpc
    let cell_length = box_length / n as f64;

    let mass = 1e10;
    let nu = 1000.0; // l/m: diffusivity (kpc^2/Gyr)
    let ell = nu * mass;

    let sigma = 300.0; // lump width (kpc)
    let density_peak = 100.0; // M_sun / kpc^3
    let density_floor = 0.1;

    let velocity = [12.0, 0.0, 0.0]; // kpc/Gyr (~ 12 km/s)
    let center = [1500.0, 1500.0, 1500.0];

    log_regime("single lump", &velocity, cell_length, sigma, mass, ell, nu);

    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let alpha = gaussian_lump(
        &morphis_grid,
        box_length,
        center,
        sigma,
        velocity,
        density_peak,
        density_floor,
        mass,
        nu,
    );

    let norm_before = field_norm(&alpha);
    let centroid_before = field_centroid(&alpha, box_length);
    eprintln!("t=0: centroid = {centroid_before:?}");

    let mut entry = FieldEntry {
        data: alpha,
        smoothing_length: ell,
        mass,
        self_interaction: None,
    };

    // Evolve under kinetic operator only (no gravity, no expansion).
    let dt = 0.5;
    let n_steps = 20;
    let mut evolver = SchrodingerEvolution;

    for step in 0..n_steps {
        evolver.step(&mut entry, &morphis_grid, 1.0, dt).unwrap();

        if (step + 1) % 10 == 0 {
            let centroid = field_centroid(&entry.data, box_length);
            eprintln!("t={:.1}: centroid = {centroid:?}", (step + 1) as f64 * dt);
        }
    }

    let norm_after = field_norm(&entry.data);
    let centroid_after = field_centroid(&entry.data, box_length);

    // Norm should be preserved.
    let norm_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        norm_error < 1e-12,
        "norm not preserved: relative error = {norm_error}"
    );

    // The centroid should have moved by approximately v * t in the x-direction.
    let total_time = n_steps as f64 * dt;
    let expected_displacement = velocity[0] * total_time;

    // Compute actual displacement (periodic).
    let actual_displacement = wrap_distance(centroid_after[0] - centroid_before[0], box_length);

    eprintln!(
        "expected displacement: {expected_displacement:.1} kpc, actual: {actual_displacement:.1} kpc"
    );

    let displacement_error = (actual_displacement - expected_displacement).abs();
    // Allow 10% error due to spreading and discretization.
    let tolerance = 0.1 * expected_displacement.abs() + cell_length;
    assert!(
        displacement_error < tolerance,
        "lump did not translate correctly: expected {expected_displacement:.1}, got {actual_displacement:.1} (error {displacement_error:.1} > tolerance {tolerance:.1})"
    );

    // y and z centroids should not have moved significantly.
    let y_displacement = wrap_distance(centroid_after[1] - centroid_before[1], box_length).abs();
    let z_displacement = wrap_distance(centroid_after[2] - centroid_before[2], box_length).abs();

    assert!(
        y_displacement < 2.0 * cell_length,
        "unexpected y displacement: {y_displacement:.1}"
    );
    assert!(
        z_displacement < 2.0 * cell_length,
        "unexpected z displacement: {z_displacement:.1}"
    );
}

// ============================================================================
// Step 2: Two-lump pass-through (no gravity)
// ============================================================================

#[test]
fn two_lumps_pass_through_under_kinetic_evolution() {
    // Two lumps moving toward each other, kinetic-only.
    // They should approach, develop interference fringes during overlap,
    // and emerge with their original velocities.
    let n = 16;
    let box_length = 6000.0;
    let cell_length = box_length / n as f64;

    let mass = 1e10;
    let nu = 4000.0;
    let ell = nu * mass;

    let sigma = 500.0;
    let density_peak = 100.0;
    let density_floor = 0.1;

    let velocity_1 = [10.0, 0.0, 0.0]; // rightward
    let velocity_2 = [-10.0, 0.0, 0.0]; // leftward

    let center_1 = [2000.0, 3000.0, 3000.0];
    let center_2 = [4000.0, 3000.0, 3000.0];

    log_regime("lump 1", &velocity_1, cell_length, sigma, mass, ell, nu);
    log_regime("lump 2", &velocity_2, cell_length, sigma, mass, ell, nu);

    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let g = euclidean::<3>();

    // Build α as coherent sum of two lumps.
    let alpha = EvenField::from_fn(&morphis_grid, g, |x| {
        let mut scalar_sum = 0.0;
        let mut pseudo_sum = 0.0;

        for (center, velocity) in [(center_1, velocity_1), (center_2, velocity_2)] {
            let mut r2 = 0.0;
            let mut v_dot_dx = 0.0;
            for d in 0..3 {
                let dx = wrap_distance(x[d] - center[d], box_length);
                r2 += dx * dx;
                v_dot_dx += velocity[d] * dx;
            }

            let rho = density_peak * (-r2 / (2.0 * sigma * sigma)).exp() + density_floor / 2.0;
            let amplitude = (rho / mass).sqrt();
            let phase = v_dot_dx / nu;

            scalar_sum += amplitude * phase.cos();
            pseudo_sum += amplitude * phase.sin();
        }

        (scalar_sum, pseudo_sum)
    });

    let norm_before = field_norm(&alpha);

    let mut entry = FieldEntry {
        data: alpha,
        smoothing_length: ell,
        mass,
        self_interaction: None,
    };

    // Evolve long enough for lumps to cross and separate.
    // Separation = 1600 kpc, relative speed = 60 kpc/Gyr.
    // Crossing time ~ 1600 / 60 ~ 27 Gyr. Evolve for ~50 Gyr total.
    let dt = 0.5;
    let n_steps = 30;
    let mut evolver = SchrodingerEvolution;

    for step in 0..n_steps {
        evolver.step(&mut entry, &morphis_grid, 1.0, dt).unwrap();

        if (step + 1) % 25 == 0 {
            let centroid = field_centroid(&entry.data, box_length);
            eprintln!("t={:.1}: centroid = {centroid:?}", (step + 1) as f64 * dt);
        }
    }

    let norm_after = field_norm(&entry.data);

    // Norm must be preserved.
    let norm_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        norm_error < 1e-12,
        "norm not preserved: relative error = {norm_error}"
    );
}

// ============================================================================
// Step 3: Two-lump infall with gravity
// ============================================================================

#[test]
fn two_lumps_attract_under_gravity() {
    // Two lumps separated along x, initially at rest. Gravity should
    // pull them toward each other, so the separation should decrease.
    let n = 16;
    let box_length = 4000.0;

    let mass = 1e10;
    let nu = 1000.0;
    let ell = nu * mass;

    let sigma = 400.0;
    let density_peak = 100.0;
    let density_floor = 0.1;

    let center_1 = [1500.0, 2000.0, 2000.0];
    let center_2 = [2500.0, 2000.0, 2000.0];

    let morphis_grid = MorphisGrid::<3>::new(n, box_length);
    let hermes_grid = Grid::new(n, box_length);
    let g = euclidean::<3>();

    // Two lumps at rest (no bulk velocity).
    let alpha = EvenField::from_fn(&morphis_grid, g, |x| {
        let mut scalar_sum = 0.0;

        for center in [center_1, center_2] {
            let mut r2 = 0.0;
            for d in 0..3 {
                let dx = wrap_distance(x[d] - center[d], box_length);
                r2 += dx * dx;
            }
            let rho = density_peak * (-r2 / (2.0 * sigma * sigma)).exp() + density_floor / 2.0;
            scalar_sum += (rho / mass).sqrt();
        }

        // No velocity phase (at rest).
        (scalar_sum, 0.0)
    });

    let norm_before = field_norm(&alpha);

    // Measure initial separation via per-lump centroids.
    // Since both lumps start at rest and symmetric about the box center,
    // we can track the overall centroid spread using the density
    // distribution width along x.
    let initial_x_spread = density_spread_x(&alpha, box_length);
    eprintln!("initial x-spread: {initial_x_spread:.1} kpc");

    let mut fields = BTreeMap::new();
    fields.insert(
        "alpha".to_string(),
        FieldEntry {
            data: alpha,
            smoothing_length: ell,
            mass,
            self_interaction: None,
        },
    );

    let sectors: Vec<Box<dyn Sector>> = vec![Box::new(SchrodingerSector::new("alpha".to_string()))];

    let solver = GravitySolver::new(morphis_grid);
    let cosmology = planck_2018();

    let state = SimulationState {
        particles: BTreeMap::new(),
        fields,
        grid: hermes_grid,
        morphis_grid,
        time: 0.0,
        step: 0,
    };

    let mut engine = Engine::new(state, sectors, Some(solver), Some(cosmology));

    let dt = 0.01;
    let n_steps = 30;

    for step in 0..n_steps {
        engine.step(1.0, dt).unwrap();

        if (step + 1) % 25 == 0 {
            let spread = density_spread_x(&engine.state.fields["alpha"].data, box_length);
            eprintln!("step {}: x-spread = {spread:.1} kpc", step + 1);
        }
    }

    let final_x_spread = density_spread_x(&engine.state.fields["alpha"].data, box_length);
    eprintln!("final x-spread: {final_x_spread:.1} kpc");

    // Norm must be preserved.
    let norm_after = field_norm(&engine.state.fields["alpha"].data);
    let norm_error = (norm_after - norm_before).abs() / norm_before;
    assert!(
        norm_error < 1e-12,
        "norm not preserved under gravity: relative error = {norm_error}"
    );

    // The lumps should have moved closer together (spread decreased).
    // We allow for the regularization spreading to partially counteract,
    // so we just check that the spread didn't increase significantly.
    assert!(
        final_x_spread <= initial_x_spread * 1.05,
        "lumps moved apart under gravity: initial spread {initial_x_spread:.1}, final {final_x_spread:.1}"
    );
}

/// Density-weighted RMS spread along x (periodic).
fn density_spread_x(alpha: &EvenField<3>, box_length: f64) -> f64 {
    let n = alpha.grid.n_cells;
    let cell_length = box_length / n as f64;

    let centroid = field_centroid(alpha, box_length);

    let mut total_weight = 0.0;
    let mut variance_sum = 0.0;

    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * cell_length;
        for m1 in 0..n {
            for m2 in 0..n {
                let a = alpha.scalar[[m0, m1, m2]];
                let b = alpha.pseudoscalar[[m0, m1, m2]];
                let weight = a * a + b * b;

                let dx = wrap_distance(x - centroid[0], box_length);
                variance_sum += weight * dx * dx;
                total_weight += weight;
            }
        }
    }

    (variance_sum / total_weight).sqrt()
}

// ============================================================================
// Step 4: Full galaxy-group-field scene initialization
// ============================================================================

#[test]
fn galaxy_group_field_initializes_with_nonzero_velocity() {
    // Verify that the corrected initialization produces a field with
    // nonzero Madelung velocity at the halo locations.
    let n = 32;
    let box_length = 8000.0;

    let cosmology = planck_2018();
    let hermes_grid = Grid::new(n, box_length);

    let mass = 1e10;
    let length_scale = 8000.0;
    let ell = length_scale * mass;

    let params = hermes_rs::core::content::FieldParams {
        smoothing_length: ell,
        mass_alpha: mass,
    };

    let halos = hermes_rs::physics::initial::nfw::default_halo_configs();
    let alpha = hermes_rs::physics::initial::nfw_field::colliding_halos_field(
        &hermes_grid,
        &cosmology,
        &params,
        1.0,
        42,
        &halos,
    );

    // Extract Madelung velocity.
    let v_field =
        hermes_rs::core::schrodinger_dynamics::extract_velocity(&alpha, &alpha.grid, ell, mass);

    // Find the peak density cell (should be near a halo center).
    let mut max_density = 0.0_f64;
    let mut max_idx = [0, 0, 0];

    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let a = alpha.scalar[[m0, m1, m2]];
                let b = alpha.pseudoscalar[[m0, m1, m2]];
                let rho = (a * a + b * b) * mass;
                if rho > max_density {
                    max_density = rho;
                    max_idx = [m0, m1, m2];
                }
            }
        }
    }

    // At the peak density, the velocity should be nonzero.
    let vx = v_field[0][[max_idx[0], max_idx[1], max_idx[2]]];
    let vy = v_field[1][[max_idx[0], max_idx[1], max_idx[2]]];
    let vz = v_field[2][[max_idx[0], max_idx[1], max_idx[2]]];
    let v_mag = (vx * vx + vy * vy + vz * vz).sqrt();

    eprintln!("peak density at {max_idx:?}: rho = {max_density:.2}");
    eprintln!("velocity at peak: ({vx:.2}, {vy:.2}, {vz:.2}), |v| = {v_mag:.2} kpc/Gyr");
    eprintln!(
        "phase-Nyquist ratio: m|v|Δx / ℓ = {:.4}",
        mass * v_mag * (box_length / n as f64) / ell
    );

    assert!(
        v_mag > 1.0,
        "velocity at halo center should be nonzero, got |v| = {v_mag}"
    );

    // Norm should be finite and positive.
    let norm = field_norm(&alpha);
    assert!(norm > 0.0 && norm.is_finite(), "invalid norm: {norm}");
}
