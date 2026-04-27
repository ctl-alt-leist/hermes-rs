//! Benchmark: timing, scaling, and timestep analysis.
//!
//! Measures per-step cost across grid sizes, reports CFL-like timestep
//! constraints, and estimates what fits in a 1-minute wall-clock budget.
//!
//! Run with:
//!   cargo run --example benchmark --release

use std::time::Instant;

use hermes_rs::config::build_configuration;
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::initial::zeldovich_init;
use hermes_rs::physics::integrator::{midpoint, step_kdk};
use hermes_rs::physics::poisson::PoissonSolver;
use hermes_rs::physics::simulation::Simulation;

fn main() {
    println!("Hermes Benchmark");
    println!("================\n");

    // ========================================================================
    // 1. Per-step timing across grid sizes
    // ========================================================================
    println!("1. Per-step timing (5 steps each, release mode)");
    println!(
        "   {:>6}  {:>10}  {:>12}  {:>10}",
        "N", "Particles", "ms/step", "steps/min"
    );
    println!("   {}", "-".repeat(48));

    for n in [8, 16, 32, 48, 64] {
        let time_per_step = benchmark_steps(n, 5);
        let steps_per_minute = 60_000.0 / time_per_step;
        println!(
            "   {:>6}  {:>10}  {:>10.1} ms  {:>10.0}",
            n,
            n * n * n,
            time_per_step,
            steps_per_minute
        );
    }

    // ========================================================================
    // 2. Timestep analysis
    // ========================================================================
    println!("\n2. Timestep constraints");
    println!();

    let cosmology = planck_2018();

    // The dynamical time at scale factor a is t_dyn ~ 1 / sqrt(G ρ(a)).
    // For a PM code, the relevant constraint is that particles shouldn't
    // cross more than one cell per step (CFL-like condition).
    //
    // The step in scale factor Δa maps to a time interval Δt = Δa / (a H(a)).
    // The displacement per step is ~ v × Δt, where v ~ p / (m a²).
    // The CFL condition is: v × Δt < h (cell size).
    //
    // We can also express this as: Δa < h × a² × H(a) × m / p_max
    //
    // For a typical simulation at z=0 (a=1), peculiar velocities are
    // ~300 km/s ≈ 307 kpc/Gyr. With a 50 Mpc box and N=32 cells,
    // h = 1562 kpc, so CFL gives Δa < h × H(a) / v_pec ≈ 0.5.
    // This is very loose — the dynamical time is the real constraint.

    println!("   Dynamical time t_dyn = 1/sqrt(4πG ρ̄) at different redshifts:");
    println!("   {:>6}  {:>8}  {:>12}", "z", "a", "t_dyn (Gyr)");
    println!("   {}", "-".repeat(32));

    for &a in &[0.02, 0.05, 0.1, 0.2, 0.5, 1.0] {
        let redshift = 1.0 / a - 1.0;
        let density = cosmology.density_matter();
        let g = hermes_rs::physics::constants::G;
        let dynamical_time = 1.0 / (4.0 * std::f64::consts::PI * g * density / (a * a * a)).sqrt();
        println!("   {:>6.1}  {:>8.3}  {:>10.4}", redshift, a, dynamical_time);
    }

    println!();
    println!("   For accurate integration, Δt should be < 0.1 × t_dyn.");
    println!("   With log-a stepping over 100 steps from a=0.02 to a=1.0:");

    let n_steps = 100;
    let log_start = 0.02_f64.ln();
    let log_end = 1.0_f64.ln();
    let d_log = (log_end - log_start) / n_steps as f64;

    println!("   Δ(ln a) = {:.4}", d_log);
    println!();
    println!(
        "   {:>6}  {:>8}  {:>10}  {:>12}  {:>10}",
        "step", "a", "Δa", "Δt (Gyr)", "Δt/t_dyn"
    );
    println!("   {}", "-".repeat(56));

    for &step in &[0, 10, 25, 50, 75, 99] {
        let a = (log_start + step as f64 * d_log).exp();
        let a_next = (log_start + (step + 1) as f64 * d_log).exp();
        let da = a_next - a;
        let dt = da / (a * cosmology.hubble_parameter(a));
        let density = cosmology.density_matter();
        let g = hermes_rs::physics::constants::G;
        let dynamical_time = 1.0 / (4.0 * std::f64::consts::PI * g * density / (a * a * a)).sqrt();
        let ratio = dt / dynamical_time;
        println!(
            "   {:>6}  {:>8.4}  {:>10.6}  {:>10.4}  {:>10.4}",
            step, a, da, dt, ratio
        );
    }

    // ========================================================================
    // 3. What fits in 1 minute
    // ========================================================================
    println!("\n3. What fits in 1 minute (wall clock)");
    println!();

    for n in [16, 32, 48, 64] {
        let time_per_step = benchmark_steps(n, 3);
        let steps_per_minute = (60_000.0 / time_per_step) as usize;

        // With log-a stepping, what range of a can we cover?
        // Use 100 steps as the reference for z=49 → z=0.
        let fraction = steps_per_minute as f64 / 100.0;
        let a_final = if fraction >= 1.0 {
            1.0
        } else {
            (log_start + fraction * (log_end - log_start)).exp()
        };
        let z_final = 1.0 / a_final - 1.0;

        // Cosmic time spanned
        let cosmic_time_start = cosmology.cosmic_time(0.02);
        let cosmic_time_end = cosmology.cosmic_time(a_final.min(1.0));
        let cosmic_time_span = cosmic_time_end - cosmic_time_start;

        println!(
            "   N={n:>2}³ ({:>6} particles): {:.1} ms/step → {steps_per_minute} steps/min",
            n * n * n,
            time_per_step,
        );

        if fraction >= 1.0 {
            let n_full_runs = steps_per_minute / 100;
            println!(
                "          → {n_full_runs} full runs (z=49→0) in 1 min, {:.1} Gyr simulated",
                cosmic_time_span
            );
        } else {
            println!(
                "          → reaches z={z_final:.1} (a={a_final:.3}), {:.2} Gyr simulated",
                cosmic_time_span
            );
        }
    }

    println!();

    // ========================================================================
    // 4. Velocity / CFL analysis from a real run
    // ========================================================================
    println!("4. Particle velocity statistics from a 32³ run");
    println!();

    let grid = Grid::new(32, 50_000.0);
    let cosmology = planck_2018();
    let _solver = PoissonSolver::new(&grid);
    let particles = zeldovich_init(32, &grid, &cosmology, 0.02, 42).unwrap();

    // Initial velocities: v = p / (m a²)
    let scale_factor = 0.02;
    let mass = particles.mass_particle;
    let h = grid.cell_length;

    let velocity_max_initial = (0..particles.count())
        .map(|p| {
            let mom = particles.momentum_of(p);
            mom.norm() / (mass * scale_factor * scale_factor)
        })
        .fold(0.0_f64, f64::max);

    println!(
        "   At z=49 (a=0.02): v_max = {:.2} kpc/Gyr = {:.2} km/s",
        velocity_max_initial,
        velocity_max_initial * hermes_rs::physics::constants::KPC_GYR_TO_KMS
    );
    println!("   Cell size h = {:.1} kpc", h);
    println!(
        "   CFL: v_max × Δt < h → Δt < {:.4} Gyr",
        h / velocity_max_initial
    );

    // Run 50 steps and check final velocities
    let config_val: toml::Value = toml::from_str(
        r#"
        [simulation]
        n_cells = 32
        n_particles = 32
        box_length = 50000.0
        [time]
        scale_factor_initial = 0.02
        scale_factor_final = 1.0
        n_steps = 100
        stepping = "log_a"
        [output]
        snapshot_interval = 100
    "#,
    )
    .unwrap();
    let config = build_configuration(None, Some(&config_val)).unwrap();
    let mut sim = Simulation::from_config(config, 42).unwrap();
    sim.run(&mut []).unwrap();

    let velocity_max_final = (0..sim.particles.count())
        .map(|p| {
            let mom = sim.particles.momentum_of(p);
            mom.norm() / (mass * 1.0 * 1.0)
        })
        .fold(0.0_f64, f64::max);

    println!(
        "   At z=0  (a=1.0):  v_max = {:.2} kpc/Gyr = {:.2} km/s",
        velocity_max_final,
        velocity_max_final * hermes_rs::physics::constants::KPC_GYR_TO_KMS
    );
    println!(
        "   CFL: v_max × Δt < h → Δt < {:.4} Gyr",
        h / velocity_max_final
    );
}

fn benchmark_steps(n: usize, n_steps: usize) -> f64 {
    let grid = Grid::new(n, 50_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);
    let mut particles = zeldovich_init(n, &grid, &cosmology, 0.02, 42).unwrap();

    let log_start = 0.02_f64.ln();
    let log_end = 1.0_f64.ln();
    let d_log = (log_end - log_start) / 100.0;

    let start = Instant::now();
    let mut forces_prev = None;

    for step in 0..n_steps {
        let a_prev = (log_start + step as f64 * d_log).exp();
        let a_next = (log_start + (step + 1) as f64 * d_log).exp();
        let a_mid = midpoint(a_prev, a_next);

        let forces = step_kdk(
            &mut particles,
            &mut solver,
            &grid,
            &cosmology,
            a_prev,
            a_mid,
            a_next,
            forces_prev.as_ref(),
        )
        .unwrap();

        forces_prev = Some(forces);
    }

    let elapsed = start.elapsed();

    elapsed.as_secs_f64() * 1000.0 / n_steps as f64
}
