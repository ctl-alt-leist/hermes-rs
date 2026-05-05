//! Third-pass diagnostics: linear growth deficit isolation.
//!
//! Traces the gravitational coupling step-by-step in the linear regime
//! where the analytic answer is known exactly.
//!
//! Run with:
//!   cargo run --example diagnostics_pass3 --release

use std::f64::consts::PI;
use std::fmt::Write;

use morphis::even_field::EvenField;
use morphis::field::Field;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric;
use ndarray::{Array3, IxDyn};
use num_complex::Complex64;

use hermes_rs::core::content::{Content, FieldParams, FieldState};
use hermes_rs::core::dynamics::Dynamics;
use hermes_rs::core::schrodinger_dynamics::{SchrodingerPoissonDynamics, kinetic_step};
use hermes_rs::engine::coupling::poisson::PoissonGravity;
use hermes_rs::physics::constants::G as GRAV;
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid as HermesGrid;

fn main() {
    let mut report = String::with_capacity(8192);
    writeln!(report, "# Third-Pass Diagnostic: Linear Growth Isolation\n").unwrap();
    writeln!(
        report,
        "Generated: {}\n",
        chrono::Local::now().format("%Y-%m-%d %H:%M:%S")
    )
    .unwrap();

    let n = 64_usize;
    let box_length = 10000.0;
    let ell_over_m = 7500.0;
    let mass = 1e10;
    let ell = ell_over_m * mass;
    let cosmology = planck_2018();
    let rho_bar = cosmology.density_matter();
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let dx = box_length / n as f64;
    let n3 = (n * n * n) as f64;

    let a_init = 0.5;
    let a_final = 1.0;
    let n_total_steps = 2000;
    let da = (a_final - a_init) / n_total_steps as f64;
    let n_test_steps = 10;

    let amp_init = 1e-3;
    let mode_idx: [usize; 3] = [1, 0, 0];
    let k0 = 2.0 * PI / box_length;
    let k2 = k0 * k0;

    let k_j = (16.0 * PI * GRAV * rho_bar / (ell_over_m * ell_over_m)).powf(0.25);
    let lambda_j = 2.0 * PI / k_j;

    writeln!(report, "## Parameters\n").unwrap();
    writeln!(report, "| Parameter | Value |").unwrap();
    writeln!(report, "|---|---|").unwrap();
    writeln!(report, "| rho_bar | {:.6e} Msun/kpc^3 |", rho_bar).unwrap();
    writeln!(
        report,
        "| rho_bar (analytic: Omega_m * rho_crit) | {:.6e} |",
        cosmology.omega_m * cosmology.density_critical(1.0)
    )
    .unwrap();
    writeln!(report, "| ell/m | {:.1} kpc^2/Gyr |", ell_over_m).unwrap();
    writeln!(report, "| k0 | {:.6e} kpc^-1 |", k0).unwrap();
    writeln!(report, "| k_J | {:.6e} kpc^-1 |", k_j).unwrap();
    writeln!(report, "| k0 / k_J | {:.4} |", k0 / k_j).unwrap();
    writeln!(report, "| lambda_J | {:.0} kpc |", lambda_j).unwrap();
    writeln!(report, "| a_init | {:.2} |", a_init).unwrap();
    writeln!(report, "| A_init | {:.1e} |", amp_init).unwrap();
    writeln!(report, "| da | {:.6e} |", da).unwrap();
    writeln!(
        report,
        "| D+(a_init) | {:.6} |",
        cosmology.growth_factor(a_init)
    )
    .unwrap();

    // ========================================================================
    // Part 1: Gravity OFF — verify density preservation
    // ========================================================================

    writeln!(report, "\n## Part 1: Gravity OFF (kinetic only)\n").unwrap();
    writeln!(
        report,
        "Verifying that a single-mode wavefunction with Zeldovich velocity"
    )
    .unwrap();
    writeln!(
        report,
        "preserves its density under free-particle evolution.\n"
    )
    .unwrap();

    // Construct wavefunction: rho = rho_bar (1 + A cos(k0 x)), v from Zeldovich.
    // Zeldovich velocity for a growing mode: v = -H a f delta / k
    // For simplicity, use zero velocity first to test pure kinetic dispersion.
    {
        let mut psi = EvenField::from_fn(&grid, g, |x| {
            let rho = rho_bar * (1.0 + amp_init * (k0 * x[0]).cos());
            let amp = (rho / mass).sqrt();
            (amp, 0.0) // zero velocity
        });

        // Extract initial amplitude at the mode.
        let amp_before = extract_mode_amplitude(&psi, mass, rho_bar, &mode_idx, n);

        writeln!(
            report,
            "| Step | a | mean |alpha|^2 | A_sim | drift from A_init |"
        )
        .unwrap();
        writeln!(report, "|---|---|---|---|---|").unwrap();

        let mean_mod_sq_0 = mean_mod_squared(&psi, n);
        writeln!(
            report,
            "| 0 | {:.4} | {:.10e} | {:.6e} | {:.4e} |",
            a_init,
            mean_mod_sq_0,
            amp_before,
            (amp_before - amp_init).abs()
        )
        .unwrap();

        for step in 1..=n_test_steps {
            let a = a_init + (step - 1) as f64 * da;
            let dt = da / (a * cosmology.hubble_parameter(a));
            // Full kinetic step (not half-step) to match what the integrator does in KDK.
            kinetic_step(&mut psi, &grid, ell, mass, a, dt);

            let mean_ms = mean_mod_squared(&psi, n);
            let amp_now = extract_mode_amplitude(&psi, mass, rho_bar, &mode_idx, n);

            writeln!(
                report,
                "| {step} | {:.4} | {:.10e} | {:.6e} | {:.4e} |",
                a + da,
                mean_ms,
                amp_now,
                (amp_now - amp_init).abs()
            )
            .unwrap();
        }
    }

    // ========================================================================
    // Part 2: Gravity ON — step-by-step trace
    // ========================================================================

    writeln!(report, "\n## Part 2: Gravity ON (full integrator)\n").unwrap();

    let params = FieldParams {
        smoothing_length: ell,
        mass_alpha: mass,
    };

    // Construct with zero velocity (simplest case).
    let psi = EvenField::from_fn(&grid, g, |x| {
        let rho = rho_bar * (1.0 + amp_init * (k0 * x[0]).cos());
        let amp = (rho / mass).sqrt();
        (amp, 0.0)
    });

    let field_state = FieldState {
        grid: grid.clone(),
        alpha: Some(psi),
        beta: None,
        gamma: None,
        params,
    };
    let mut content = Content::Fields(field_state);
    let hermes_grid = HermesGrid::new(n, box_length);
    let mut dynamics = SchrodingerPoissonDynamics::new(PoissonGravity::new(hermes_grid));

    let growth_init = cosmology.growth_factor(a_init);

    writeln!(report, "| Step | a | dt (Gyr) | mean |alpha|^2 | A_sim | A_lin | A_sim/A_lin | Phi_hat actual | Phi_hat analytic | Phi ratio | theta_kin | theta_pot |").unwrap();
    writeln!(report, "|---|---|---|---|---|---|---|---|---|---|---|---|").unwrap();

    // Print initial state.
    {
        let psi_ref = content.fields().unwrap().alpha.as_ref().unwrap();
        let mean_ms = mean_mod_squared(psi_ref, n);
        let amp_now = extract_mode_amplitude(psi_ref, mass, rho_bar, &mode_idx, n);

        writeln!(
            report,
            "| 0 | {:.4} | — | {:.10e} | {:.6e} | {:.6e} | {:.6} | — | — | — | — | — |",
            a_init,
            mean_ms,
            amp_now,
            amp_init,
            amp_now / amp_init
        )
        .unwrap();
    }

    for step in 1..=n_test_steps {
        let a_prev = a_init + (step - 1) as f64 * da;
        let a_next = a_prev + da;
        let a_mid = (a_prev + a_next) / 2.0;
        let dt = da / (a_mid * cosmology.hubble_parameter(a_mid));

        // Before stepping, extract the Poisson source to see what goes in.
        let psi_ref = content.fields().unwrap().alpha.as_ref().unwrap();
        let delta_hat_coeff = extract_delta_hat(psi_ref, mass, rho_bar, &mode_idx, n);

        // Analytic Phi_hat at this mode.
        let prefactor = 4.0 * PI * GRAV * rho_bar * a_mid * a_mid;
        let phi_hat_analytic = -prefactor * delta_hat_coeff / k2;

        // Actually run the Poisson solver to get the real Phi_hat.
        let rho_field = psi_ref.density(mass);
        let mut delta_grid = Array3::<f64>::zeros((n, n, n));
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    delta_grid[[m0, m1, m2]] =
                        rho_field.at(&[m0, m1, m2]).component(&[]) / rho_bar - 1.0;
                }
            }
        }
        let delta_field = Field::scalar_field(&grid, g, |x| {
            let m0 = ((x[0] / dx) as usize).min(n - 1);
            let m1 = ((x[1] / dx) as usize).min(n - 1);
            let m2 = ((x[2] / dx) as usize).min(n - 1);
            delta_grid[[m0, m1, m2]]
        });
        let source = &delta_field * prefactor;
        let phi = source.laplacian_inverse();

        // Extract Phi at the mode from the actual solver output.
        let mut phi_grid = Array3::<f64>::zeros((n, n, n));
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    phi_grid[[m0, m1, m2]] = phi.at(&[m0, m1, m2]).component(&[]);
                }
            }
        }
        let phi_hat = fft_3d_static(&phi_grid, n);
        let phi_hat_actual = phi_hat[[mode_idx[0], mode_idx[1], mode_idx[2]]];

        let phi_ratio = phi_hat_actual.norm() / phi_hat_analytic.norm();

        // Phase rotations at this mode.
        let theta_kin = ell_over_m * k2 * dt / (2.0 * a_mid * a_mid);
        let theta_pot = mass * phi_hat_actual.norm() / n3 * dt / ell;

        // Step the integrator.
        dynamics
            .step(&mut content, &cosmology, a_prev, a_next)
            .unwrap();

        // Extract state after step.
        let psi_after = content.fields().unwrap().alpha.as_ref().unwrap();
        let mean_ms = mean_mod_squared(psi_after, n);
        let amp_now = extract_mode_amplitude(psi_after, mass, rho_bar, &mode_idx, n);
        let amp_lin = amp_init * cosmology.growth_factor(a_next) / growth_init;
        let ratio = amp_now / amp_lin;

        writeln!(report,
            "| {step} | {:.4} | {:.5e} | {:.10e} | {:.6e} | {:.6e} | {:.6} | {:.6e} | {:.6e} | {:.8} | {:.4e} | {:.4e} |",
            a_next, dt, mean_ms, amp_now, amp_lin, ratio,
            phi_hat_actual.norm(), phi_hat_analytic.norm(), phi_ratio,
            theta_kin, theta_pot
        ).unwrap();
    }

    // ========================================================================
    // Part 3: Direct Poisson coefficient check at the mode
    // ========================================================================

    writeln!(report, "\n## Part 3: Poisson Source Verification\n").unwrap();
    writeln!(
        report,
        "Checking whether the density that enters the Poisson solver matches"
    )
    .unwrap();
    writeln!(report, "the density extracted from the wavefunction.\n").unwrap();

    // Fresh wavefunction.
    let psi_fresh = EvenField::from_fn(&grid, g, |x| {
        let rho = rho_bar * (1.0 + amp_init * (k0 * x[0]).cos());
        ((rho / mass).sqrt(), 0.0)
    });

    // Method A: Extract density via psi.density(mass), then compute delta.
    let rho_a = psi_fresh.density(mass);
    let mut delta_a = Array3::<f64>::zeros((n, n, n));
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                delta_a[[m0, m1, m2]] = rho_a.at(&[m0, m1, m2]).component(&[]) / rho_bar - 1.0;
            }
        }
    }
    let delta_hat_a = fft_3d_static(&delta_a, n);
    let coeff_a = delta_hat_a[[1, 0, 0]];

    // Method B: Compute delta directly from |alpha|^2.
    let mut delta_b = Array3::<f64>::zeros((n, n, n));
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let s = psi_fresh.scalar[IxDyn(&[m0, m1, m2])];
                let ps = psi_fresh.pseudoscalar[IxDyn(&[m0, m1, m2])];
                let mod_sq = s * s + ps * ps;
                delta_b[[m0, m1, m2]] = mod_sq / (rho_bar / mass) - 1.0;
            }
        }
    }
    let delta_hat_b = fft_3d_static(&delta_b, n);
    let coeff_b = delta_hat_b[[1, 0, 0]];

    // Method C: What the analytic value should be.
    // delta_hat at (1,0,0) for A cos(k0 x) with un-normalized FFT = A * N^3 / 2
    let coeff_analytic = amp_init * n3 / 2.0;

    writeln!(
        report,
        "| Method | delta_hat at (1,0,0) | ratio to analytic |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|").unwrap();
    writeln!(
        report,
        "| A (via density()) | ({:.6e}, {:.6e}), mag {:.6e} | {:.8} |",
        coeff_a.re,
        coeff_a.im,
        coeff_a.norm(),
        coeff_a.norm() / coeff_analytic
    )
    .unwrap();
    writeln!(
        report,
        "| B (via |alpha|^2) | ({:.6e}, {:.6e}), mag {:.6e} | {:.8} |",
        coeff_b.re,
        coeff_b.im,
        coeff_b.norm(),
        coeff_b.norm() / coeff_analytic
    )
    .unwrap();
    writeln!(
        report,
        "| Analytic (A N^3 / 2) | {:.6e} | 1.0 |",
        coeff_analytic
    )
    .unwrap();

    // Also check: is the Poisson solver receiving the same delta as what we compute?
    // The potential step constructs delta as: rho / rho_bar - 1, using morphis Field arithmetic.
    // Let's replicate that exactly.
    writeln!(
        report,
        "\n### Replicating potential_step's delta construction\n"
    )
    .unwrap();

    let ones = Field::scalar_field(&grid, g, |_| 1.0);
    let mut overdensity = &rho_a * (1.0 / rho_bar);
    overdensity = &overdensity - &ones;

    let mut delta_c = Array3::<f64>::zeros((n, n, n));
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                delta_c[[m0, m1, m2]] = overdensity.at(&[m0, m1, m2]).component(&[]);
            }
        }
    }
    let delta_hat_c = fft_3d_static(&delta_c, n);
    let coeff_c = delta_hat_c[[1, 0, 0]];

    writeln!(
        report,
        "| C (via Field arithmetic) | ({:.6e}, {:.6e}), mag {:.6e} | {:.8} |",
        coeff_c.re,
        coeff_c.im,
        coeff_c.norm(),
        coeff_c.norm() / coeff_analytic
    )
    .unwrap();

    // ========================================================================
    // Part 4: rho_bar consistency
    // ========================================================================

    writeln!(report, "\n## Part 4: rho_bar Consistency\n").unwrap();

    let rho_bar_code = cosmology.density_matter();
    let rho_crit = cosmology.density_critical(1.0);
    let rho_bar_manual = cosmology.omega_m * rho_crit;

    writeln!(report, "| Source | Value |").unwrap();
    writeln!(report, "|---|---|").unwrap();
    writeln!(
        report,
        "| cosmology.density_matter() | {:.10e} |",
        rho_bar_code
    )
    .unwrap();
    writeln!(
        report,
        "| Omega_m * density_critical(1) | {:.10e} |",
        rho_bar_manual
    )
    .unwrap();
    writeln!(
        report,
        "| Difference | {:.4e} |",
        (rho_bar_code - rho_bar_manual).abs()
    )
    .unwrap();

    // Mean |alpha|^2 * mass vs rho_bar.
    let mean_rho = mean_mod_squared(&psi_fresh, n) * mass;
    writeln!(report, "| mean(|alpha|^2) * m | {:.10e} |", mean_rho).unwrap();
    writeln!(
        report,
        "| mean_rho / rho_bar | {:.10} |",
        mean_rho / rho_bar
    )
    .unwrap();

    // ========================================================================
    // Write report
    // ========================================================================

    writeln!(report, "\n---\n## Summary\n").unwrap();
    writeln!(report, "See sections above. Key numbers to check:").unwrap();
    writeln!(
        report,
        "- Part 1: Does density drift with kinetic-only evolution?"
    )
    .unwrap();
    writeln!(
        report,
        "- Part 2: Does A_sim/A_lin converge to 1.0 or drift away?"
    )
    .unwrap();
    writeln!(report, "- Part 2: Is Phi ratio consistently 1.0?").unwrap();
    writeln!(
        report,
        "- Part 3: Do all methods give the same delta_hat coefficient?"
    )
    .unwrap();
    writeln!(report, "- Part 4: Is rho_bar self-consistent?").unwrap();

    std::fs::create_dir_all("data").ok();
    std::fs::write("data/diagnostics-pass3-report.md", &report).unwrap();
    println!("Report written to data/diagnostics-pass3-report.md");
}

// ============================================================================
// Helpers
// ============================================================================

fn mean_mod_squared(psi: &EvenField<3>, n: usize) -> f64 {
    let mut sum = 0.0;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let s = psi.scalar[IxDyn(&[m0, m1, m2])];
                let ps = psi.pseudoscalar[IxDyn(&[m0, m1, m2])];
                sum += s * s + ps * ps;
            }
        }
    }
    sum / (n * n * n) as f64
}

fn extract_mode_amplitude(
    psi: &EvenField<3>,
    mass: f64,
    rho_bar: f64,
    mode_idx: &[usize; 3],
    n: usize,
) -> f64 {
    let n3 = (n * n * n) as f64;

    // Compute delta = rho / rho_bar - 1 on the grid.
    let mut delta = Array3::<f64>::zeros((n, n, n));
    let target = rho_bar / mass;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let s = psi.scalar[IxDyn(&[m0, m1, m2])];
                let ps = psi.pseudoscalar[IxDyn(&[m0, m1, m2])];
                let mod_sq = s * s + ps * ps;
                delta[[m0, m1, m2]] = mod_sq / target - 1.0;
            }
        }
    }

    let delta_hat = fft_3d_static(&delta, n);
    // For a real cosine, the amplitude is 2 * |coeff| / N^3.
    delta_hat[[mode_idx[0], mode_idx[1], mode_idx[2]]].norm() * 2.0 / n3
}

fn extract_delta_hat(
    psi: &EvenField<3>,
    mass: f64,
    rho_bar: f64,
    mode_idx: &[usize; 3],
    n: usize,
) -> Complex64 {
    let mut delta = Array3::<f64>::zeros((n, n, n));
    let target = rho_bar / mass;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let s = psi.scalar[IxDyn(&[m0, m1, m2])];
                let ps = psi.pseudoscalar[IxDyn(&[m0, m1, m2])];
                delta[[m0, m1, m2]] = (s * s + ps * ps) / target - 1.0;
            }
        }
    }
    let delta_hat = fft_3d_static(&delta, n);
    delta_hat[[mode_idx[0], mode_idx[1], mode_idx[2]]]
}

fn fft_3d_static(data: &Array3<f64>, n: usize) -> Array3<Complex64> {
    let n_c = n / 2 + 1;
    let handler_r2c = ndrustfft::R2cFftHandler::new(n);
    let handler_c2c_1 = ndrustfft::FftHandler::new(n);
    let handler_c2c_0 = ndrustfft::FftHandler::new(n);

    let mut complex = Array3::<Complex64>::zeros((n, n, n_c));
    ndrustfft::ndfft_r2c(data, &mut complex, &handler_r2c, 2);
    let mut scratch = complex.clone();
    ndrustfft::ndfft(&complex, &mut scratch, &handler_c2c_1, 1);
    complex.assign(&scratch);
    ndrustfft::ndfft(&complex, &mut scratch, &handler_c2c_0, 0);
    complex.assign(&scratch);
    complex
}
