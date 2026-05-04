//! Diagnostic test battery for the Schrodinger-Poisson integrator.
//!
//! Runs 8 independent tests isolating different failure modes,
//! produces a single markdown report at data/diagnostics-report.md.
//!
//! Run with:
//!   cargo run --example diagnostics --release

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
use hermes_rs::physics::cosmology::{Cosmology, planck_2018};
use hermes_rs::physics::grid::Grid as HermesGrid;
use hermes_rs::physics::spectral::{fft_3d_dyn as fft_3d, ifft_3d_dyn as ifft_3d};

// ============================================================================
// Parameters
// ============================================================================

struct DiagParams {
    n: usize,
    box_length: f64,
    ell_over_m: f64,
    mass: f64,
    ell: f64,
    a_init: f64,
    a_final: f64,
    cosmology: Cosmology,
    rho_bar: f64,
    morphis_grid: MorphisGrid<3>,
}

impl DiagParams {
    fn new() -> Self {
        let n = 64;
        let box_length = 10000.0;
        let ell_over_m = 7500.0;
        let mass = 1e10;
        let ell = ell_over_m * mass;
        let a_init = 0.33;
        let a_final = 1.0;
        let cosmology = planck_2018();
        let rho_bar = cosmology.density_matter();
        let morphis_grid = MorphisGrid::<3>::new(n, box_length);

        Self {
            n,
            box_length,
            ell_over_m,
            mass,
            ell,
            a_init,
            a_final,
            cosmology,
            rho_bar,
            morphis_grid,
        }
    }

    fn dx(&self) -> f64 {
        self.box_length / self.n as f64
    }
    fn k_fundamental(&self) -> f64 {
        2.0 * PI / self.box_length
    }
    fn k_nyquist(&self) -> f64 {
        PI * self.n as f64 / self.box_length
    }

    fn dt_at(&self, a: f64) -> f64 {
        let da = (self.a_final - self.a_init) / 2000.0;
        da / (a * self.cosmology.hubble_parameter(a))
    }

    fn lambda_jeans(&self) -> f64 {
        let k_j4 = 16.0 * PI * GRAV * self.rho_bar / (self.ell_over_m * self.ell_over_m);
        2.0 * PI / k_j4.powf(0.25)
    }
}

fn wavenum(m: usize, n: usize, box_length: f64) -> f64 {
    let freq = if m <= n / 2 {
        m as f64
    } else {
        m as f64 - n as f64
    };
    2.0 * PI * freq / box_length
}

// ============================================================================
// Test 1: Dimensional Sanity Check
// ============================================================================

fn test_1(report: &mut String, p: &DiagParams) {
    writeln!(report, "## Test 1: Dimensional Sanity Check\n").unwrap();

    let dt = p.dt_at(p.a_init);
    let sigma_v = 307.0;
    let lambda_db = 2.0 * PI * p.ell_over_m / sigma_v;
    let lambda_j = p.lambda_jeans();
    let courant = p.ell_over_m * p.k_nyquist().powi(2) * dt / (2.0 * p.a_init * p.a_init);

    writeln!(report, "### Parameters\n").unwrap();
    writeln!(report, "| Parameter | Value | Units |").unwrap();
    writeln!(report, "|-----------|-------|-------|").unwrap();
    writeln!(report, "| G | {:.4e} | kpc^3 Msun^-1 Gyr^-2 |", GRAV).unwrap();
    writeln!(
        report,
        "| H0 | {:.6} | Gyr^-1 |",
        p.cosmology.hubble_constant()
    )
    .unwrap();
    writeln!(report, "| rho_bar | {:.4} | Msun kpc^-3 |", p.rho_bar).unwrap();
    writeln!(report, "| ell/m | {:.1} | kpc^2 Gyr^-1 |", p.ell_over_m).unwrap();
    writeln!(report, "| ell | {:.4e} | kpc^2 Msun Gyr^-1 |", p.ell).unwrap();
    writeln!(report, "| m | {:.1e} | Msun |", p.mass).unwrap();
    writeln!(report, "| box_length | {:.0} | kpc |", p.box_length).unwrap();
    writeln!(report, "| n_cells | {} | |", p.n).unwrap();
    writeln!(report, "| dx | {:.2} | kpc |", p.dx()).unwrap();
    writeln!(report, "| a_init | {:.2} | |", p.a_init).unwrap();
    writeln!(report, "| dt(a_init) | {:.5} | Gyr |", dt).unwrap();
    writeln!(report, "| k_f | {:.6e} | kpc^-1 |", p.k_fundamental()).unwrap();
    writeln!(report, "| k_Ny | {:.6e} | kpc^-1 |", p.k_nyquist()).unwrap();
    writeln!(report, "| lambda_dB(300 km/s) | {:.1} | kpc |", lambda_db).unwrap();
    writeln!(
        report,
        "| lambda_J | {:.0} kpc ({:.2} Mpc) | |",
        lambda_j,
        lambda_j / 1000.0
    )
    .unwrap();
    writeln!(
        report,
        "| lambda_J / dx | {:.1} | cells |",
        lambda_j / p.dx()
    )
    .unwrap();
    writeln!(
        report,
        "| Courant (k_Ny) | {:.6} | (must be < pi) |",
        courant
    )
    .unwrap();

    writeln!(report, "\n### Phase Rotations at a = {:.2}\n", p.a_init).unwrap();
    writeln!(
        report,
        "| k | theta_kin (rad) | theta_pot (delta=1, rad) | pot/kin |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|").unwrap();

    let a2 = p.a_init * p.a_init;
    for (label, k) in [
        ("k_f", p.k_fundamental()),
        ("3 k_f", 3.0 * p.k_fundamental()),
        ("k_Ny", p.k_nyquist()),
    ] {
        let k2 = k * k;
        let theta_kin = p.ell_over_m * k2 * dt / (2.0 * a2);
        let phi_amp = 4.0 * PI * GRAV * p.rho_bar * a2 / k2;
        let theta_pot = phi_amp * dt / p.ell_over_m;
        writeln!(
            report,
            "| {label} | {theta_kin:.6e} | {theta_pot:.6e} | {:.4} |",
            theta_pot / theta_kin
        )
        .unwrap();
    }

    let pass = courant < PI;
    writeln!(
        report,
        "\n**TEST 1: {}**\n",
        if pass {
            "PASSED"
        } else {
            "FAILED — CFL violation"
        }
    )
    .unwrap();
}

// ============================================================================
// Test 2: Poisson Solver Verification
// ============================================================================

fn test_2(report: &mut String, p: &DiagParams) {
    writeln!(report, "## Test 2: Poisson Solver Verification\n").unwrap();

    let g = metric::euclidean::<3>();
    let a = p.a_init;
    let prefactor = 4.0 * PI * GRAV * p.rho_bar * a * a;
    let amp = 0.1;
    let n = p.n;

    writeln!(
        report,
        "| Mode | k (kpc^-1) | max |error| | num_max/analytic_max | sign at delta_max |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|---|").unwrap();

    let mut all_pass = true;

    for (label, mode_n) in [
        ("low (2 k_f)", 2_usize),
        ("mid (N/4 k_f)", n / 4),
        ("high (N/2 k_f)", n / 2),
    ] {
        let k0 = 2.0 * PI * mode_n as f64 / p.box_length;
        let k2 = k0 * k0;
        let phi_analytic_max = prefactor * amp / k2;

        let delta = Field::scalar_field(&p.morphis_grid, g, |x| amp * (k0 * x[0]).cos());
        let source = &delta * prefactor;
        let phi_num = source.laplacian_inverse();

        let mut max_error = 0.0_f64;
        let mut phi_num_at_origin = 0.0_f64;

        for m0 in 0..n {
            let x = (m0 as f64 + 0.5) * p.dx();
            let phi_analytic = -phi_analytic_max * (k0 * x).cos();
            let phi_numerical = phi_num.at(&[m0, 0, 0]).component(&[]);
            let error = (phi_numerical - phi_analytic).abs();
            if error > max_error {
                max_error = error;
            }
            if m0 == 0 {
                phi_num_at_origin = phi_numerical;
            }
        }

        let phi_analytic_at_origin = -phi_analytic_max * (k0 * 0.5 * p.dx()).cos();
        let ratio = phi_num_at_origin / phi_analytic_at_origin;

        let sign_str = if phi_num_at_origin < 0.0 {
            "negative (correct)"
        } else {
            "POSITIVE (WRONG)"
        };
        if phi_num_at_origin >= 0.0 || max_error > 1e-6 {
            all_pass = false;
        }

        writeln!(
            report,
            "| {label} | {k0:.6e} | {max_error:.4e} | {ratio:.8} | {sign_str} |"
        )
        .unwrap();
    }

    // Gaussian blob test.
    writeln!(report, "\n### Gaussian Blob Test\n").unwrap();
    let sigma = 3.0 * p.dx();
    let center = p.box_length / 2.0;
    let blob = Field::scalar_field(&p.morphis_grid, g, |x| {
        let r2 = (x[0] - center).powi(2) + (x[1] - center).powi(2) + (x[2] - center).powi(2);
        amp * (-r2 / (2.0 * sigma * sigma)).exp()
    });
    let phi_blob = (&blob * prefactor).laplacian_inverse();

    let center_idx = n / 2;
    let phi_at_center = phi_blob
        .at(&[center_idx, center_idx, center_idx])
        .component(&[]);

    let mut phi_min = f64::MAX;
    let mut min_idx = [0_usize; 3];
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let val = phi_blob.at(&[m0, m1, m2]).component(&[]);
                if val < phi_min {
                    phi_min = val;
                    min_idx = [m0, m1, m2];
                }
            }
        }
    }

    writeln!(report, "- Phi at center: {phi_at_center:.6e}").unwrap();
    writeln!(
        report,
        "- Phi minimum: {phi_min:.6e} at [{},{},{}]",
        min_idx[0], min_idx[1], min_idx[2]
    )
    .unwrap();
    writeln!(
        report,
        "- Sign: {}",
        if phi_at_center < 0.0 {
            "negative (correct)"
        } else {
            "POSITIVE (WRONG)"
        }
    )
    .unwrap();
    if phi_at_center >= 0.0 {
        all_pass = false;
    }

    writeln!(
        report,
        "\n**TEST 2: {}**\n",
        if all_pass { "PASSED" } else { "FAILED" }
    )
    .unwrap();
}

// ============================================================================
// Test 3: Kinetic Step Verification
// ============================================================================

fn test_3(report: &mut String, p: &DiagParams) {
    writeln!(report, "## Test 3: Kinetic Step Verification\n").unwrap();

    let g = metric::euclidean::<3>();
    let n = p.n;
    let a = p.a_init;
    let dt = p.dt_at(a);

    writeln!(
        report,
        "| Mode | analytic theta | mean theta | std dev | |mean-analytic| | max |d|alpha|^2| |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|---|---|").unwrap();

    let mut all_pass = true;

    for (label, mode_n) in [("low (2)", 2_usize), ("mid (8)", 8), ("high (16)", 16)] {
        let k0 = 2.0 * PI * mode_n as f64 / p.box_length;
        let theta_analytic = -p.ell * k0 * k0 * dt / (2.0 * p.mass * a * a);

        let mut psi = EvenField::from_fn(&p.morphis_grid, g, |x| {
            let phase = k0 * x[0];
            (phase.cos(), phase.sin())
        });

        // Record before.
        let before: Vec<(f64, f64, f64)> = (0..n)
            .map(|m0| {
                let s = psi.scalar[IxDyn(&[m0, 0, 0])];
                let ps = psi.pseudoscalar[IxDyn(&[m0, 0, 0])];
                (ps.atan2(s), s * s + ps * ps, 0.0)
            })
            .collect();

        kinetic_step(&mut psi, &p.morphis_grid, p.ell, p.mass, a, dt);

        // Record after.
        let mut dphases = Vec::with_capacity(n);
        let mut max_mod_diff = 0.0_f64;

        for m0 in 0..n {
            let s = psi.scalar[IxDyn(&[m0, 0, 0])];
            let ps = psi.pseudoscalar[IxDyn(&[m0, 0, 0])];
            let phase_after = ps.atan2(s);
            let mod_sq_after = s * s + ps * ps;

            let mut dp = phase_after - before[m0].0;
            while dp > PI {
                dp -= 2.0 * PI;
            }
            while dp < -PI {
                dp += 2.0 * PI;
            }
            dphases.push(dp);

            let md = (mod_sq_after - before[m0].1).abs();
            if md > max_mod_diff {
                max_mod_diff = md;
            }
        }

        let mean = dphases.iter().sum::<f64>() / n as f64;
        let std = (dphases.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
        let err = (mean - theta_analytic).abs();

        if std > 1e-10 || err > 1e-10 || max_mod_diff > 1e-12 {
            all_pass = false;
        }

        writeln!(report, "| {label} | {theta_analytic:.6e} | {mean:.6e} | {std:.2e} | {err:.2e} | {max_mod_diff:.2e} |").unwrap();
    }

    writeln!(
        report,
        "\n**TEST 3: {}**\n",
        if all_pass { "PASSED" } else { "FAILED" }
    )
    .unwrap();
}

// ============================================================================
// Test 4: Madelung Round Trip
// ============================================================================

fn test_4(report: &mut String, p: &DiagParams) {
    writeln!(report, "## Test 4: Madelung Round Trip\n").unwrap();

    let g = metric::euclidean::<3>();
    let n = p.n;
    let k0 = 2.0 * PI * 3.0 / p.box_length;
    let k1 = 2.0 * PI * 5.0 / p.box_length;

    let psi = EvenField::from_fn(&p.morphis_grid, g, |x| {
        let rho = p.rho_bar * (1.0 + 0.1 * (k1 * x[0]).cos());
        let amp = (rho / p.mass).sqrt();
        let phase = k0 * x[0];
        (amp * phase.cos(), amp * phase.sin())
    });

    // Density round-trip.
    let rho_field = psi.density(p.mass);
    let mut max_rho_err = 0.0_f64;
    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * p.dx();
        let expected = p.rho_bar * (1.0 + 0.1 * (k1 * x).cos());
        let actual = rho_field.at(&[m0, 0, 0]).component(&[]);
        let err = (actual - expected).abs();
        if err > max_rho_err {
            max_rho_err = err;
        }
    }

    // Velocity round-trip via spectral derivative of phase.
    let v_expected = p.ell_over_m * k0;
    let n_c = n / 2 + 1;
    let handler = ndrustfft::R2cFftHandler::new(n);
    let mut phase_arr = ndarray::Array1::<f64>::zeros(n);
    for m0 in 0..n {
        let s = psi.scalar[IxDyn(&[m0, 0, 0])];
        let ps = psi.pseudoscalar[IxDyn(&[m0, 0, 0])];
        phase_arr[m0] = ps.atan2(s);
    }
    let mut phase_hat = ndarray::Array1::<Complex64>::zeros(n_c);
    ndrustfft::ndfft_r2c(&phase_arr, &mut phase_hat, &handler, 0);

    let mut deriv_hat = ndarray::Array1::<Complex64>::zeros(n_c);
    for m in 0..n_c {
        let k = wavenum(m, n, p.box_length);
        deriv_hat[m] = phase_hat[m] * Complex64::new(0.0, k);
    }
    let mut deriv = ndarray::Array1::<f64>::zeros(n);
    ndrustfft::ndifft_r2c(&deriv_hat, &mut deriv, &handler, 0);

    let mut max_v_err = 0.0_f64;
    for m0 in 0..n {
        let v_actual = p.ell_over_m * deriv[m0];
        let err = (v_actual - v_expected).abs();
        if err > max_v_err {
            max_v_err = err;
        }
    }

    let rel_rho = max_rho_err / p.rho_bar;
    let rel_v = max_v_err / v_expected.abs();

    writeln!(
        report,
        "- Max |rho_out - rho_in|: {max_rho_err:.4e} (relative: {rel_rho:.4e})"
    )
    .unwrap();
    writeln!(
        report,
        "- Max |v_out - v_in|: {max_v_err:.4e} (relative: {rel_v:.4e})"
    )
    .unwrap();

    let pass = rel_rho < 1e-10 && rel_v < 1e-4;
    writeln!(
        report,
        "\n**TEST 4: {}**\n",
        if pass { "PASSED" } else { "FAILED" }
    )
    .unwrap();
}

// ============================================================================
// Test 5: Energy Conservation Without Gravity
// ============================================================================

fn test_5(report: &mut String, p: &DiagParams) {
    writeln!(report, "## Test 5: Energy Conservation Without Gravity\n").unwrap();

    let g = metric::euclidean::<3>();
    let n = p.n;
    let n_c = n / 2 + 1;
    let a = p.a_init;
    let dt = p.dt_at(a);

    let sigma = 5.0 * p.dx();
    let center = p.box_length / 2.0;
    let mut psi = EvenField::from_fn(&p.morphis_grid, g, |x| {
        let r2 = (x[0] - center).powi(2) + (x[1] - center).powi(2) + (x[2] - center).powi(2);
        ((-r2 / (2.0 * sigma * sigma)).exp(), 0.0)
    });

    let compute_ekin = |psi: &EvenField<3>| -> f64 {
        let s_hat = fft_3d(&psi.scalar, n);
        let p_hat = fft_3d(&psi.pseudoscalar, n);
        let mut e = 0.0;
        for m0 in 0..n {
            let kx = wavenum(m0, n, p.box_length);
            for m1 in 0..n {
                let ky = wavenum(m1, n, p.box_length);
                for m2 in 0..n_c {
                    let kz = wavenum(m2, n, p.box_length);
                    let k2 = kx * kx + ky * ky + kz * kz;
                    let pw = s_hat[[m0, m1, m2]].norm_sqr() + p_hat[[m0, m1, m2]].norm_sqr();
                    let w = if m2 == 0 || m2 == n / 2 { 1.0 } else { 2.0 };
                    e += w * k2 * pw;
                }
            }
        }
        let n3 = (n * n * n) as f64;
        e * p.ell * p.ell / (2.0 * p.mass * a * a) / (n3 * n3) * p.dx().powi(3)
    };

    let e0 = compute_ekin(&psi);
    for _ in 0..100 {
        kinetic_step(&mut psi, &p.morphis_grid, p.ell, p.mass, a, dt);
    }
    let e_final = compute_ekin(&psi);

    let drift = (e_final - e0).abs() / e0.abs().max(1e-30);

    // Mass conservation.
    let rho = psi.density(p.mass);
    let mut total_mass = 0.0;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                total_mass += rho.at(&[m0, m1, m2]).component(&[]);
            }
        }
    }
    total_mass *= p.dx().powi(3);

    let psi0 = EvenField::from_fn(&p.morphis_grid, g, |x| {
        let r2 = (x[0] - center).powi(2) + (x[1] - center).powi(2) + (x[2] - center).powi(2);
        ((-r2 / (2.0 * sigma * sigma)).exp(), 0.0)
    });
    let rho0 = psi0.density(p.mass);
    let mut m0_total = 0.0;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                m0_total += rho0.at(&[m0, m1, m2]).component(&[]);
            }
        }
    }
    m0_total *= p.dx().powi(3);

    let mass_drift = (total_mass - m0_total).abs() / m0_total.abs().max(1e-30);

    writeln!(report, "- E_kin(initial): {e0:.10e}").unwrap();
    writeln!(report, "- E_kin(final, 100 steps): {e_final:.10e}").unwrap();
    writeln!(report, "- E_kin relative drift: {drift:.4e}").unwrap();
    writeln!(report, "- Mass(initial): {m0_total:.10e}").unwrap();
    writeln!(report, "- Mass(final): {total_mass:.10e}").unwrap();
    writeln!(report, "- Mass relative drift: {mass_drift:.4e}").unwrap();

    let pass = drift < 1e-10 && mass_drift < 1e-10;
    writeln!(
        report,
        "\n**TEST 5: {}**\n",
        if pass { "PASSED" } else { "FAILED" }
    )
    .unwrap();
}

// ============================================================================
// Test 6: Static Background Test
// ============================================================================

fn test_6(report: &mut String, p: &DiagParams) {
    writeln!(report, "## Test 6: Static Background Test\n").unwrap();

    let g = metric::euclidean::<3>();
    let n = p.n;
    let uniform = (p.rho_bar / p.mass).sqrt();

    let params = FieldParams {
        smoothing_length: p.ell,
        mass_alpha: p.mass,
    };
    let psi = EvenField::from_fn(&p.morphis_grid, g, |_| (uniform, 0.0));
    let field_state = FieldState {
        grid: p.morphis_grid.clone(),
        alpha: Some(psi),
        beta: None,
        gamma: None,
        params,
    };

    let mut content = Content::Fields(field_state);
    let hermes_grid = HermesGrid::new(p.n, p.box_length);
    let mut dynamics = SchrodingerPoissonDynamics::new(PoissonGravity::new(hermes_grid));
    let da = (p.a_final - p.a_init) / 2000.0;
    let target = uniform * uniform;
    let mut max_dev = 0.0_f64;

    for step in 0..100 {
        let a_prev = p.a_init + step as f64 * da;
        dynamics
            .step(&mut content, &p.cosmology, a_prev, a_prev + da)
            .unwrap();

        let psi_ref = content.fields().unwrap().alpha.as_ref().unwrap();
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    let s = psi_ref.scalar[IxDyn(&[m0, m1, m2])];
                    let ps = psi_ref.pseudoscalar[IxDyn(&[m0, m1, m2])];
                    let dev = (s * s + ps * ps - target).abs();
                    if dev > max_dev {
                        max_dev = dev;
                    }
                }
            }
        }
    }

    let rel = max_dev / target;
    writeln!(report, "- |alpha|^2 target: {target:.10e}").unwrap();
    writeln!(
        report,
        "- Max deviation: {max_dev:.4e} (relative: {rel:.4e})"
    )
    .unwrap();

    let pass = rel < 1e-10;
    writeln!(
        report,
        "\n**TEST 6: {}**\n",
        if pass { "PASSED" } else { "FAILED" }
    )
    .unwrap();
}

// ============================================================================
// Test 7: Linear Growth Test
// ============================================================================

fn test_7(report: &mut String, p: &DiagParams) {
    writeln!(report, "## Test 7: Linear Growth Test\n").unwrap();

    let g = metric::euclidean::<3>();
    let n = p.n;
    let n_c = n / 2 + 1;
    let n_steps = 200;
    let amp_init = 0.01;
    let da = (p.a_final - p.a_init) / 2000.0;
    let growth_init = p.cosmology.growth_factor(p.a_init);
    let k_j = (16.0 * PI * GRAV * p.rho_bar / (p.ell_over_m * p.ell_over_m)).powf(0.25);

    writeln!(
        report,
        "| Mode | k/k_J | A_sim(final) | A_lin(final) | rel error |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|---|").unwrap();

    for mode_n in [1_usize, 2, 4, 8, 16] {
        let k0 = 2.0 * PI * mode_n as f64 / p.box_length;

        let params = FieldParams {
            smoothing_length: p.ell,
            mass_alpha: p.mass,
        };
        let psi = EvenField::from_fn(&p.morphis_grid, g, |x| {
            let rho = p.rho_bar * (1.0 + amp_init * (k0 * x[0]).cos());
            ((rho.max(1e-30) / p.mass).sqrt(), 0.0)
        });
        let field_state = FieldState {
            grid: p.morphis_grid.clone(),
            alpha: Some(psi),
            beta: None,
            gamma: None,
            params,
        };
        let mut content = Content::Fields(field_state);
        let hermes_grid = HermesGrid::new(p.n, p.box_length);
        let mut dynamics = SchrodingerPoissonDynamics::new(PoissonGravity::new(hermes_grid));

        let mut a_current = p.a_init;
        for step in 0..n_steps {
            let a_prev = p.a_init + step as f64 * da;
            dynamics
                .step(&mut content, &p.cosmology, a_prev, a_prev + da)
                .unwrap();
            a_current = a_prev + da;
        }

        // Extract Fourier amplitude at k0.
        let psi_f = content.fields().unwrap().alpha.as_ref().unwrap();
        let rho_f = psi_f.density(p.mass);
        let mut delta_grid = Array3::<f64>::zeros((n, n, n));
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    delta_grid[[m0, m1, m2]] =
                        rho_f.at(&[m0, m1, m2]).component(&[]) / p.rho_bar - 1.0;
                }
            }
        }

        let delta_hat = {
            let h_r = ndrustfft::R2cFftHandler::new(n);
            let h1 = ndrustfft::FftHandler::new(n);
            let h0 = ndrustfft::FftHandler::new(n);
            let mut c = Array3::<Complex64>::zeros((n, n, n_c));
            ndrustfft::ndfft_r2c(&delta_grid, &mut c, &h_r, 2);
            let mut s = c.clone();
            ndrustfft::ndfft(&c, &mut s, &h1, 1);
            c.assign(&s);
            ndrustfft::ndfft(&c, &mut s, &h0, 0);
            c.assign(&s);
            c
        };

        let n3 = (n * n * n) as f64;
        let amp_sim = delta_hat[[mode_n, 0, 0]].norm() / n3 * 2.0;
        let amp_lin = amp_init * p.cosmology.growth_factor(a_current) / growth_init;
        let rel = (amp_sim - amp_lin).abs() / amp_lin.abs().max(1e-30);

        writeln!(
            report,
            "| k={mode_n} k_f | {:.3} | {amp_sim:.6e} | {amp_lin:.6e} | {rel:.4e} |",
            k0 / k_j
        )
        .unwrap();
    }

    writeln!(report, "\nModes with k/k_J < 1 should track linear theory.").unwrap();
    writeln!(
        report,
        "Modes with k/k_J > 1 should show suppressed growth.\n"
    )
    .unwrap();
    writeln!(report, "**TEST 7: see table for interpretation**\n").unwrap();
}

// ============================================================================
// Test 8: Cosmic Web Run
// ============================================================================

fn test_8(report: &mut String, p: &DiagParams) {
    writeln!(
        report,
        "## Test 8: Cosmic Web Run (ell/m = {:.0})\n",
        p.ell_over_m
    )
    .unwrap();

    let g = metric::euclidean::<3>();
    let n = p.n;
    let n_c = n / 2 + 1;
    let n_steps = 200;
    let da = (p.a_final - p.a_init) / 2000.0;
    let delta_rms_init = 0.1;

    use rand::Rng;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    let k_f = p.k_fundamental();
    let k_min = 1.5 * k_f;
    let k_max = 0.5 * p.k_nyquist();

    let mut rng = ChaCha20Rng::seed_from_u64(42);
    let mut delta_hat = Array3::<Complex64>::zeros((n, n, n_c));
    for m0 in 0..n {
        let kx = wavenum(m0, n, p.box_length);
        for m1 in 0..n {
            let ky = wavenum(m1, n, p.box_length);
            for m2 in 0..n_c {
                let kz = wavenum(m2, n, p.box_length);
                let k = (kx * kx + ky * ky + kz * kz).sqrt();
                if k < k_min || k > k_max {
                    continue;
                }
                let a_k = k_f / k;
                let u1: f64 = rng.random_range(1e-10..1.0);
                let u2: f64 = rng.random_range(0.0..2.0 * PI);
                let gr = (-2.0 * u1.ln()).sqrt() * u2.cos();
                let gi = (-2.0 * u1.ln()).sqrt() * u2.sin();
                delta_hat[[m0, m1, m2]] = Complex64::new(gr * a_k, gi * a_k);
            }
        }
    }

    let delta_raw = ifft_3d(&delta_hat, n);
    let rms = (delta_raw.iter().map(|d| d * d).sum::<f64>() / delta_raw.len() as f64).sqrt();
    let norm = delta_rms_init / rms.max(1e-30);

    let params = FieldParams {
        smoothing_length: p.ell,
        mass_alpha: p.mass,
    };
    let psi = EvenField::from_fn(&p.morphis_grid, g, |x| {
        let m0 = ((x[0] / p.dx()) as usize).min(n - 1);
        let m1 = ((x[1] / p.dx()) as usize).min(n - 1);
        let m2 = ((x[2] / p.dx()) as usize).min(n - 1);
        let delta = norm * delta_raw[[m0, m1, m2]];
        let rho = p.rho_bar * (1.0_f64 + delta).max(1e-10);
        ((rho / p.mass).sqrt(), 0.0)
    });

    let field_state = FieldState {
        grid: p.morphis_grid.clone(),
        alpha: Some(psi),
        beta: None,
        gamma: None,
        params,
    };
    let mut content = Content::Fields(field_state);
    let hermes_grid = HermesGrid::new(p.n, p.box_length);
    let mut dynamics = SchrodingerPoissonDynamics::new(PoissonGravity::new(hermes_grid));

    writeln!(
        report,
        "| step | a | delta_std | delta_max | E_kin | E_pot | E_kin/E_pot |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|---|---|---|").unwrap();

    for step in 0..=n_steps {
        if step > 0 {
            let a_prev = p.a_init + (step - 1) as f64 * da;
            dynamics
                .step(&mut content, &p.cosmology, a_prev, a_prev + da)
                .unwrap();
        }
        if step % 20 != 0 {
            continue;
        }

        let a_now = p.a_init + step as f64 * da;
        let psi_r = content.fields().unwrap().alpha.as_ref().unwrap();
        let rho_f = psi_r.density(p.mass);

        let mut deltas: Vec<f64> = Vec::with_capacity(n * n * n);
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    deltas.push(rho_f.at(&[m0, m1, m2]).component(&[]) / p.rho_bar - 1.0);
                }
            }
        }

        let d_mean = deltas.iter().sum::<f64>() / deltas.len() as f64;
        let d_std =
            (deltas.iter().map(|d| (d - d_mean).powi(2)).sum::<f64>() / deltas.len() as f64).sqrt();
        let d_max = deltas.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        // Kinetic energy.
        let s_hat = fft_3d(&psi_r.scalar, n);
        let p_hat = fft_3d(&psi_r.pseudoscalar, n);
        let a2 = a_now * a_now;
        let n3 = (n * n * n) as f64;
        let cv = p.dx().powi(3);
        let mut e_kin = 0.0;
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n_c {
                    let k2 = wavenum(m0, n, p.box_length).powi(2)
                        + wavenum(m1, n, p.box_length).powi(2)
                        + wavenum(m2, n, p.box_length).powi(2);
                    let pw = s_hat[[m0, m1, m2]].norm_sqr() + p_hat[[m0, m1, m2]].norm_sqr();
                    let w = if m2 == 0 || m2 == n / 2 { 1.0 } else { 2.0 };
                    e_kin += w * k2 * pw;
                }
            }
        }
        e_kin *= p.ell * p.ell / (2.0 * p.mass * a2) / (n3 * n3) * cv;

        // Potential energy.
        let delta_field = Field::scalar_field(&p.morphis_grid, g, |x| {
            let m0 = ((x[0] / p.dx()) as usize).min(n - 1);
            let m1 = ((x[1] / p.dx()) as usize).min(n - 1);
            let m2 = ((x[2] / p.dx()) as usize).min(n - 1);
            deltas[m0 * n * n + m1 * n + m2]
        });
        let phi = (&delta_field * (4.0 * PI * GRAV * p.rho_bar * a2)).laplacian_inverse();
        let mut e_pot = 0.0;
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    e_pot += rho_f.at(&[m0, m1, m2]).component(&[])
                        * phi.at(&[m0, m1, m2]).component(&[]);
                }
            }
        }
        e_pot *= 0.5 * cv;

        let ratio = if e_pot.abs() > 1e-30 {
            e_kin / e_pot
        } else {
            f64::NAN
        };

        writeln!(report, "| {step:4} | {a_now:.4} | {d_std:.4e} | {d_max:.4e} | {e_kin:.4e} | {e_pot:.4e} | {ratio:.4} |").unwrap();
    }

    writeln!(report, "\n**TEST 8: see table for interpretation**\n").unwrap();
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    let p = DiagParams::new();
    let mut report = String::with_capacity(16384);

    writeln!(report, "# Schrodinger-Poisson Diagnostic Report\n").unwrap();
    writeln!(
        report,
        "Generated: {}\n",
        chrono::Local::now().format("%Y-%m-%d %H:%M:%S")
    )
    .unwrap();

    for (num, name, func) in [
        (
            1,
            "Dimensional Sanity Check",
            test_1 as fn(&mut String, &DiagParams),
        ),
        (2, "Poisson Solver Verification", test_2),
        (3, "Kinetic Step Verification", test_3),
        (4, "Madelung Round Trip", test_4),
        (5, "Energy Conservation Without Gravity", test_5),
        (6, "Static Background Test", test_6),
        (7, "Linear Growth Test", test_7),
        (8, "Cosmic Web Run", test_8),
    ] {
        println!("Running Test {num}: {name}...");
        func(&mut report, &p);
    }

    writeln!(report, "---\n## Summary\n").unwrap();
    writeln!(
        report,
        "See individual test sections above for PASSED/FAILED status."
    )
    .unwrap();

    std::fs::create_dir_all("data").ok();
    std::fs::write("data/diagnostics-report.md", &report).unwrap();
    println!("\nReport written to data/diagnostics-report.md");
}
