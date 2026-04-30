//! Second-pass diagnostics: Poisson normalization and Madelung velocity.
//!
//! Run with:
//!   cargo run --example diagnostics_pass2 --release

use std::f64::consts::PI;
use std::fmt::Write;

use morphis::even_field::EvenField;
use morphis::field::Field;
use morphis::grid::Grid as MorphisGrid;
use morphis::metric;
use ndarray::{Array3, IxDyn};
use num_complex::Complex64;

// hermes fft_3d/ifft_3d not used directly; local static versions below.
use hermes_rs::physics::constants::G as GRAV;
use hermes_rs::physics::cosmology::planck_2018;

// ============================================================================
// Shared helpers
// ============================================================================

fn wavenum(m: usize, n: usize, box_length: f64) -> f64 {
    let freq = if m <= n / 2 {
        m as f64
    } else {
        m as f64 - n as f64
    };
    2.0 * PI * freq / box_length
}

// ============================================================================
// Section A: Poisson Solver Single-Mode Probe
// ============================================================================

fn section_a(report: &mut String) {
    writeln!(report, "## Section A: Poisson Solver Single-Mode Probe\n").unwrap();

    let n = 64_usize;
    let box_length = 10000.0;
    let cosmology = planck_2018();
    let rho_bar = cosmology.density_matter();
    let a = 0.33;
    let prefactor = 4.0 * PI * GRAV * rho_bar * a * a;
    let amp = 0.1;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();

    writeln!(report, "Prefactor 4 pi G rho_bar a^2 = {prefactor:.6e}\n").unwrap();

    let modes: [(&str, [usize; 3]); 3] = [
        ("(1,0,0)", [1, 0, 0]),
        ("(2,0,0)", [2, 0, 0]),
        ("(2,1,1)", [2, 1, 1]),
    ];

    writeln!(
        report,
        "| Mode | k^2 | delta_hat | Phi_hat actual | Phi_hat analytic | ratio |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|---|---|").unwrap();

    for (label, mode_idx) in &modes {
        let kx = wavenum(mode_idx[0], n, box_length);
        let ky = wavenum(mode_idx[1], n, box_length);
        let kz = wavenum(mode_idx[2], n, box_length);
        let k2 = kx * kx + ky * ky + kz * kz;

        // Construct delta = A cos(k0 . x).
        let delta = Field::scalar_field(&grid, g, |x| {
            amp * (kx * x[0] + ky * x[1] + kz * x[2]).cos()
        });

        // FFT delta to check its Fourier amplitude.
        // We need the raw data from the Field. Extract via .at().
        let mut delta_grid = Array3::<f64>::zeros((n, n, n));
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    delta_grid[[m0, m1, m2]] = delta.at(&[m0, m1, m2]).component(&[]);
                }
            }
        }

        // FFT using the same routine as the Poisson solver uses internally.
        let delta_hat = fft_3d_static(&delta_grid, n);
        let delta_coeff = delta_hat[[mode_idx[0], mode_idx[1], mode_idx[2]]];

        // Run the actual Poisson solver.
        let source = &delta * prefactor;
        let phi = source.laplacian_inverse();

        // Extract phi on the grid.
        let mut phi_grid = Array3::<f64>::zeros((n, n, n));
        for m0 in 0..n {
            for m1 in 0..n {
                for m2 in 0..n {
                    phi_grid[[m0, m1, m2]] = phi.at(&[m0, m1, m2]).component(&[]);
                }
            }
        }

        // FFT phi to get its Fourier amplitude at the mode.
        let phi_hat = fft_3d_static(&phi_grid, n);
        let phi_coeff = phi_hat[[mode_idx[0], mode_idx[1], mode_idx[2]]];

        // Analytic: Phi_hat = -source_hat / k^2 = -prefactor * delta_hat / k^2
        let phi_analytic = -prefactor * delta_coeff / k2;

        let ratio = phi_coeff.norm() / phi_analytic.norm();

        writeln!(
            report,
            "| {label} | {k2:.6e} | {:.6e} | {:.6e} | {:.6e} | {ratio:.8} |",
            delta_coeff.norm(),
            phi_coeff.norm(),
            phi_analytic.norm()
        )
        .unwrap();

        // Also report the raw complex values for debugging.
        writeln!(
            report,
            "|   | | delta_hat = ({:.6e}, {:.6e}) | phi_hat = ({:.6e}, {:.6e}) | phi_analytic = ({:.6e}, {:.6e}) | |",
            delta_coeff.re, delta_coeff.im,
            phi_coeff.re, phi_coeff.im,
            phi_analytic.re, phi_analytic.im
        )
        .unwrap();
    }

    writeln!(report).unwrap();
}

// ============================================================================
// Section B: FFT Convention Inspection
// ============================================================================

fn section_b(report: &mut String) {
    writeln!(report, "## Section B: FFT Convention Inspection\n").unwrap();

    let n = 64_usize;
    let n_c = n / 2 + 1;
    let box_length = 10000.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let n3 = (n * n * n) as f64;

    // Test 1: constant field f(x) = 1.
    writeln!(report, "### Constant Field f(x) = 1\n").unwrap();

    let f_const = Field::scalar_field(&grid, g, |_| 1.0);
    let mut f_grid = Array3::<f64>::zeros((n, n, n));
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                f_grid[[m0, m1, m2]] = f_const.at(&[m0, m1, m2]).component(&[]);
            }
        }
    }

    let f_hat = fft_3d_static(&f_grid, n);
    let zero_mode = f_hat[[0, 0, 0]];

    let mut max_nonzero = 0.0_f64;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n_c {
                if m0 == 0 && m1 == 0 && m2 == 0 {
                    continue;
                }
                let mag = f_hat[[m0, m1, m2]].norm();
                if mag > max_nonzero {
                    max_nonzero = mag;
                }
            }
        }
    }

    writeln!(
        report,
        "- Zero-mode coefficient: ({:.6e}, {:.6e}), magnitude = {:.6e}",
        zero_mode.re,
        zero_mode.im,
        zero_mode.norm()
    )
    .unwrap();
    writeln!(report, "- N^3 = {n3:.0}").unwrap();
    writeln!(report, "- N = {n}").unwrap();
    writeln!(report, "- Zero-mode / N^3 = {:.6e}", zero_mode.re / n3).unwrap();
    writeln!(report, "- Zero-mode / N = {:.6e}", zero_mode.re / n as f64).unwrap();
    writeln!(report, "- Max non-zero mode magnitude: {max_nonzero:.4e}\n").unwrap();

    if (zero_mode.re - n3).abs() < 1.0 {
        writeln!(
            report,
            "**Convention: forward FFT is un-normalized (zero mode = N^3)**\n"
        )
        .unwrap();
    } else if (zero_mode.re - n as f64).abs() < 1.0 {
        writeln!(
            report,
            "**Convention: forward FFT normalizes by N^2 (zero mode = N)**\n"
        )
        .unwrap();
    } else if (zero_mode.re - 1.0).abs() < 1e-10 {
        writeln!(
            report,
            "**Convention: forward FFT is fully normalized (zero mode = 1)**\n"
        )
        .unwrap();
    } else {
        writeln!(
            report,
            "**Convention: UNKNOWN (zero mode = {:.6e})**\n",
            zero_mode.re
        )
        .unwrap();
    }

    // Test 2: cosine field f(x) = cos(k0 . x) with k0 = (2,0,0) k_f.
    writeln!(report, "### Cosine Field f(x) = cos(2 k_f x)\n").unwrap();

    let k0 = 2.0 * 2.0 * PI / box_length;
    let f_cos = Field::scalar_field(&grid, g, |x| (k0 * x[0]).cos());
    let mut fc_grid = Array3::<f64>::zeros((n, n, n));
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                fc_grid[[m0, m1, m2]] = f_cos.at(&[m0, m1, m2]).component(&[]);
            }
        }
    }

    let fc_hat = fft_3d_static(&fc_grid, n);

    let coeff_plus = fc_hat[[2, 0, 0]];
    // Negative frequency for R2C: k0 = (N-2, 0, 0) maps to (-2, 0, 0).
    // But for R2C, only m2 in [0, N/2] is stored. The (N-2, 0, 0) mode IS stored.
    let coeff_minus = fc_hat[[n - 2, 0, 0]];
    let coeff_zero = fc_hat[[0, 0, 0]];

    writeln!(
        report,
        "- Coeff at (2,0,0): ({:.6e}, {:.6e}), magnitude = {:.6e}",
        coeff_plus.re,
        coeff_plus.im,
        coeff_plus.norm()
    )
    .unwrap();
    writeln!(
        report,
        "- Coeff at (N-2,0,0): ({:.6e}, {:.6e}), magnitude = {:.6e}",
        coeff_minus.re,
        coeff_minus.im,
        coeff_minus.norm()
    )
    .unwrap();
    writeln!(
        report,
        "- Zero-mode: ({:.6e}, {:.6e})",
        coeff_zero.re, coeff_zero.im
    )
    .unwrap();
    writeln!(
        report,
        "- Coeff / (N^3 / 2) = {:.6e}",
        coeff_plus.norm() / (n3 / 2.0)
    )
    .unwrap();
    writeln!(
        report,
        "- Coeff / (N / 2) = {:.6e}",
        coeff_plus.norm() / (n as f64 / 2.0)
    )
    .unwrap();
    writeln!(
        report,
        "- Coeff / (1/2) = {:.6e}\n",
        coeff_plus.norm() / 0.5
    )
    .unwrap();

    // Also test morphis's own FFT via laplacian_inverse's internal path.
    // We can't directly access morphis's FFT, but we can infer its convention
    // from the laplacian_inverse result.
    writeln!(report, "### Morphis laplacian_inverse convention probe\n").unwrap();
    writeln!(
        report,
        "Testing: source = cos(k0 x), laplacian_inverse should give -cos(k0 x) / k0^2\n"
    )
    .unwrap();

    let source = Field::scalar_field(&grid, g, |x| (k0 * x[0]).cos());
    let result = source.laplacian_inverse();

    // Sample at a few points.
    let k2 = k0 * k0;
    writeln!(
        report,
        "| x_index | source | result | analytic (-1/k^2 * source) | ratio |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|---|").unwrap();
    for m0 in [0, 8, 16, 32] {
        let s = source.at(&[m0, 0, 0]).component(&[]);
        let r = result.at(&[m0, 0, 0]).component(&[]);
        let analytic = -s / k2;
        let ratio = if analytic.abs() > 1e-30 {
            r / analytic
        } else {
            f64::NAN
        };
        writeln!(
            report,
            "| {m0} | {s:.6e} | {r:.6e} | {analytic:.6e} | {ratio:.8} |"
        )
        .unwrap();
    }

    writeln!(report).unwrap();
}

// ============================================================================
// Section C: Round-Trip FFT
// ============================================================================

fn section_c(report: &mut String) {
    writeln!(report, "## Section C: Round-Trip FFT\n").unwrap();

    let n = 64_usize;

    // Random field with fixed seed.
    use rand::Rng;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    let mut rng = ChaCha20Rng::seed_from_u64(12345);

    let mut original = Array3::<f64>::zeros((n, n, n));
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                original[[m0, m1, m2]] = rng.random_range(-1.0..1.0);
            }
        }
    }

    // Forward FFT (using hermes's fft_3d which wraps ndrustfft).
    let forward = fft_3d_static(&original, n);

    // Inverse FFT.
    let roundtrip = ifft_3d_static(&forward, n);

    let mut max_diff = 0.0_f64;
    for m0 in 0..n {
        for m1 in 0..n {
            for m2 in 0..n {
                let diff = (roundtrip[[m0, m1, m2]] - original[[m0, m1, m2]]).abs();
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }
    }

    writeln!(report, "- Max |roundtrip - original|: {max_diff:.4e}").unwrap();
    let pass = max_diff < 1e-10;
    writeln!(
        report,
        "\n**Section C: {}**\n",
        if pass { "PASSED" } else { "FAILED" }
    )
    .unwrap();
}

// ============================================================================
// Section D: Madelung Velocity Extraction Comparison
// ============================================================================

fn section_d(report: &mut String) {
    writeln!(report, "## Section D: Madelung Velocity Extraction\n").unwrap();

    let n = 64_usize;
    let n_c = n / 2 + 1;
    let box_length = 10000.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let dx = box_length / n as f64;

    let cosmology = planck_2018();
    let rho_bar = cosmology.density_matter();
    let ell_over_m = 7500.0;
    let mass = 1e10;
    let ell = ell_over_m * mass;

    let k0 = 2.0 * PI * 3.0 / box_length;
    let v0 = 100.0; // kpc/Gyr

    // Phase: S(x) = -(m v0) / (ell k0) * cos(k0 x)
    // Velocity: v_x(x) = (ell/m) dS/dx = v0 sin(k0 x)
    let phase_amp = mass * v0 / (ell * k0);

    writeln!(report, "- v0 = {v0} kpc/Gyr").unwrap();
    writeln!(report, "- k0 = {k0:.6e} kpc^-1").unwrap();
    writeln!(
        report,
        "- phase amplitude m*v0/(ell*k0) = {phase_amp:.6e} rad\n"
    )
    .unwrap();

    // Construct alpha.
    let psi = EvenField::from_fn(&grid, g, |x| {
        let amp = (rho_bar / mass).sqrt();
        let phase = -phase_amp * (k0 * x[0]).cos();
        (amp * phase.cos(), amp * phase.sin())
    });

    // Method 1: Phase-gradient extraction.
    // v = (ell/m) * d(arg(alpha))/dx via spectral derivative.
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
        let k = wavenum(m, n, box_length);
        deriv_hat[m] = phase_hat[m] * Complex64::new(0.0, k);
    }
    let mut deriv_phase = ndarray::Array1::<f64>::zeros(n);
    ndrustfft::ndifft_r2c(&deriv_hat, &mut deriv_phase, &handler, 0);

    // Method 2: Wavefunction-gradient extraction.
    // v = (ell / (m |alpha|^2)) * Im(alpha_bar * d(alpha)/dx)
    // alpha = scalar + pseudoscalar * I
    // alpha_bar = scalar - pseudoscalar * I
    // d(alpha)/dx = d(scalar)/dx + d(pseudoscalar)/dx * I
    // alpha_bar * d(alpha)/dx = (scalar * d(scalar)/dx + pseudo * d(pseudo)/dx)
    //                         + (scalar * d(pseudo)/dx - pseudo * d(scalar)/dx) * I
    // Im part = scalar * d(pseudo)/dx - pseudo * d(scalar)/dx

    // Spectral derivatives of scalar and pseudoscalar along x.
    let mut s_arr = ndarray::Array1::<f64>::zeros(n);
    let mut p_arr = ndarray::Array1::<f64>::zeros(n);
    for m0 in 0..n {
        s_arr[m0] = psi.scalar[IxDyn(&[m0, 0, 0])];
        p_arr[m0] = psi.pseudoscalar[IxDyn(&[m0, 0, 0])];
    }

    let mut s_hat = ndarray::Array1::<Complex64>::zeros(n_c);
    let mut p_hat = ndarray::Array1::<Complex64>::zeros(n_c);
    ndrustfft::ndfft_r2c(&s_arr, &mut s_hat, &handler, 0);
    ndrustfft::ndfft_r2c(&p_arr, &mut p_hat, &handler, 0);

    let mut ds_hat = ndarray::Array1::<Complex64>::zeros(n_c);
    let mut dp_hat = ndarray::Array1::<Complex64>::zeros(n_c);
    for m in 0..n_c {
        let k = wavenum(m, n, box_length);
        ds_hat[m] = s_hat[m] * Complex64::new(0.0, k);
        dp_hat[m] = p_hat[m] * Complex64::new(0.0, k);
    }

    let mut ds_real = ndarray::Array1::<f64>::zeros(n);
    let mut dp_real = ndarray::Array1::<f64>::zeros(n);
    ndrustfft::ndifft_r2c(&ds_hat, &mut ds_real, &handler, 0);
    ndrustfft::ndifft_r2c(&dp_hat, &mut dp_real, &handler, 0);

    // Compare both methods against analytic.
    writeln!(
        report,
        "| x_index | v_analytic | v_phase | v_wf | err_phase | err_wf |"
    )
    .unwrap();
    writeln!(report, "|---|---|---|---|---|---|").unwrap();

    let mut max_err_phase = 0.0_f64;
    let mut max_err_wf = 0.0_f64;

    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * dx;
        let v_analytic = v0 * (k0 * x).sin();

        let v_phase = ell_over_m * deriv_phase[m0];

        let s = s_arr[m0];
        let ps = p_arr[m0];
        let mod_sq = s * s + ps * ps;
        let im_part = s * dp_real[m0] - ps * ds_real[m0];
        let v_wf = ell_over_m * im_part / mod_sq;

        let err_p = (v_phase - v_analytic).abs();
        let err_w = (v_wf - v_analytic).abs();
        if err_p > max_err_phase {
            max_err_phase = err_p;
        }
        if err_w > max_err_wf {
            max_err_wf = err_w;
        }

        // Only print every 8th point.
        if m0 % 8 == 0 {
            writeln!(report, "| {m0} | {v_analytic:.4e} | {v_phase:.4e} | {v_wf:.4e} | {err_p:.4e} | {err_w:.4e} |").unwrap();
        }
    }

    writeln!(
        report,
        "\n- Max error (phase-gradient): {max_err_phase:.4e}"
    )
    .unwrap();
    writeln!(
        report,
        "- Max error (wavefunction-gradient): {max_err_wf:.4e}"
    )
    .unwrap();
    writeln!(
        report,
        "- Relative to v0={v0}: phase {:.4e}, wf {:.4e}\n",
        max_err_phase / v0,
        max_err_wf / v0
    )
    .unwrap();
}

// ============================================================================
// Section E: Density Round Trip
// ============================================================================

fn section_e(report: &mut String) {
    writeln!(report, "## Section E: Density Round Trip (Zero Velocity)\n").unwrap();

    let n = 64_usize;
    let box_length = 10000.0;
    let grid = MorphisGrid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();
    let dx = box_length / n as f64;

    let cosmology = planck_2018();
    let rho_bar = cosmology.density_matter();
    let mass = 1e10;

    let k1 = 2.0 * PI * 5.0 / box_length;

    // Zero velocity: alpha = sqrt(rho / m), purely real.
    let psi = EvenField::from_fn(&grid, g, |x| {
        let rho = rho_bar * (1.0 + 0.1 * (k1 * x[0]).cos());
        ((rho / mass).sqrt(), 0.0)
    });

    let rho_out = psi.density(mass);

    let mut max_abs = 0.0_f64;
    let mut max_rel = 0.0_f64;

    for m0 in 0..n {
        let x = (m0 as f64 + 0.5) * dx;
        let rho_expected = rho_bar * (1.0 + 0.1 * (k1 * x).cos());
        let rho_actual = rho_out.at(&[m0, 0, 0]).component(&[]);
        let abs_err = (rho_actual - rho_expected).abs();
        let rel_err = abs_err / rho_expected;
        if abs_err > max_abs {
            max_abs = abs_err;
        }
        if rel_err > max_rel {
            max_rel = rel_err;
        }
    }

    writeln!(report, "- Max absolute error: {max_abs:.4e}").unwrap();
    writeln!(report, "- Max relative error: {max_rel:.4e}").unwrap();

    let pass = max_rel < 1e-10;
    writeln!(
        report,
        "\n**Section E: {}**{}\n",
        if pass { "PASSED" } else { "FAILED" },
        if !pass {
            format!(" — relative error {max_rel:.4e}")
        } else {
            String::new()
        }
    )
    .unwrap();
}

// ============================================================================
// FFT helpers (operating on Array3, not ArrayD)
// ============================================================================

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

fn ifft_3d_static(complex: &Array3<Complex64>, n: usize) -> Array3<f64> {
    let handler_c2c_0 = ndrustfft::FftHandler::new(n);
    let handler_c2c_1 = ndrustfft::FftHandler::new(n);
    let handler_r2c = ndrustfft::R2cFftHandler::new(n);

    let mut work = complex.clone();
    let mut scratch = work.clone();
    ndrustfft::ndifft(&work, &mut scratch, &handler_c2c_0, 0);
    work.assign(&scratch);
    ndrustfft::ndifft(&work, &mut scratch, &handler_c2c_1, 1);
    work.assign(&scratch);

    let mut real = Array3::<f64>::zeros((n, n, n));
    ndrustfft::ndifft_r2c(&work, &mut real, &handler_r2c, 2);
    real
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    let mut report = String::with_capacity(8192);

    writeln!(report, "# Second-Pass Diagnostic Report\n").unwrap();
    writeln!(
        report,
        "Generated: {}\n",
        chrono::Local::now().format("%Y-%m-%d %H:%M:%S")
    )
    .unwrap();

    println!("Section A: Poisson Solver Single-Mode Probe...");
    section_a(&mut report);

    println!("Section B: FFT Convention Inspection...");
    section_b(&mut report);

    println!("Section C: Round-Trip FFT...");
    section_c(&mut report);

    println!("Section D: Madelung Velocity Extraction...");
    section_d(&mut report);

    println!("Section E: Density Round Trip...");
    section_e(&mut report);

    // Summary.
    writeln!(report, "---\n## Summary\n").unwrap();
    writeln!(
        report,
        "See individual sections above for results and verdicts."
    )
    .unwrap();

    std::fs::create_dir_all("data").ok();
    std::fs::write("data/diagnostics-pass2-report.md", &report).unwrap();
    println!("\nReport written to data/diagnostics-pass2-report.md");
}
