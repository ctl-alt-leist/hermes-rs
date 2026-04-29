//! Fuzzy dark matter initial conditions.
//!
//! Zel'dovich initialization converted to wavefunction form via the
//! inverse Madelung transformation:
//!
//!   alpha(x) = sqrt(rho(x) / m) * exp(I * m * v(x) . x / ell)

use morphis::even_field::EvenField;
use morphis::metric;
use ndarray::Array3;
use ndrustfft::{FftHandler, R2cFftHandler, ndfft, ndfft_r2c, ndifft, ndifft_r2c};
use num_complex::Complex64;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::core::content::FieldParams;
use crate::error::HermesError;
use crate::physics::cosmology::Cosmology;
use crate::physics::grid::Grid as HermesGrid;
use crate::scenes::cosmic_web::init::power_spectrum;

/// Generate Zel'dovich initial conditions as a wavefunction.
pub fn zeldovich_wavefunction(
    grid: &HermesGrid,
    cosmology: &Cosmology,
    params: &FieldParams,
    scale_factor_initial: f64,
    seed: u64,
) -> Result<EvenField<3>, HermesError> {
    let n = grid.n_cells;
    let n_complex = n / 2 + 1;
    let box_length = grid.box_length;
    let cell_length = grid.cell_length;

    let ell = params.smoothing_length;
    let mass = params.mass_alpha;

    let density_mean = cosmology.density_matter();
    let growth = cosmology.growth_factor(scale_factor_initial);
    let growth_rate = cosmology.growth_rate(scale_factor_initial);
    let hubble_a = cosmology.hubble_parameter(scale_factor_initial);

    // Generate overdensity in Fourier space.
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let mut delta_hat = Array3::<Complex64>::zeros((n, n, n_complex));

    let k_to_hmpc = 1000.0 / cosmology.hubble;
    let volume_box = box_length * box_length * box_length;

    for m0 in 0..n {
        let kx = grid.wavevector_component(m0);
        for m1 in 0..n {
            let ky = grid.wavevector_component(m1);
            for m2 in 0..n_complex {
                let kz = grid.wavevector_component(m2);
                let k2 = kx * kx + ky * ky + kz * kz;

                if k2 < 1e-30 {
                    delta_hat[[m0, m1, m2]] = Complex64::new(0.0, 0.0);
                    continue;
                }

                let k_hmpc = k2.sqrt() * k_to_hmpc;
                let power = power_spectrum(k_hmpc, cosmology);
                let volume_hmpc3 = volume_box / k_to_hmpc.powi(3);
                let sigma = (power / volume_hmpc3).sqrt();

                let re: f64 = rng.random_range(-1.0..1.0);
                let im: f64 = rng.random_range(-1.0..1.0);
                delta_hat[[m0, m1, m2]] = Complex64::new(re * sigma, im * sigma);
            }
        }
    }

    // Compute displacement field: Psi_d(k) = i k_d / k^2 * delta(k).
    let mut displacement: [Array3<f64>; 3] = std::array::from_fn(|_| Array3::zeros((n, n, n)));

    for (d, displacement_d) in displacement.iter_mut().enumerate() {
        let mut psi_hat = Array3::<Complex64>::zeros((n, n, n_complex));

        for m0 in 0..n {
            let kx = grid.wavevector_component(m0);
            for m1 in 0..n {
                let ky = grid.wavevector_component(m1);
                for m2 in 0..n_complex {
                    let kz = grid.wavevector_component(m2);
                    let k2 = kx * kx + ky * ky + kz * kz;
                    if k2 < 1e-30 {
                        continue;
                    }

                    let kd = match d {
                        0 => kx,
                        1 => ky,
                        _ => kz,
                    };
                    psi_hat[[m0, m1, m2]] = Complex64::new(0.0, kd / k2) * delta_hat[[m0, m1, m2]];
                }
            }
        }

        *displacement_d = ifft_3d(&psi_hat, n);
    }

    // Compute overdensity: delta = -D+(a) * div(Psi).
    let mut div_psi = Array3::<f64>::zeros((n, n, n));
    #[allow(clippy::needless_range_loop)]
    for d in 0..3 {
        let psi_d_hat = fft_3d(&displacement[d], n);
        let mut deriv_hat = Array3::<Complex64>::zeros((n, n, n_complex));

        for m0 in 0..n {
            let kx = grid.wavevector_component(m0);
            for m1 in 0..n {
                let ky = grid.wavevector_component(m1);
                for m2 in 0..n_complex {
                    let kz = grid.wavevector_component(m2);
                    let kd = match d {
                        0 => kx,
                        1 => ky,
                        _ => kz,
                    };
                    deriv_hat[[m0, m1, m2]] = psi_d_hat[[m0, m1, m2]] * Complex64::new(0.0, kd);
                }
            }
        }

        let deriv_real = ifft_3d(&deriv_hat, n);
        div_psi += &deriv_real;
    }

    // Build wavefunction via inverse Madelung transform.
    let velocity_factor = growth * growth_rate * hubble_a;
    let morphis_grid = morphis::grid::Grid::<3>::new(n, box_length);
    let g = metric::euclidean::<3>();

    let result = EvenField::from_fn(&morphis_grid, g, |x| {
        let m0 = ((x[0] / cell_length) as usize).min(n - 1);
        let m1 = ((x[1] / cell_length) as usize).min(n - 1);
        let m2 = ((x[2] / cell_length) as usize).min(n - 1);

        let delta = -growth * div_psi[[m0, m1, m2]];
        let rho = density_mean * (1.0 + delta).max(1e-10);

        let vx = velocity_factor * displacement[0][[m0, m1, m2]];
        let vy = velocity_factor * displacement[1][[m0, m1, m2]];
        let vz = velocity_factor * displacement[2][[m0, m1, m2]];

        // Phase: S = m * v . x / ell
        let phase = mass * (vx * x[0] + vy * x[1] + vz * x[2]) / ell;
        let amplitude = (rho / mass).sqrt();

        (amplitude * phase.cos(), amplitude * phase.sin())
    });

    Ok(result)
}

fn fft_3d(data: &Array3<f64>, n: usize) -> Array3<Complex64> {
    let n_complex = n / 2 + 1;
    let handler_r2c = R2cFftHandler::new(n);
    let handler_c2c_1 = FftHandler::new(n);
    let handler_c2c_0 = FftHandler::new(n);

    let mut complex = Array3::<Complex64>::zeros((n, n, n_complex));
    ndfft_r2c(data, &mut complex, &handler_r2c, 2);
    let mut scratch = complex.clone();
    ndfft(&complex, &mut scratch, &handler_c2c_1, 1);
    complex.assign(&scratch);
    ndfft(&complex, &mut scratch, &handler_c2c_0, 0);
    complex.assign(&scratch);
    complex
}

fn ifft_3d(complex: &Array3<Complex64>, n: usize) -> Array3<f64> {
    let handler_c2c_0 = FftHandler::new(n);
    let handler_c2c_1 = FftHandler::new(n);
    let handler_r2c = R2cFftHandler::new(n);

    let mut work = complex.clone();
    let mut scratch = work.clone();
    ndifft(&work, &mut scratch, &handler_c2c_0, 0);
    work.assign(&scratch);
    ndifft(&work, &mut scratch, &handler_c2c_1, 1);
    work.assign(&scratch);

    let mut real = Array3::<f64>::zeros((n, n, n));
    ndifft_r2c(&work, &mut real, &handler_r2c, 2);
    real
}
