//! Shared R2C FFT helpers for 3D real fields.
//!
//! Forward: real [N,N,N] → complex [N,N,N/2+1] (unnormalized).
//! Inverse: complex [N,N,N/2+1] → real [N,N,N] (includes 1/N³ normalization).

use ndarray::{Array3, ArrayD};
use ndrustfft::{FftHandler, R2cFftHandler, ndfft, ndfft_r2c, ndifft, ndifft_r2c};
use num_complex::Complex64;

/// Forward 3D R2C FFT: real [N,N,N] → complex [N,N,N/2+1].
pub fn fft_3d(data: &Array3<f64>, n: usize) -> Array3<Complex64> {
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

/// Inverse 3D C2R FFT: complex [N,N,N/2+1] → real [N,N,N].
pub fn ifft_3d(complex: &Array3<Complex64>, n: usize) -> Array3<f64> {
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

/// Forward 3D R2C FFT from ArrayD (reshapes internally).
pub fn fft_3d_dyn(data: &ArrayD<f64>, n: usize) -> Array3<Complex64> {
    let data_3d = data
        .view()
        .into_shape_with_order((n, n, n))
        .expect("data shape mismatch for fft_3d_dyn")
        .to_owned();
    fft_3d(&data_3d, n)
}

/// Inverse 3D C2R FFT returning ArrayD.
pub fn ifft_3d_dyn(complex: &Array3<Complex64>, n: usize) -> ArrayD<f64> {
    ifft_3d(complex, n).into_dyn()
}
