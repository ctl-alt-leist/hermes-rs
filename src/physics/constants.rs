/// Physical constants in working units: kpc, M_☉, Gyr, eV (k_B = 1).
///
/// Ported from the plexis constants module to maintain cross-project consistency.
use std::f64::consts::PI;

// ============================================================================
// Fundamental constants
// ============================================================================

/// Gravitational constant (kpc³ M_☉⁻¹ Gyr⁻²).
pub const G: f64 = 4.4996e-6;

/// Speed of light (kpc / Gyr).
pub const C: f64 = 3.0660e5;

/// Proton mass in eV / (kpc / Gyr)².
///
/// Bridges thermal (eV) and kinetic (kpc / Gyr) scales.
pub const PROTON_MASS: f64 = 9.9809e-3;

/// Proton mass in M_☉.
pub const PROTON_MASS_MSUN: f64 = 8.4096e-58;

// ============================================================================
// Unit conversions
// ============================================================================

/// 1 eV = 11604.5 K.
pub const EV_TO_KELVIN: f64 = 11604.5;

/// 1 Mpc = 1000 kpc.
pub const KPC_PER_MPC: f64 = 1000.0;

/// 1 kpc / Gyr ≈ 0.978 km / s.
pub const KPC_GYR_TO_KMS: f64 = 0.97779;

// ============================================================================
// Cosmological constants
// ============================================================================

/// 100 km / s / Mpc in Gyr⁻¹ (≈ 0.1023). Multiply by h to get H₀.
pub const H100: f64 = 100.0 / (KPC_GYR_TO_KMS * KPC_PER_MPC);

/// ρ_c = CRITICAL_DENSITY_FACTOR × H². Units: M_☉ Gyr² / kpc³.
pub const CRITICAL_DENSITY_FACTOR: f64 = 3.0 / (8.0 * PI * G);

// ============================================================================
// Astrophysical reference values
// ============================================================================

/// Virial overdensity parameter (spherical overdensity definition).
pub const DELTA_VIR: f64 = 200.0;

/// Mean molecular weight for fully ionized primordial gas (H + He).
pub const MU_PRIMORDIAL: f64 = 0.59;

// ============================================================================
// CGS boundary constants
//
// Used at the boundary between working units and CGS when interfacing
// with tabulated data (cooling functions, cross sections).
// ============================================================================

/// 1 Gyr in seconds.
pub const GYR_IN_SECONDS: f64 = 3.15576e16;

/// 1 M_☉ in grams.
pub const MSUN_IN_GRAMS: f64 = 1.98892e33;

/// Proton mass in grams.
pub const MP_IN_GRAMS: f64 = 1.6726e-24;

/// Boltzmann constant in CGS (erg / K).
pub const KB_CGS: f64 = 1.3807e-16;

/// 1 kpc in cm.
pub const KPC_IN_CM: f64 = 3.08568e21;
