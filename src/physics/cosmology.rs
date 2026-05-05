/// FLRW cosmological background.
///
/// Standard ΛCDM cosmology with the Friedmann parameterization:
///
/// ```text
/// H(a) = H₀ √[ Ω_Λ + Ω_k a⁻² + Ω_m a⁻³ + Ω_r a⁻⁴ ]
/// ```
///
/// with the constraint Ω_Λ + Ω_k + Ω_m + Ω_r = 1.
///
/// All methods take the scale factor `a` as the time variable. The scale
/// factor is related to redshift by a = 1/(1+z).
use serde::Deserialize;

use crate::error::HermesError;
use crate::physics::constants::{CRITICAL_DENSITY_FACTOR, H100};

// ============================================================================
// Cosmology struct
// ============================================================================

/// ΛCDM cosmological parameters.
///
/// Parameterized by the four density fractions (Ω_Λ, Ω_m, Ω_k, Ω_r),
/// the dimensionless Hubble parameter h, the baryon density Ω_b, the
/// amplitude of matter fluctuations σ₈, and the spectral index n_s.
#[derive(Debug, Clone, Deserialize)]
pub struct Cosmology {
    /// Dimensionless Hubble parameter h (H₀ = 100h km/s/Mpc).
    pub hubble: f64,
    /// Total matter density parameter Ω_m (CDM + baryons).
    pub omega_m: f64,
    /// Baryon density parameter Ω_b (subset of Ω_m).
    pub omega_b: f64,
    /// Radiation density parameter Ω_r (photons + relativistic neutrinos).
    pub omega_r: f64,
    /// Spatial curvature parameter Ω_k.
    pub omega_k: f64,
    /// Vacuum energy density parameter Ω_v (cosmological constant).
    pub omega_v: f64,
    /// RMS matter fluctuation amplitude in 8 h⁻¹ Mpc spheres.
    pub sigma_8: f64,
    /// Spectral index of the primordial power spectrum.
    pub spectral_index: f64,
}

impl Cosmology {
    /// Validate that cosmological parameters are physically consistent.
    pub fn validate(&self) -> Result<(), HermesError> {
        if self.hubble <= 0.0 {
            return Err(HermesError::Cosmology(format!(
                "hubble must be positive, got {}",
                self.hubble
            )));
        }
        if self.omega_m <= 0.0 {
            return Err(HermesError::Cosmology(format!(
                "omega_m must be positive, got {}",
                self.omega_m
            )));
        }
        if self.omega_b <= 0.0 {
            return Err(HermesError::Cosmology(format!(
                "omega_b must be positive, got {}",
                self.omega_b
            )));
        }
        if self.omega_b > self.omega_m {
            return Err(HermesError::Cosmology(format!(
                "omega_b ({}) exceeds omega_m ({})",
                self.omega_b, self.omega_m
            )));
        }
        if self.omega_v < 0.0 {
            return Err(HermesError::Cosmology(format!(
                "omega_v must be non-negative, got {}",
                self.omega_v
            )));
        }
        if self.omega_r < 0.0 {
            return Err(HermesError::Cosmology(format!(
                "omega_r must be non-negative, got {}",
                self.omega_r
            )));
        }
        if self.sigma_8 <= 0.0 {
            return Err(HermesError::Cosmology(format!(
                "sigma_8 must be positive, got {}",
                self.sigma_8
            )));
        }

        let total = self.omega_v + self.omega_k + self.omega_m + self.omega_r;
        if (total - 1.0).abs() > 1e-6 {
            return Err(HermesError::Cosmology(format!(
                "density parameters must sum to 1: Ω_Λ + Ω_k + Ω_m + Ω_r = {total:.10}"
            )));
        }

        Ok(())
    }

    // ========================================================================
    // Derived parameters
    // ========================================================================

    /// Hubble constant H₀ in Gyr⁻¹.
    pub fn hubble_constant(&self) -> f64 {
        self.hubble * H100
    }

    /// Hubble time t_H = 1/H₀ in Gyr.
    pub fn hubble_time(&self) -> f64 {
        1.0 / self.hubble_constant()
    }

    /// Cold dark matter density parameter Ω_cdm = Ω_m − Ω_b.
    pub fn omega_cdm(&self) -> f64 {
        self.omega_m - self.omega_b
    }

    /// Cosmic baryon fraction f_b = Ω_b / Ω_m.
    pub fn baryon_fraction(&self) -> f64 {
        self.omega_b / self.omega_m
    }

    /// True if the universe is spatially flat (|Ω_k| < 10⁻⁶).
    pub fn is_flat(&self) -> bool {
        self.omega_k.abs() < 1e-6
    }

    // ========================================================================
    // Scale-factor-dependent quantities
    // ========================================================================

    /// Dimensionless Hubble parameter E(a) = H(a) / H₀.
    ///
    /// From the Friedmann equation:
    ///
    /// ```text
    /// E(a) = √[ Ω_Λ + Ω_k a⁻² + Ω_m a⁻³ + Ω_r a⁻⁴ ]
    /// ```
    pub fn expansion_rate(&self, a: f64) -> f64 {
        let a2 = a * a;
        let a3 = a2 * a;
        let a4 = a3 * a;

        (self.omega_v + self.omega_k / a2 + self.omega_m / a3 + self.omega_r / a4).sqrt()
    }

    /// Hubble parameter H(a) in Gyr⁻¹.
    pub fn hubble_parameter(&self, a: f64) -> f64 {
        self.hubble_constant() * self.expansion_rate(a)
    }

    /// Critical density ρ_c(a) = 3H(a)² / (8πG) in M_☉ / kpc³.
    pub fn density_critical(&self, a: f64) -> f64 {
        let h_a = self.hubble_parameter(a);

        CRITICAL_DENSITY_FACTOR * h_a * h_a
    }

    /// Mean comoving matter density ρ̄_m(a) in M_☉ / kpc³.
    ///
    /// In comoving coordinates this is constant: Ω_m × ρ_c(a=1).
    pub fn density_matter(&self) -> f64 {
        self.omega_m * self.density_critical(1.0)
    }

    /// Cosmic age at scale factor `a` in Gyr.
    ///
    /// Numerical integration via midpoint quadrature:
    ///
    /// ```text
    /// t(a) = (1/H₀) ∫₀ᵃ da' / [a' E(a')]
    /// ```
    pub fn cosmic_time(&self, a: f64) -> f64 {
        let a_floor = 1e-8;
        if a <= a_floor {
            return 0.0;
        }

        let n_steps = 500;
        let da = (a - a_floor) / n_steps as f64;
        let mut integral = 0.0;
        for n in 0..n_steps {
            let a_mid = a_floor + (n as f64 + 0.5) * da;
            integral += da / (a_mid * self.expansion_rate(a_mid));
        }

        self.hubble_time() * integral
    }

    /// Linear growth factor D₊(a), normalized so D₊(1) = 1.
    ///
    /// Uses the Carroll, Press & Turner (1992) approximation, which is
    /// accurate to better than 1% for flat and open ΛCDM models.
    pub fn growth_factor(&self, a: f64) -> f64 {
        let d_unnorm = growth_factor_unnormalized(self, a);
        let d_0 = growth_factor_unnormalized(self, 1.0);

        d_unnorm / d_0
    }

    /// Logarithmic growth rate f(a) = d ln D₊ / d ln a.
    ///
    /// Evaluated by finite difference on the normalized growth factor.
    pub fn growth_rate(&self, a: f64) -> f64 {
        let epsilon = 1e-4 * a;
        let d_plus = self.growth_factor(a + epsilon);
        let d_minus = self.growth_factor(a - epsilon);

        (a / self.growth_factor(a)) * (d_plus - d_minus) / (2.0 * epsilon)
    }

    /// Kick factor for symplectic integration: ∫ da / [a² H(a)].
    ///
    /// Converts a momentum impulse from the force at scale factor `a`
    /// into the corresponding change in canonical momentum over the
    /// interval [a_start, a_end].
    pub fn kick_factor(&self, a_start: f64, a_end: f64) -> f64 {
        integrate(a_start, a_end, |a| 1.0 / (a * a * self.hubble_parameter(a)))
    }

    /// Drift factor for symplectic integration: ∫ da / [a³ H(a)].
    ///
    /// Converts a canonical momentum into the corresponding position
    /// displacement over the interval [a_start, a_end].
    pub fn drift_factor(&self, a_start: f64, a_end: f64) -> f64 {
        integrate(a_start, a_end, |a| {
            1.0 / (a * a * a * self.hubble_parameter(a))
        })
    }
}

// ============================================================================
// Conversion from EngineConfig
// ============================================================================

impl Cosmology {
    /// Build a Cosmology from the new EngineConfig spacetime parameters.
    ///
    /// The EngineConfig stores H₀ in km/s/Mpc (e.g. 67.4); we convert
    /// to the dimensionless h (0.674) used internally.
    pub fn from_engine_config(config: &crate::config::EngineConfig) -> Result<Self, HermesError> {
        let spacetime = &config.ontology.spacetime;

        let hubble_kms = spacetime.hubble.ok_or_else(|| {
            HermesError::Config("cosmology requires hubble parameter".to_string())
        })?;

        Ok(Self {
            hubble: hubble_kms / 100.0,
            omega_m: spacetime.omega_m.unwrap_or(0.315),
            omega_b: spacetime.omega_b.unwrap_or(0.0493),
            omega_r: spacetime.omega_r.unwrap_or(9.15e-5),
            omega_k: spacetime.omega_k.unwrap_or(0.0),
            omega_v: spacetime.omega_v.unwrap_or(0.685),
            sigma_8: spacetime.sigma_8.unwrap_or(0.811),
            spectral_index: spacetime.spectral_index.unwrap_or(0.965),
        })
    }
}

// ============================================================================
// Factory functions
// ============================================================================

/// Planck 2018 best-fit ΛCDM cosmology.
pub fn planck_2018() -> Cosmology {
    Cosmology {
        hubble: 0.674,
        omega_m: 0.315,
        omega_b: 0.0493,
        omega_r: 9.15e-5,
        omega_k: 0.0,
        omega_v: 0.6849085,
        sigma_8: 0.811,
        spectral_index: 0.965,
    }
}

/// Einstein-de Sitter cosmology (Ω_m = 1, flat, no Λ).
///
/// All quantities have closed-form analytic solutions, making this
/// the primary test case for numerical methods:
///
///   - E(a) = a⁻³ᐟ²
///   - t(a) = (2/3) a³ᐟ² / H₀
///   - D₊(a) = a
///   - f(a) = 1
pub fn einstein_de_sitter() -> Cosmology {
    Cosmology {
        hubble: 0.7,
        omega_m: 1.0,
        omega_b: 0.05,
        omega_r: 0.0,
        omega_k: 0.0,
        omega_v: 0.0,
        sigma_8: 0.811,
        spectral_index: 1.0,
    }
}

// ============================================================================
// Internal helpers
// ============================================================================

/// Carroll, Press & Turner (1992) growth factor (unnormalized).
fn growth_factor_unnormalized(cosmo: &Cosmology, a: f64) -> f64 {
    let omega_a = cosmo.omega_m / (cosmo.omega_m + cosmo.omega_v * a * a * a);
    let lambda_a = 1.0 - omega_a;

    (5.0 / 2.0) * omega_a * a
        / (omega_a.powf(4.0 / 7.0) - lambda_a + (1.0 + omega_a / 2.0) * (1.0 + lambda_a / 70.0))
}

/// Midpoint-rule numerical integration of f(a) over [a_start, a_end].
fn integrate<F: Fn(f64) -> f64>(a_start: f64, a_end: f64, f: F) -> f64 {
    let n_steps = 500;
    let da = (a_end - a_start) / n_steps as f64;
    let mut result = 0.0;
    for n in 0..n_steps {
        let a = a_start + (n as f64 + 0.5) * da;
        result += da * f(a);
    }

    result
}
