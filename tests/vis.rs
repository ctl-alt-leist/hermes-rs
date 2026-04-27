#![cfg(feature = "vis")]

use hermes_rs::cosmology::planck_2018;
use hermes_rs::grid::Grid;
use hermes_rs::initial::zeldovich_init;
use hermes_rs::vis::{colormap_hot, particle_density_colors};

// ============================================================================
// Colormap
// ============================================================================

#[test]
fn colormap_black_at_zero() {
    let [r, g, b] = colormap_hot(0.0);
    assert!(
        r < 0.01 && g < 0.01 && b < 0.01,
        "colormap(0) should be near black"
    );
}

#[test]
fn colormap_bright_at_one() {
    let [r, g, b] = colormap_hot(1.0);
    assert!(r > 0.9 && g > 0.9, "colormap(1) should be near white");
}

#[test]
fn colormap_monotone_brightness() {
    let brightness = |t: f64| -> f32 {
        let [r, g, b] = colormap_hot(t);
        r + g + b
    };

    let brightness_low = brightness(0.1);
    let brightness_mid = brightness(0.5);
    let brightness_high = brightness(0.9);

    assert!(
        brightness_mid > brightness_low,
        "brightness should increase with value"
    );
    assert!(
        brightness_high > brightness_mid,
        "brightness should increase with value"
    );
}

#[test]
fn colormap_clamps_out_of_range() {
    let below = colormap_hot(-0.5);
    let above = colormap_hot(1.5);

    assert_eq!(
        below,
        colormap_hot(0.0),
        "negative values should clamp to 0"
    );
    assert_eq!(above, colormap_hot(1.0), "values > 1 should clamp to 1");
}

// ============================================================================
// Particle density colors
// ============================================================================

#[test]
fn density_colors_correct_count() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();

    let colors = particle_density_colors(&particles, &grid);

    assert_eq!(colors.len(), particles.count());
}

#[test]
fn density_colors_valid_rgb() {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();

    let colors = particle_density_colors(&particles, &grid);

    for (p, [r, g, b]) in colors.iter().enumerate() {
        assert!(
            (0.0..=1.0).contains(r) && (0.0..=1.0).contains(g) && (0.0..=1.0).contains(b),
            "particle {p} has invalid color: [{r}, {g}, {b}]"
        );
    }
}

#[test]
fn density_colors_vary_for_nonuniform() {
    let grid = Grid::new(16, 100_000.0);
    let cosmology = planck_2018();
    let particles = zeldovich_init(16, &grid, &cosmology, 0.02, 42).unwrap();

    let colors = particle_density_colors(&particles, &grid);

    // Zel'dovich ICs are non-uniform — colors should not all be identical.
    let first = colors[0];
    let any_different = colors.iter().any(|c| (c[0] - first[0]).abs() > 0.01);

    assert!(
        any_different,
        "density colors should vary for non-uniform particles"
    );
}
