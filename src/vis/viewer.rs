//! Interactive 3D particle viewer via kiss3d.

use crate::colormap::particle_density_colors;
use crate::physics::grid::Grid;
use crate::physics::particles::Particles;

/// Open an interactive 3D window showing the particle distribution.
///
/// Particles are rendered as colored points on a dark background.
/// Color encodes local density on a hot colormap (log-scaled).
/// The camera orbits freely with mouse controls.
pub fn render_particles_3d(particles: &Particles, grid: &Grid) {
    use kiss3d::light::Light;
    use kiss3d::nalgebra::Point3;
    use kiss3d::window::Window;

    let mut window = Window::new_with_size("hermes — particle viewer", 1200, 900);
    window.set_background_color(0.0, 0.0, 0.0);
    window.set_light(Light::StickToCamera);
    window.set_point_size(2.0);

    let colors = particle_density_colors(particles, grid);

    let scale = 1.0 / grid.box_length as f32;

    while window.render() {
        for (p, color) in colors.iter().enumerate() {
            let pos = particles.position_components(p);
            let point = Point3::new(
                pos[0] as f32 * scale - 0.5,
                pos[1] as f32 * scale - 0.5,
                pos[2] as f32 * scale - 0.5,
            );
            window.draw_point(&point, &Point3::new(color[0], color[1], color[2]));
        }
    }
}
