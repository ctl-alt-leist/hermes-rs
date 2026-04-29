mod plots;
mod viewer;

pub use crate::colormap::{colormap_hot, particle_density_colors};
pub use plots::{plot_conservation, plot_power_spectrum, render_density_slice};
pub use viewer::render_particles_3d;
