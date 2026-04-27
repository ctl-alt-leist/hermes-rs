//! Periodic cubic grid for the simulation domain.
//!
//! A comoving box of side `box_length` divided into `n_cells`³ equal cubic
//! cells. All spatial operations respect periodic boundary conditions.

/// Periodic cubic grid geometry.
#[derive(Debug, Clone)]
pub struct Grid {
    /// Number of cells per side (total cells = n_cells³).
    pub n_cells: usize,
    /// Comoving box side length (kpc).
    pub box_length: f64,
    /// Cell side length h = box_length / n_cells (kpc).
    pub cell_length: f64,
}

impl Grid {
    /// Create a grid with `n_cells` per side in a box of `box_length` kpc.
    pub fn new(n_cells: usize, box_length: f64) -> Self {
        let cell_length = box_length / n_cells as f64;

        Self {
            n_cells,
            box_length,
            cell_length,
        }
    }

    /// Total number of cells in the grid.
    pub fn total_cells(&self) -> usize {
        self.n_cells * self.n_cells * self.n_cells
    }

    /// Volume of a single cell (kpc³).
    pub fn cell_volume(&self) -> f64 {
        self.cell_length * self.cell_length * self.cell_length
    }

    /// Total volume of the box (kpc³).
    pub fn box_volume(&self) -> f64 {
        self.box_length * self.box_length * self.box_length
    }

    /// Wrap a cell index into [0, n_cells) with periodic boundaries.
    pub fn wrap_index(&self, m: isize) -> usize {
        let n = self.n_cells as isize;

        ((m % n + n) % n) as usize
    }

    /// Wrap a position coordinate into [0, box_length) with periodic boundaries.
    pub fn wrap_position(&self, x: f64) -> f64 {
        let wrapped = x % self.box_length;
        if wrapped < 0.0 {
            wrapped + self.box_length
        } else {
            wrapped
        }
    }

    /// Wrap a 3D position vector into the periodic box.
    pub fn wrap_position_3d(&self, position: &mut [f64; 3]) {
        for component in position.iter_mut() {
            *component = self.wrap_position(*component);
        }
    }

    /// Cell center position for integer index triple (m0, m1, m2).
    pub fn cell_center(&self, m0: usize, m1: usize, m2: usize) -> [f64; 3] {
        [
            (m0 as f64 + 0.5) * self.cell_length,
            (m1 as f64 + 0.5) * self.cell_length,
            (m2 as f64 + 0.5) * self.cell_length,
        ]
    }

    /// Cell center position as a morphis grade-1 vector.
    pub fn cell_center_vector(
        &self,
        m0: usize,
        m1: usize,
        m2: usize,
    ) -> morphis::vector::Vector<3> {
        let c = self.cell_center(m0, m1, m2);

        crate::algebra::vector_from_components(c[0], c[1], c[2])
    }

    /// Wavevector component for Fourier index m along one axis.
    ///
    /// Returns k = 2π m / box_length, with the convention that indices
    /// above n_cells/2 wrap to negative frequencies.
    pub fn wavevector_component(&self, m: usize) -> f64 {
        let n = self.n_cells;
        let freq = if m <= n / 2 {
            m as f64
        } else {
            m as f64 - n as f64
        };

        2.0 * std::f64::consts::PI * freq / self.box_length
    }
}
