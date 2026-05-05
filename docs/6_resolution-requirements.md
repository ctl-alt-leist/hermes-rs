# Resolution Requirements

The Schrödinger-Poisson formulation introduces a free parameter $\nu = l/m$ with dimensions kpc$^2$/Gyr that controls the regularization of the density field. At finite $\nu$, the kinetic operator acts like a pressure that prevents density features from collapsing below a characteristic scale. For the formulation to reproduce classical (pressureless) gravitational dynamics, this regularization must be invisible at the scales of interest. This document derives the resolution constraints and tabulates them for the simulation scenarios in hermes.

## The two constraints

The Madelung decomposition of the Schrödinger equation produces an Euler equation with an additional term:

$$
\frac{\partial \mathbf{v}}{\partial t} + (\mathbf{v} \cdot \nabla)\mathbf{v} = -\nabla\Phi + \frac{\nu^2}{2} \ \nabla\left(\frac{\nabla^2\sqrt{\rho}}{\sqrt{\rho}}\right)
$$

The last term is the regularization pressure. It acts on scales comparable to the de Broglie wavelength $\lambda_{dB} = 2\pi\nu / v$, where $v$ is the characteristic velocity. For this pressure to be negligible at a structure of size $r$, we need $\lambda_{dB} \ll r$:

$$
\nu \ll \frac{v \ r}{2\pi}
$$

The grid must simultaneously resolve the de Broglie wavelength. The phase of α changes by $v \ \Delta x / \nu$ per cell, and this must stay below $\pi$ for the discrete representation to be faithful:

$$
\Delta x < \frac{\pi\nu}{v}, \quad n > \frac{L \ v}{\pi\nu}
$$

## The cancellation

Setting $\nu_{max} = vr / (4\pi)$ (a factor of 2 below the borderline) and substituting into the Nyquist constraint:

$$
n > \frac{L \ v}{\pi \ vr / (4\pi)} = \frac{4L}{r}
$$

The velocity $v$ and diffusivity $\nu$ cancel. The grid requirement reduces to the geometric condition that the smallest structure of interest is sampled by at least 4 cells per side. This is the Nyquist sampling theorem in disguise: the field formulation, when configured to suppress its artificial pressure, requires the same spatial resolution as any direct discretization. The de Broglie machinery constrains $\nu$ to be small enough that the regularization is invisible, but it does not impose an additional resolution penalty beyond what the geometry already demands.

The factor of 4 reflects the choice of pressure-suppression strictness. Using $\nu_{max} = vr/(2\pi)$ exactly gives $n = 2L/r$; a more conservative $\nu_{max} = vr/(8\pi)$ gives $n = 8L/r$. Four cells per feature is a reasonable middle ground for structures that are approximately Gaussian.

## Scenario table

For each scenario, $L$ is the box size (chosen large enough to contain the structure with buffer against periodic images), $r$ is the smallest feature to resolve, and $v$ is the characteristic bulk velocity. The grid size $n$ is rounded up to the next power of 2 for FFT efficiency. Memory assumes 16 bytes per complex grid value (double precision).

| Scenario | $r$ | $v$ (kpc/Gyr) | $L$ | $\nu$ (kpc$^2$/Gyr) | $n$ | $n^3$ | Memory |
|---|---|---|---|---|---|---|---|
| Solar system | 1 AU | 31 | 100 AU | $1.2 \times 10^{-8}$ | 512 | $1.3 \times 10^8$ | 2.1 GB |
| Galaxy (spiral arms) | 1 kpc | 205 | 50 kpc | 16 | 256 | $1.7 \times 10^7$ | 270 MB |
| Galaxy group (halo cores) | 50 kpc | 307 | 1 Mpc | 1,220 | 128 | $2.1 \times 10^6$ | 34 MB |
| Cluster (subhalos) | 100 kpc | 1,023 | 5 Mpc | 8,140 | 256 | $1.7 \times 10^7$ | 270 MB |
| Cosmic web (filaments) | 1 Mpc | 511 | 100 Mpc | 40,700 | 512 | $1.3 \times 10^8$ | 2.1 GB |
| Cosmic web tendril | 0.5 Mpc | 205 | 20 Mpc | 8,160 | 256 | $1.7 \times 10^7$ | 270 MB |

All scenarios are feasible at $n \leq 512$ with sub-GB to few-GB memory. The field formulation does not impose a punishing resolution cost for these scales once $\nu$ is chosen correctly.

## Implications for hermes scenes

### Galaxy group

The current galaxy-group-field scene uses $L = 8$ Mpc, $n = 64$, $\nu = 8000$ kpc$^2$/Gyr. The table above calls for $L = 1$ Mpc, $n = 128$, $\nu \approx 1200$ kpc$^2$/Gyr. The current configuration violates both constraints: $\Delta x = 125$ kpc is too coarse to resolve 50 kpc halo cores, and $\nu$ is about 7$\times$ the upper bound. The halos are deep in the regularization-dominated regime, which is why they exhibit solitonic bouncing rather than classical merger dynamics.

### Cosmic web

The cosmic-web-field scene resolves Mpc-scale filaments, where $\nu \approx 40{,}000$ kpc$^2$/Gyr is acceptable and $n = 64$ on a 10 Mpc box gives $\Delta x \approx 160$ kpc. This is adequate for filament-scale structure but would not resolve individual halo cores within the web. For the cosmic web, the current parameterization is in a reasonable regime.

### General rule

The relationship $n \sim 4L/r$ depends only on the geometric ratio of box size to smallest feature. Every continuum method faces this same Nyquist requirement. The field formulation makes the constraint explicit through $\nu$, where finite-difference fluid solvers hide their numerical dissipation in implicit constants. When $\nu_{max}$ and $\nu_{min}$ (from the Nyquist constraint) become incompatible at affordable grid sizes, the Schrödinger formulation is structurally wrong for that problem and a classical transport approach should be used instead.
