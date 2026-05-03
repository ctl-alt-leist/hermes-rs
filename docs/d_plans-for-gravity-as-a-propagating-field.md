# Gravity as a Propagating Field

CLAUDE: This whole set of notes is just a mathematical set based on something we hope to do. We aren't going to do this
in our work now. Only once we have a unified framework will we incorporate this field. What you should know, is that we
will want to be able to either use (a) the gravitational method we use now, (b) use this field for gravity. So, this
will need to be a configuration input of the physics engine as well as data that could be stored in the data write. It
is just worth keeping in mind that our new framework should be capable of accepting this toggle between gravitational
approaches. Of course, this isn't really anything new in principle, we want to be able to toggle all sorts of methods
used (or one way to think of it, which Lagrangians we are turning on and off).

The gravitational sector of the framework, as developed in the earlier notes, treats the gravitational potential $Φ$ as a quantity determined at each instant by the matter distribution through the comoving Poisson equation. This is the Newtonian formulation: gravity is an instantaneous response to matter, with no dynamical content of its own. The matter wavefunctions $α$ and $β$ are proper fields with state and evolution; $Φ$ is a derived bookkeeping device.

This document develops the next layer of the framework, in which the gravitational sector is itself a field. The instantaneous potential $Φ(\mathbf{x}, t)$ is replaced by a dynamical scalar $φ(\mathbf{x}, t)$ that obeys a wave equation, propagates at a finite speed $c_g$, carries its own energy and momentum, and reduces to the Newtonian potential in the appropriate static limit. The matter sectors are unaffected at the level of their dynamical equations; what changes is the equation for $φ$ itself, which is no longer an elliptic constraint but a hyperbolic evolution. The framework gains a fourth dynamical sector with the same algorithmic character as the matter sectors: a per-cell state, a spectral propagator, and a Lagrangian-derived equation of motion.

The treatment proceeds from the internal geometry of the gravitational field, through its Lagrangian density and equation of motion, through the free wave behavior and the spectral structure that governs it, through the source coupling and the recovery of retarded Newtonian potentials, through the limits in which the parameter $c_g$ approaches infinity (recovering the original instantaneous formulation) and the speed of light (recovering Nordström scalar gravity), and finally through what the propagating formulation captures and what it does not.

## The Internal Geometry of $φ$

The gravitational field is a scalar:

$$
φ(\mathbf{x}, t) \;\in\; \mathcal{G}^0
$$

This is the simplest grade choice consistent with the physics. The earlier $Φ$ was already a scalar, and elevating it to a dynamical field does not change its grade. What changes is that the field now has a conjugate momentum, also scalar:

$$
π_φ(\mathbf{x}, t) \;=\; ∂_t φ(\mathbf{x}, t) \;\in\; \mathcal{G}^0
$$

The dynamical state of the gravitational sector is therefore the canonical pair $(φ, π_φ)$, both grade-zero, totaling two real-valued fields per cell. This matches the data count of a matter wavefunction in the even subalgebra, but the structure is different. Where $α \in \mathcal{G}^+$ encodes density and phase as scalar plus pseudoscalar in a single complex-like object, the gravitational pair $(φ, π_φ)$ encodes value and time-derivative as two scalars in a symplectic relationship. The matter pair is $(\text{modulus}, \text{phase})$; the gravity pair is $(\text{position}, \text{velocity})$ in the field-configuration sense.

There is a substantive choice not made here. One could form the combination

$$
ψ_φ(\mathbf{x}, t) \;=\; ω_0 \, φ(\mathbf{x}, t) \;+\; \mathbb{1} \, π_φ(\mathbf{x}, t)
$$

for some reference frequency $ω_0$, placing gravity in the even subalgebra alongside $α$ and $β$. This is algebraically possible and exhibits a structural parallel to the matter sectors. It is not adopted here because $ω_0$ has no physical meaning for the gravitational field (gravity has no de Broglie scale and no preferred frequency), and the combination is most naturally diagonal in Fourier space rather than in real space. The clean expression of gravity's even-subalgebra character is per-mode, where $ω_k = c_g |\mathbf{k}|$ provides a natural mode-dependent frequency, and the spectral propagator developed below will display this structure explicitly. In configuration space, the canonical pair is the more transparent representation.

## The Lagrangian

The gravitational Lagrangian density takes the standard scalar-field form, with a kinetic term, a gradient term, and a source coupling to the matter sectors. Adopting the conventions of the earlier notes:

$$
\mathcal{L}_φ \;=\; \frac{1}{8 π G} \left[
	\frac{1}{2 \, c_g^2} \, (∂_t φ)^2
	\;-\; \frac{1}{2} \, (∇ φ) \rfloor (∇ φ)
\right]
\;-\; ρ \, φ
$$

with $ρ$ the total mass density sourcing the field

$$
ρ \;=\; m_α \, \bar{α} α \;+\; m_β \, \bar{β} β \;+\; \mathcal{E}_γ / c^2
$$

including all matter species and the energy-density contribution of the electromagnetic sector. The prefactor $1 / (8 π G)$ normalizes the kinetic and gradient terms to give the Newtonian limit with the conventional Poisson coefficient. The relative sign between the kinetic and gradient terms is the relativistic Lagrangian sign convention: time-derivatives squared positive, spatial-derivatives squared negative, producing the d'Alembertian operator with the right causal structure.

The structural parallel to the dark-matter Lagrangian is exact in form. Where $\mathcal{L}_α$ has an antisymmetric first-order time-derivative term (the Schrödinger structure), $\mathcal{L}_φ$ has a symmetric second-order time-derivative term (the wave-equation structure). Where $\mathcal{L}_α$ has the gradient inner product weighted by $ℓ^2 / (2 m_α)$, $\mathcal{L}_φ$ has the gradient inner product weighted by $1 / (16 π G)$. The two sectors are different in the order of their time evolution, but they share the gradient-term geometry.

The source coupling $-ρ \, φ$ links gravity to all matter sectors simultaneously. It is universal: every species that carries mass contributes to $ρ$ on equal footing, and the resulting $φ$ couples back to each species through the gravitational interaction term in its own equation of motion.

## The Equation of Motion

Variation of the Lagrangian with respect to $φ$ yields the wave equation:

$$
-\frac{1}{c_g^2} \, ∂_t^2 φ \;+\; ∇^2 φ \;=\; 4 π G \, ρ
$$

In the comoving cosmological setting where the source is the overdensity rather than the absolute density, and the spatial Laplacian acts in comoving coordinates, the equation reads

$$
-\frac{1}{c_g^2} \, ∂_t^2 φ \;+\; ∇^2 φ \;=\; 4 π G \, a^2 \, (ρ - \bar{ρ})
$$

with the same factor of $a^2$ that appears in the standard comoving Poisson equation, and with $\bar{ρ}$ the cosmic mean density that gets absorbed into the background metric.

This is the inhomogeneous wave equation. In the limit $c_g \to \infty$, the time-derivative term vanishes and the equation reduces to the ordinary Poisson equation $∇^2 φ = 4 π G \, a^2 \, (ρ - \bar{ρ})$, recovering the original instantaneous formulation as the static limit of a propagating field. For finite $c_g$, the time-derivative term encodes the fact that gravitational signals do not propagate instantaneously: a change in the source at one location takes a finite time to influence the field at another location.

The wave operator $-c_g^{-2} \, ∂_t^2 + ∇^2$ is the d'Alembertian with metric signature $(-, +, +, +)$ and signal speed $c_g$. Its causal structure is the past light cone of the chosen $c_g$: the value of $φ$ at $(\mathbf{x}, t)$ depends on the source on the surface

$$
|\mathbf{x} - \mathbf{x}'| \;=\; c_g \, (t - t')
$$

with $t' < t$, weighted by the retarded Green's function. Information about source changes propagates outward at $c_g$.

## Free Propagation: The Spectral Picture

Set the source to zero and consider the homogeneous wave equation. In Fourier space at each wavenumber $\mathbf{k}$, it becomes a one-dimensional harmonic-oscillator equation:

$$
∂_t^2 \hat{φ}(\mathbf{k}, t) \;+\; ω_\mathbf{k}^2 \, \hat{φ}(\mathbf{k}, t) \;=\; 0
$$

with the per-mode frequency

$$
ω_\mathbf{k} \;=\; c_g \, |\mathbf{k}|
$$

This is the linear dispersion relation, a defining feature of the wave equation as distinct from the Schrödinger equation. The matter sector has $ω(\mathbf{k}) = ν \, |\mathbf{k}|^2 / 2$ (quadratic dispersion); the gravitational sector has $ω(\mathbf{k}) = c_g \, |\mathbf{k}|$ (linear dispersion). Linear dispersion means signals propagate at a fixed speed regardless of wavelength: the group velocity $∂ ω / ∂ k = c_g$ is the same for every mode. This is the property that makes wave-equation propagation a clean model of causality. Quadratic dispersion would produce wavelength-dependent propagation speeds and the spreading of localized disturbances seen in the matter sector.

The general homogeneous solution at each mode is

$$
\hat{φ}(\mathbf{k}, t) \;=\; A_\mathbf{k} \, e^{\mathbb{1} \, ω_\mathbf{k} \, t} \;+\; B_\mathbf{k} \, e^{-\mathbb{1} \, ω_\mathbf{k} \, t}
$$

a superposition of two counter-rotating phases. Equivalently, the canonical pair $(\hat{φ}_\mathbf{k}, \hat{π}_{φ, \mathbf{k}})$ rotates rigidly in its symplectic plane:

$$
\begin{pmatrix} ω_\mathbf{k} \, \hat{φ} \\ \hat{π}_φ \end{pmatrix}_{t + Δt}
\;=\;
\begin{pmatrix}
	\cos(ω_\mathbf{k} Δt) & \sin(ω_\mathbf{k} Δt) \\
	-\sin(ω_\mathbf{k} Δt) & \cos(ω_\mathbf{k} Δt)
\end{pmatrix}
\begin{pmatrix} ω_\mathbf{k} \, \hat{φ} \\ \hat{π}_φ \end{pmatrix}_t
$$

This is the structural payoff of the formulation. Per Fourier mode, free gravitational evolution is a unit-modulus rotation in the canonical plane, exactly analogous to the unit-modulus phase rotation that the matter sectors undergo in their kinetic step. Both sectors have a kinetic operator that is diagonal in Fourier space and unitary, and both can be integrated exactly in the spectral representation by precomputing the per-mode rotation coefficients and applying them as a pointwise multiplication in $\mathbf{k}$-space. The cost per timestep is one R2C FFT pair on each of $φ$ and $π_φ$, plus the pointwise rotation. The gravitational sector inherits the algorithmic shape of the matter sectors despite being a different equation.

The price for this clean structure is a Courant condition. The wave equation has finite signal speed $c_g$, so a timestep that propagates information further than one cell in real space is unstable:

$$
c_g \, Δt \;\leq\; Δx
$$

This is the standard CFL constraint for hyperbolic equations and is qualitatively different from the matter sector, where split-step is unconditionally stable for any $Δt$. For $c_g$ chosen as the speed of light and $Δx \sim$ Mpc, the CFL bound gives $Δt \lesssim 3 \, \text{Myr}$, far below the Gyr-scale cosmological timestep, so subcycling of the gravity solver inside each cosmological step is unavoidable when $c_g$ is set to a physical speed. For $c_g$ chosen artificially as a numerical regulator, the bound is whatever the user picks.

## The Source Term and Retarded Potentials

The inhomogeneous wave equation with a localized source $ρ(\mathbf{x}', t')$ is solved by the retarded Green's function

$$
G_{\text{ret}}(\mathbf{x}, t; \mathbf{x}', t') \;=\; -\frac{δ\!\big( t - t' - |\mathbf{x} - \mathbf{x}'| / c_g \big)}{|\mathbf{x} - \mathbf{x}'|}
$$

giving the explicit solution

$$
φ(\mathbf{x}, t) \;=\; -G \int d^3x' \, \frac{ρ(\mathbf{x}', t - |\mathbf{x} - \mathbf{x}'| / c_g)}{|\mathbf{x} - \mathbf{x}'|}
$$

This is the relativistic retarded Newtonian potential. The field at a spacetime point $(\mathbf{x}, t)$ is the Newtonian-style integral of the source over space, but the source is evaluated at the retarded time $t' = t - |\mathbf{x} - \mathbf{x}'| / c_g$, on the past light cone of the field point. Sources at distance $r$ contribute to the present field with a time delay $r / c_g$.

In the limit $c_g \to \infty$, the retarded time becomes the present time, the integral becomes the instantaneous Newtonian potential, and the framework recovers the original Poisson formulation exactly. For finite $c_g$, the retardation is automatic: changes in the source propagate outward at speed $c_g$, the gravitational influence of a moving mass is the present mass distribution combined with the retardation kernel, and a binary system or a collapsing structure naturally emits gravitational signals that travel through the box at $c_g$.

Two qualitative features of retarded gravity follow. First, the gravitational field of a moving mass is not centered on its instantaneous position but on its retarded position, the position it occupied $r / c_g$ ago. For non-relativistic motion this correction is small, but it is structurally present and accumulates over long propagation distances. Second, accelerating masses radiate. A mass that oscillates produces an oscillating $φ$ in its vicinity, and that oscillation propagates outward as a gravitational wave at speed $c_g$. The wave carries energy and momentum away from the source.

## How Matter Feels the Field

The coupling of $φ$ to the matter sectors is the same as the coupling of $Φ$ to the matter sectors in the earlier formulation: the matter equations contain a term $m \, φ$ in place of $m \, Φ$, and nothing else changes. For dark matter,

$$
ℓ \, \mathbb{1} \, ∂_t α \;=\; -\frac{ℓ^2}{2 \, m_α \, a^2} \, ∇^2 α \;+\; m_α \, φ \, α
$$

and analogously for the baryon and electromagnetic sectors. The matter wavefunctions feel gravity through pointwise multiplication by $φ$, exactly as before. What is different is that $φ$ is now a dynamical field with its own evolution, so the matter at $(\mathbf{x}, t)$ sees a $φ$ that reflects the source distribution on the past light cone rather than the instantaneous source distribution. The matter equation form is unchanged; the $φ$ that appears in it is structurally different.

This is why the propagating-gravity formulation is a clean addition to the framework rather than a rewrite. The matter sectors are unaltered. The Poisson solver is replaced by a wave-equation propagator. The interface between gravity and matter — the substitution of $φ$ for $Φ$ in each matter equation — is identical to before.

A consistency check: in the static limit, the matter sees $φ = Φ$ exactly, the matter equations reduce to the original equations, and the dynamics is unchanged. The propagating formulation is a refinement that activates only when $∂_t ρ$ is large enough for retardation to matter, which on cosmological scales means only when something genuinely time-dependent is happening (collapse, oscillation, structure formation). For slowly-evolving configurations, the propagating and instantaneous formulations agree to high precision.

## Gravitational Radiation

The wave equation has propagating modes, and these modes carry energy. The energy density of the gravitational field follows from the standard scalar-field stress-energy:

$$
\mathcal{E}_φ \;=\; \frac{1}{8 π G} \left[
	\frac{1}{2 \, c_g^2} \, (∂_t φ)^2
	\;+\; \frac{1}{2} \, (∇ φ) \rfloor (∇ φ)
\right]
$$

with the relative signs flipped from the Lagrangian (standard Legendre transform). The first term is the kinetic energy density of the field; the second is the gradient energy density, which is the gravitational analog of an elastic potential energy.

In the absence of sources, this energy is conserved: the wave equation's free solutions carry energy density at rate $\mathcal{E}_φ$ that propagates outward at $c_g$, and the integrated total $\int d^3x \, \mathcal{E}_φ$ is exactly preserved. With sources present, the gravitational field exchanges energy with the matter sectors through the source coupling, and the joint matter-plus-gravity total energy is conserved.

The radiation from accelerating matter is real. A spherical pulsation of an overdensity emits gravitational waves at $c_g$ that carry energy away from the source. This is a feature of scalar relativistic gravity that genuine general relativity does not have: in GR, a spherical pulsation does not radiate (Birkhoff's theorem), because the tensor structure of gravitational radiation excludes monopole modes. Scalar gravity has monopole radiation; tensor gravity does not. This is the central reason scalar relativistic gravity is the wrong theory for the real universe.

For the framework, the consequence is that the propagating gravitational sector includes scalar gravitational waves that real gravity would not. On the timescales and length scales relevant to cosmological structure formation, these effects are negligible (the energy in scalar gravitational radiation is small compared to the energy in matter motion), and the propagating-gravity formulation is a faithful representation of retarded Newtonian gravity. For specific small-scale problems involving rapidly accelerating masses, the scalar radiation is a representation artifact, not a physical signal.

## The Speed Parameter $c_g$

The wave equation's propagation speed $c_g$ is a free parameter. Three regimes are physically meaningful.

**Newtonian limit ($c_g \to \infty$).** The time-derivative term in the wave equation vanishes, the equation reduces to Poisson, and the framework recovers the original instantaneous formulation. This is the parameter regime where the propagating-gravity sector is computationally redundant with the elliptic Poisson solver. It is useful as a verification check: any propagating-gravity simulation should agree with the Poisson formulation in the limit of large $c_g$.

**Physical relativistic limit ($c_g = c$).** The propagation speed is the speed of light. Gravitational signals travel at $c$, retardation is on the same timescale as electromagnetic retardation, and the framework agrees with Nordström scalar gravity, the historical scalar precursor to general relativity. This regime is physically motivated and corresponds to a well-defined classical field theory, although it does not match real gravity in detail (light bending wrong by a factor of two, no Mercury perihelion, no tensor radiation modes). For cosmological-scale structure-formation calculations on boxes of $\sim 10$ Mpc, the predictions of $c_g = c$ scalar gravity are indistinguishable from those of full general relativity, because gravity is everywhere weak and the matter is non-relativistic.

**Numerical regulator regime ($c_g$ moderate but finite).** Choosing $c_g$ such that $c_g \, Δt_{\text{cosmo}} \sim L_{\text{box}}$ makes a gravitational disturbance cross the box in roughly one cosmological timestep. This regime is faster than the true relativistic limit (which would require many cosmological substeps for one box-crossing of light) and is a sensible engineering choice when the goal is to incorporate causality structurally without committing to the actual speed of light. The simulation has finite gravity speed, retarded potentials work correctly, energy conservation includes the field, but the specific value of $c_g$ is chosen for computational convenience rather than for physical accuracy. This is a model choice, not an approximation.

The relationship between these regimes is monotonic. As $c_g$ decreases from infinity, retardation effects appear progressively, structure formation timing shifts slightly (because gravitational signals between distant overdensities arrive late), and gravitational radiation becomes more energetically significant. As $c_g$ increases, all of these effects shrink to zero and the formulation becomes operationally equivalent to instantaneous Poisson.

## Symmetries and Conservation Laws

The Lagrangian $\mathcal{L}_φ$ is invariant under the same translations, rotations, and time-translation that the matter Lagrangians are. Noether's theorem produces three corresponding conservation laws.

Translation invariance gives momentum conservation. The momentum density of the gravitational field is

$$
\mathbf{p}_φ \;=\; \frac{1}{4 π G \, c_g^2} \, (∂_t φ) \, ∇ φ
$$

and the integrated total $\int d^3x \, \mathbf{p}_φ$ is conserved when the integration is taken over a region where no momentum flows across the boundary. Combined with the matter momenta from the matter sectors, the joint total is conserved exactly.

Time translation gives energy conservation. The integrated total energy

$$
E_{\text{tot}} \;=\; \int d^3x \, (\mathcal{E}_α + \mathcal{E}_β + \mathcal{E}_γ + \mathcal{E}_φ)
$$

is conserved exactly when the cosmological background is held fixed, and conserved up to the Layzer-Irvine background-evolution terms when the FLRW expansion is included. The gravitational field's contribution to the total energy is essential here: in the instantaneous-Poisson formulation, the gravitational potential energy is not a separately-stored quantity and is recomputed at each step from the matter distribution, with no clean conservation diagnostic. In the propagating-gravity formulation, $\mathcal{E}_φ$ is a real field-theoretic quantity that the dynamics exactly conserves (modulo source exchange with matter), and total energy conservation is structurally manifest.

Rotation invariance gives angular momentum conservation, and the gravitational contribution to the angular momentum density is the wedge of position with the field momentum density.

These conservation laws are the structural reason the propagating formulation is cleaner than the instantaneous formulation, beyond any considerations of causality or relativistic accuracy. In the instantaneous formulation, the gravitational field is a passive response to matter, with no state of its own and no separately-conserved quantities. In the propagating formulation, gravity is a full dynamical sector with its own energy, momentum, and angular momentum, and the conservation laws unify cleanly with those of the matter sectors.

## What the Formulation Captures and Does Not

The propagating gravitational sector captures the structural features of relativistic gravity that come from the wave equation: finite propagation speed, retarded potentials, gravitational radiation, and a self-conserved field energy. These features become manifest in the simulation rather than requiring post-hoc reconstruction. A binary system radiates; a collapsing region sends a gravitational signal outward at $c_g$; a moving mass produces a field centered on its retarded position. None of these effects requires special handling; all of them follow from solving the wave equation.

The formulation does not capture the features of relativistic gravity that come from the tensor structure of the metric. Light bending is wrong by a factor of two relative to general relativity, because the spatial part of the metric is missing. The Mercury perihelion shift is absent. Frame dragging is absent. The gravitational radiation from a spherical pulsation is nonzero, where general relativity forbids it (Birkhoff's theorem). Black holes do not have the right structure. The formulation is scalar relativistic gravity, not approximate general relativity. These limitations are structural and cannot be removed without elevating the gravitational sector from a scalar field to a tensor field, which is a much larger commitment.

For the regime in which the framework operates — non-relativistic matter on cosmological scales, weak gravitational fields, and large boxes where gravity propagates across the box in a small fraction of the simulation duration — the missing features of full general relativity are negligible. The propagating formulation is a faithful relativistic-causal upgrade of Newtonian gravity, and it interpolates correctly between the Newtonian limit at large $c_g$ and the Nordström limit at $c_g = c$. The framework is honest about what it is: a classical field theory of gravity that respects causality and conserves field energy, but does not capture the full geometric content of general relativity.

The choice to make this commitment rather than going to full GR is deliberate. Tensor gravity requires a metric field with the structure $g_{μν}(\mathbf{x}, t)$ as the dynamical object, gauge-fixing to handle the diffeomorphism freedom, constraint propagation to maintain the Hamiltonian and momentum constraints, and a substantially more complex integration scheme. The cost is high and the payoff is small for the regime of interest. Scalar relativistic gravity captures the structural features that matter for cosmological structure formation while keeping the algorithmic shape uniform with the rest of the framework: a Lagrangian, a per-cell field state, a spectral propagator, and a clean coupling to matter. This is the right level of model.

## Summary

The gravitational sector is promoted from an instantaneous elliptic constraint to a dynamical scalar field $φ \in \mathcal{G}^0$ with conjugate momentum $π_φ \in \mathcal{G}^0$, governed by an inhomogeneous wave equation with propagation speed $c_g$ and source the total mass density. The Lagrangian has the standard scalar-field structure with kinetic, gradient, and source-coupling terms, and the equation of motion is the d'Alembertian wave equation with the cosmological matter overdensity as the source. Free propagation has linear dispersion $ω_\mathbf{k} = c_g \, |\mathbf{k}|$, mode-by-mode unitary evolution in the canonical plane, and the same algorithmic character as the matter sectors' kinetic step: spectral, exact per mode, one FFT pair per timestep. The Newtonian instantaneous formulation is recovered as $c_g \to \infty$; the Nordström scalar relativistic formulation is recovered at $c_g = c$; intermediate values of $c_g$ are sensible numerical regulators that introduce causality structurally without committing to a specific physical speed. The matter sectors are unchanged in form, with the substitution $Φ \to φ$ throughout. Conservation of energy, momentum, and angular momentum become structurally manifest, with the gravitational field carrying its own contributions to each. The formulation captures finite propagation speed, retarded potentials, and gravitational radiation, but not the tensor features of full general relativity. For non-relativistic matter on cosmological scales, the propagating gravitational sector is the natural relativistic-causal upgrade of the original Newtonian formulation, with the same computational character and a richer set of conserved quantities.
