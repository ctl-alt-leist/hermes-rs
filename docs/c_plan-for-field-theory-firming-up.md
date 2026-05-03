
CLAUDE: This is actually what we did last. It' good to have an understanding that we worked on this (and we probably learned and made changes to the plan as we went). So, we'll get rid of these notes soon (c field theory planning up), but it's worth having on our mind as we work in case anything comes up related.

# Firming Up the Field Theory

The field-theoretic engine in Hermes works. The dynamics evolves, the visuals look broadly sensible, and the pipeline produces snapshots end-to-end. What it does not yet have is the kind of foundation under which one can stop tweaking dials and start trusting results. Past instabilities — the field failing to concentrate at the right rate, phases looking pathological on the lattice, the slow drift away from the Zel'dovich initialization that worked cleanly in the particle-mesh code — are symptoms of a small set of structural questions that have not yet been answered with full clarity.

This note records what those structural questions are, what the answers should be, and how the work on the new branch should be sequenced so that each layer is verified before the next is built on top of it. The aim is firming up rather than rebuilding: the existing code captures roughly the right physics, and the goal is to put each piece of it onto solid ground rather than start over.

The treatment is organized around seven concerns, which together cover the substantive content of the field-theoretic formulation at scale zero. We begin with the dynamical equation and the conservation laws it implies. We then turn to the lattice, where most of the past pain has lived. We work out the initialization through Zel'dovich displacements in a form that is shared cleanly between the particle-mesh and field representations. We pin down what the FFT actually computes and where its conventions can quietly bite. We collapse the unit system to its minimum. We sketch what the field abstraction in Morphis should look like at this point in the project. And we close with the test suite that should be in place before any cosmological run is reported.

A note on notation. The symbol $Ψ$ is reserved across the framework for the system state — at a single scale $s$ as $Ψ^{(s)}$ or across all scales as the full state of physical reality. The earlier project notes use $\mathbf{Ψ}$ for the Zel'dovich displacement field, which collides with this reservation. In this document the displacement field is written $\mathbf{χ}$, and the same change should propagate through the field-theory notes once the firming-up work begins.

A code-review pass through the existing implementation has run alongside this analysis, and where the two converge on a concrete finding the section in question is annotated with the file and line where the relevant code lives. The findings are summarized in their own section after the test list, but the substantive fixes are folded into the discussion in place. The takeaway from the cross-check: the integrator core (signs, factors, FFT conventions, parameter mapping) is correct, the initialization has two surgical bugs that match one of the predictions in section 3 below, and the resolution constraint of section 2 manifests concretely at the current default parameters.

## 1 · The dynamical equation and its sign conventions

The dark matter wavefunction at scale zero satisfies

$$
ℓ \, \mathbb{1} \, ∂_t α \;=\; -\frac{ℓ^2}{2 \, m_α \, a^2} \, ∇^2 α \;+\; m_α \, Φ \, α
$$

with the gravitational potential determined at each instant by the comoving Poisson equation

$$
∇^2 Φ \;=\; 4 π G \, a^2 \, (ρ_α \;-\; \bar{ρ}_α)
$$

and the mass density given by $ρ_α = m_α \, \bar{α} α$. This is the standard Schrödinger-Madelung structure in cosmological form (Widrow & Kaiser 1993), with $ℓ$ playing the role of $\hbar$, $\mathbb{1}$ playing the role of $i$, and the field's own mass density acting as the source of the gravitational potential it then feels. The factor $a^{-2}$ in the kinetic term comes from $∇^2_{\text{phys}} = a^{-2} \, ∇^2_{\text{comoving}}$ and is essential whenever the box expands — it is correctly present in the Hermes implementation but is missing from the equations as written in project notes 1 and 7, which describe the static-spacetime form. The notes should be brought into alignment with the cosmological form recorded here.

The first thing to internalize about this equation is that it is unitary. The Hamiltonian $H = -ℓ^2 \, / \, (2 \, m_α \, a^2) \, ∇^2 + m_α \, Φ$ is self-adjoint when $Φ$ is real, so the integrated mass

$$
M_α \;=\; \int dV \, m_α \, \bar{α} α
$$

is exactly conserved by the continuum equation. This is not a side property; it is *the* diagnostic that the integrator is doing its job. In numerical practice $M_α$ should drift no faster than floating-point round-off — a growth or decay of $M_α$ at the percent level over a cosmological run is not "approximate conservation" but a sign that the time integrator has lost unitarity. Past episodes of fields exploding or dissipating most likely trace to this.

The integrator that respects unitarity exactly, at machine precision, is split-step Fourier. Strang splitting of one full timestep is

$$
α(t + Δt) \;=\; e^{-\mathbb{1} \, m_α \, Φ \, Δt / (2 \, ℓ)} \;\, e^{\mathbb{1} \, ℓ \, ∇^2 \, Δt / (2 \, m_α \, a^2)} \;\, e^{-\mathbb{1} \, m_α \, Φ \, Δt / (2 \, ℓ)} \, α(t)
$$

with the kinetic factor applied in Fourier space (where $∇^2$ is diagonal) and the two potential half-kicks applied as pointwise multiplications in real space. Each factor is exactly unitary in floating-point arithmetic — the kinetic phase has unit modulus per mode and the potential phase has unit modulus per cell — so the discrete propagator preserves $\| α \|^2$ to round-off regardless of $Δt$. The method is second-order accurate in $Δt$ and is the canonical choice for Schrödinger-Poisson everywhere it appears in the literature.

Two sign conventions silently determine whether the equation says what we want. The first is the sign on the time derivative. With our convention $α = \sqrt{ρ_α \, / \, m_α} \; \exp(\mathbb{1} \, S_α \, / \, ℓ)$ giving bulk velocity $\mathbf{v}_α = ∇ S_α \, / \, m_α$, the choice of $+ \, ℓ \, \mathbb{1} \, ∂_t α$ on the left is what makes gravity attractive: a static overdensity with zero initial phase develops $∂_t \mathbf{v}_α = -∇ Φ$, pointing inward. Flipping the sign on the time derivative inverts this and overdensities expel rather than accrete. The diagnostic for catching such a sign flip is a single timestep on a static Gaussian overdensity with $α$ purely real — $∂_t \, |α|^2$ should be positive at the center.

The second is the sign of the so-called "quantum pressure" that emerges from the Madelung decomposition,

$$
Q \;=\; -\frac{ℓ^2}{2 \, m_α \, a^2} \; \frac{∇^2 \sqrt{ρ_α}}{\sqrt{ρ_α}}
$$

which appears in the Bernoulli equation as $-∇ Q$ and acts as a regularization opposing concentration on small scales. We do not implement $Q$ directly; it falls out of the Laplacian acting on $α$. If $Q$ ever appears as a separately-coded term, the formulation has drifted away from the wavefunction representation and has reintroduced a closure problem we were trying to avoid.

## 2 · The lattice, where the dragons live

The grid imposes three independent constraints. Getting any of them wrong produces precisely the symptoms described in the project — fields that develop grid-scale texture, phase patterns that look pathological, concentration that runs away or dissipates without a clear physical cause.

The first constraint is Nyquist on the density. The cell size $Δx$ must resolve the smallest density features, which at $z = 49$ is trivial because the linear power spectrum is essentially flat at high $k$.

The second constraint, which is the one that most often bites in fuzzy-DM and Schrödinger-Poisson codes, is Nyquist on the *phase*. The lattice can only represent phase differences satisfying

$$
\frac{m_α \, |\mathbf{v}| \, Δx}{ℓ} \;<\; π
$$

between adjacent cells. Any larger and the phase wraps within a cell and aliases as a different velocity. Solving for the maximum representable bulk velocity gives

$$
|\mathbf{v}_{\max}| \;<\; \frac{π \, ℓ}{m_α \, Δx} \;=\; \frac{σ_v}{2}
$$

where the second equality uses the calibration $ℓ \, / \, m_α = σ_v \, Δx \, / \, (2 π)$ from the project notes. With $σ_v = 300 \, \mathrm{km / s}$, the lattice represents bulk velocities up to roughly $150 \, \mathrm{km / s}$ before phase aliasing sets in. Cosmic-web bulk flows comfortably exceed this in collapsed regions; halo internal velocities exceed it routinely. When the phase aliases, the wavefunction develops grid-scale noise that the kinetic propagator then amplifies, and from the outside this looks like the field "exploding" or "developing texture for no reason".

The same constraint, stated at the parameters currently in the Hermes default, is uncomfortably tight. With $ν = ℓ \, / \, m_α = 2000 \, \mathrm{kpc}^2 \, / \, \mathrm{Gyr}$, box length $L = 10 \, \mathrm{Mpc}$, and $N = 64$ giving $Δx \approx 156 \, \mathrm{kpc}$, the phase-Nyquist ceiling is $|\mathbf{v}_{\max}| < π \, ν \, / \, Δx \approx 40 \, \mathrm{kpc / Gyr}$. Equivalently, the de Broglie wavelength at a typical $100 \, \mathrm{kpc / Gyr}$ flow is $λ_{\text{dB}} = 2 π \, ν \, / \, |\mathbf{v}| \approx 125 \, \mathrm{kpc}$, *below* the lattice Nyquist scale $2 \, Δx \approx 312 \, \mathrm{kpc}$. The wave content that the formulation exists to represent is structurally unresolvable on this grid. Increasing $ν$, refining the grid, or shrinking the box are all valid responses; the choice should be made deliberately rather than left implicit.

The cleanest fix is to recalibrate $ℓ \, / \, m_α$ upward. Setting

$$
\frac{ℓ}{m_α} \;=\; \frac{v_{\max} \, Δx}{π}
$$

with $v_{\max}$ chosen as the largest *bulk* velocity expected to appear as a coherent flow gives the right ceiling. For $v_{\max} = 600 \, \mathrm{km / s}$ and $Δx = 0.78 \, \mathrm{Mpc}$, this gives $ℓ \, / \, m_α \approx 0.15 \, \mathrm{Mpc}^2 \, / \, \mathrm{Gyr}$, roughly four times the current calibration. The cost is a slightly larger Jeans length below which structure formation is suppressed, but at this resolution the Jeans length is already at the cell scale and a factor of two does not change the physical reach of the simulation.

The third constraint is the Courant condition on the kinetic step. Split-step is unconditionally stable, but it loses accuracy when the kinetic phase per step at the Nyquist mode exceeds order unity,

$$
\frac{ℓ \, Δt}{m_α \, Δx^2} \;\lesssim\; \frac{1}{π}
$$

For our parameters this allows $Δt$ up to a couple of Gyr from the kinetic step alone, comfortably above the Poisson and electromagnetic Courant bounds, so this is not the binding constraint. Worth a unit test to assert it explicitly, but not a tuning parameter to worry about.

A separate concern is dealiasing of nonlinear products. The Schrödinger equation itself is linear, but the Poisson coupling $Φ \, α$ is a product of two fields, and $Φ$ comes from $|α|^2$. Modes near the Nyquist scale in $α$ produce content beyond Nyquist when squared, which then aliases back into the resolved range. In practice the spectrum of $Φ$ is extremely red because of the $1 \, / \, k^2$ from the inverse Laplacian, and this aliasing is mild enough that production fuzzy-DM codes do not implement the 2/3 truncation rule. If the code ever shows checkerboard patterns in $Φ$, the truncation is the standard fix.

## 3 · Initialization that the particle-mesh and field engines genuinely share

The particle-mesh side already initializes from a Zel'dovich displacement of a Gaussian random density field, and the field side should consume the same random field through the inverse Madelung map. Sharing the random field is structural — it ensures that the two engines start from the *same* physical state, with the only difference being the representation. Any divergence in their evolution is then physics rather than initialization noise.

The shared input is a Gaussian random density contrast $δ_{\text{lin}}(\mathbf{x})$ on the grid, with the linear matter power spectrum scaled to the initialization redshift by the linear growth factor $D_+(z_{\text{init}}) \, / \, D_+(0) \approx 1 \, / \, 50$ at $z = 49$. From this single field, the displacement is

$$
\mathbf{χ}(\mathbf{x}) \;=\; -∇^{-2} ∇ \, δ_{\text{lin}}(\mathbf{x})
$$

with the Fourier-space form

$$
\hat{\mathbf{χ}}_{\mathbf{k}} \;=\; -\frac{\mathbb{1} \, \mathbf{k}}{k^2} \; \hat{δ}_{\mathbf{k}}
$$

and the velocity is $\mathbf{v}(\mathbf{x}) = \dot{D}_+ \, \mathbf{χ}(\mathbf{x})$.

The particle-mesh engine consumes this through Lagrangian sampling — each particle starts at a grid position $\mathbf{q}$, is displaced to $\mathbf{x} = \mathbf{q} + \mathbf{χ}(\mathbf{q})$, and given velocity $\mathbf{v}(\mathbf{q})$. The field engine consumes the same displacement and velocity through the inverse Madelung map.

This is the step where past initializations have gone wrong. The recipe in note 6 — "the linear approximation that the phase is locally a plane wave with the local velocity" — produces, in the obvious reading,

$$
α_{\text{naive}}(\mathbf{x}) \;=\; \sqrt{ρ_α(\mathbf{x}) \, / \, m_α} \; \exp\!\big( \mathbb{1} \, m_α \, \mathbf{v}(\mathbf{x}) \cdot \mathbf{x} \, / \, ℓ \big)
$$

and this is wrong globally. The phase $S(\mathbf{x}) = m_α \, \mathbf{v}(\mathbf{x}) \cdot \mathbf{x}$ has gradient

$$
∇ S \;=\; m_α \, \mathbf{v}(\mathbf{x}) \;+\; m_α \, (\mathbf{x} \cdot ∇) \, \mathbf{v}(\mathbf{x})
$$

which equals $m_α \, \mathbf{v}$ only when the second term vanishes — that is, only when $\mathbf{v}$ is constant. The "local plane wave" approximation pretends each point sees a uniform velocity equal to its local $\mathbf{v}$, but the resulting phase field has a global structure inconsistent with the velocity field it was supposed to encode. When the field is then evolved, the Madelung extraction $\mathbf{v}_α = ∇ \arg(α) \, / \, m_α$ produces a velocity that does not match the velocity that went in.

This is the bug currently sitting in `src/scenes/cosmic_web_field/init.rs:161-162` — the Zel'dovich wavefunction initialization uses precisely this $\mathbf{v}(\mathbf{x}) \cdot \mathbf{x}$ form. The random-density initialization in the same file uses the velocity-potential form correctly at line 321, and the surrounding comment is unusually direct about why: "the kinetic step immediately disperses" the $\mathbf{v} \cdot \mathbf{x}$ phase. The Zel'dovich path needs the same $φ_v$ treatment.

The correct procedure goes through the velocity *potential*. Because Zel'dovich displacement is curl-free (it is the gradient of a scalar), the velocity field is also curl-free, and there exists a scalar $φ_v$ with $\mathbf{v} = ∇ φ_v$. In Fourier space this is

$$
\hat{φ}_{v, \mathbf{k}} \;=\; -\frac{\dot{D}_+}{k^2} \; \hat{δ}_{\mathbf{k}}
$$

The inverse Madelung map then reads

$$
α(\mathbf{x}) \;=\; \sqrt{ρ_α(\mathbf{x}) \, / \, m_α} \; \exp\!\big( \mathbb{1} \, m_α \, φ_v(\mathbf{x}) \, / \, ℓ \big)
$$

where the density is $ρ_α(\mathbf{x}) = \bar{ρ}_α \, \big(1 + δ_{\text{lin}}(\mathbf{x})\big)$. This $α$ has $|α|^2 = ρ_α \, / \, m_α$ exactly, $∇ \arg(α) \, / \, m_α = \mathbf{v}$ exactly, and is consistent at the discrete level with the Madelung continuity equation that governs the dynamics. The phase wraps once or twice across the box at $z = 49$ because peculiar velocities are tiny then, and there is no aliasing to worry about at the start.

The same procedure applies to the baryon wavefunction $β$, with the only change being the mean density $\bar{ρ}_β = (Ω_b \, / \, Ω_m) \, \bar{ρ}_α$. At $z = 49$, baryons trace dark matter on all scales above the baryonic Jeans length, so the same $δ_{\text{lin}}$ and $φ_v$ apply.

### A second initialization bug, in the per-mode amplitude

A second, independent bug in the same initialization code affects the per-mode amplitude of the random Gaussian field. For an unnormalized FFT in a periodic box (the convention used here, sketched in section 4), drawing a Gaussian random field with target power spectrum $P(k)$ requires per-mode complex amplitudes with variance

$$
σ^2(\mathbf{k}) \;=\; P(k) \cdot V_{\text{box}}
$$

in matched units. The current Zel'dovich initialization in `src/scenes/cosmic_web_field/init.rs:68-69` has the inverse, $σ^2 = P(k) \, / \, V_{\text{box}}$, which gives the wrong relative weighting across $k$-modes by a factor of $V_{\text{box}}^2$. The cosmic-web particle-mesh initialization in `src/scenes/cosmic_web_pm/init.rs:167-168` has the correct form. The fix is one character.

### A note on the initialization redshift

A separate but related concern surfaces in the same code. To make linear perturbations visible without integrating from $z = 49$, the current Zel'dovich initialization applies a `perturbation_boost = 50000` factor that pushes the typical density contrast to $δ \sim 0.1$ — well into the nonlinear regime that the linear Zel'dovich approximation does not describe. The boosted initial state corresponds to no physical epoch.

This is a reasonable workaround for early visualization tests, but it should not survive the firming-up work. Two principled alternatives exist. The first is to lower the starting redshift to $z \sim 5\text{–}10$, where naturally-grown perturbations are large enough to evolve visibly without artificial boosting and the linear Zel'dovich approximation remains accurate. The second is to extend the initialization to second-order Lagrangian perturbation theory (2LPT), which captures the leading nonlinear corrections and remains valid through $δ \sim 1$. Either choice keeps the initial state physically meaningful; a $50\,000\times$ boost does not.

For a first pass through the firming-up work the lower-redshift route is simpler and sufficient. 2LPT is worth adding later as a refinement once the linear path is solid.

## 4 · What the FFT actually computes

Spectral methods feel like they should be conceptually simple, but the conventions silently differ between libraries and chasing down a factor of $2 π$ or a sign error costs days when it costs anything at all. Worth being explicit.

On a periodic torus of side $L$ with $N$ uniform cells, the sample points are $x_n = n \, L \, / \, N$ for $n = 0, 1, \ldots, N - 1$. The forward DFT in NumPy convention (which `rustfft` and most Rust crates also follow) is

$$
\hat{α}_j \;=\; \sum_{n = 0}^{N - 1} α_n \, \exp\!\big( -2 π \mathbb{1} \, j \, n \, / \, N \big)
$$

with the inverse carrying the $1 \, / \, N$ normalization. The mode index $j$ corresponds to physical wavenumber

$$
k_j \;=\; \frac{2 π \, j_{\text{signed}}}{L}
$$

where $j_{\text{signed}} = j$ for $j \leq N \, / \, 2 - 1$ and $j_{\text{signed}} = j - N$ for $j \geq N \, / \, 2$, placing negative wavenumbers in the second half of the array. The Nyquist mode is at $|k| = π \, / \, Δx$.

With these conventions the spectral derivative is

$$
\widehat{∂_x α}_j \;=\; \mathbb{1} \, k_j \, \hat{α}_j
$$

the Laplacian multiplies by $-k^2$, and the inverse Laplacian divides by $-k^2$ for $k \neq 0$ with the zero mode set to vanish (the gauge that makes $⟨Φ⟩ = 0$).

The factor of $2 π$ is the conventions trap that I have personally watched eat days of debugging. NumPy's `fftfreq` returns "cycles per unit length", not "radians per unit length"; the wavenumber as used in physics is the latter. Forgetting the factor of $2 π$ gives a Laplacian that is $(2 π)^2 \approx 39.5$ times too small, which silently makes gravity forty times weaker and the Jeans length absurdly long. The diagnostic is to take the spectral Laplacian of $\sin(2 π m \, x \, / \, L)$ for a small mode number $m$ and check that the result is $-(2 π m \, / \, L)^2 \sin(2 π m \, x \, / \, L)$ to floating-point precision.

The Nyquist mode is ambiguous in sign for real-valued fields and a derivative there is not well-defined. For even-order derivatives like the Laplacian this does not matter; for odd-order derivatives like the gradient, the standard fix is to set the Nyquist mode to zero before the multiplication by $\mathbb{1} \, k$. The kinetic step in split-step Schrödinger uses $-k^2$ and so is unaffected, but any diagnostic that computes $\mathbf{v}_α = ∇ \arg(α) \, / \, m_α$ explicitly does need to handle this.

The benefit that justifies the conventions overhead is that spectral derivatives are *exact* on any finite linear combination of grid-resolved Fourier modes, with no $\mathcal{O}(Δx^2)$ truncation error. There is no boundary stencil — the periodic boundary is structural rather than imposed — and there is one nonlocal kernel, the FFT pair, that handles every derivative we care about.

One subtlety specific to R2C transforms (which Hermes uses for field initialization to halve storage relative to a full complex DFT): when generating a Gaussian random field by sampling complex amplitudes in $k$-space and inverting, the result is real-valued only if the amplitudes satisfy the Hermitian symmetry $\hat{α}(-\mathbf{k}) = \hat{α}^*(\mathbf{k})$. The R2C convention enforces this implicitly by storing only the positive-$k_x$ half-space, but the special modes at $k_x = 0$ and $k_x = π \, / \, Δx$ have no negative-$k_x$ partner and must be drawn as real (not complex) to avoid biasing the field with a tiny imaginary residual. The current particle-mesh initialization comments acknowledge this is at noise level for particle positions; for the phase-sensitive field-theoretic initialization, where the Madelung extraction reads the phase to extract velocities, it should be tightened. This is a few lines in the random-amplitude generator, not a structural change.

## 5 · Collapsing the unit system

The Hermes code currently carries $ℓ$ and $m_α$ as separate parameters, but the dynamics depends only on the combination $ν \equiv ℓ \, / \, m_α$. This is the single physical scale of the dark-matter sector — a phase-space volume per time, with dimensions of $\mathrm{Mpc}^2 \, / \, \mathrm{Gyr}$ in our chosen units, equivalently a kinematic quantum diffusivity. Setting $m_α = 1$ in code units (which is just the choice of mass scale) and substituting, the evolution equation reduces to

$$
\mathbb{1} \, ∂_t α \;=\; -\frac{ν}{2 \, a^2} \, ∇^2 α \;+\; \frac{Φ}{ν} \, α
$$

with $|α|^2$ being the density directly, and the Poisson source becoming $∇^2 Φ = 4 π G \, a^2 \, (|α|^2 - \bar{ρ}_α)$. The de Broglie calibration becomes

$$
ν \;=\; \frac{σ_v \, Δx}{2 π}
$$

with no $m_α$ floating around. This collapse is worth doing in code: store and pass $ν$, not $(ℓ, m_α)$ separately. It eliminates a class of bugs where someone changes $m_α$ and forgets to rescale $ℓ$ to keep $ν$ fixed, and it makes the diagnostic that "the dark-matter sector has exactly one tunable parameter" structurally visible rather than merely true.

The comoving Poisson coefficient is the kind of factor that quietly determines structure-formation timing without any obvious failure mode if it is wrong. The code-review pass confirms that the existing implementation uses $∇^2 Φ = 4 π G \, \bar{ρ} \, a^2 \, δ$, which is the standard form when $\bar{ρ}$ is taken as the present-day physical mean density and the Laplacian on the left is in comoving coordinates. The verification is reassuring; the convention should nonetheless be made structurally manifest in the code, with a single named constant defining which density (physical or comoving) is stored on the grid and the Poisson coefficient derived from that choice rather than written by hand. Factor-of-$a$ errors in this part of the code do not crash anything; they just produce growth timing off by powers of $a$, which is exactly the kind of bug that survives until someone compares to linear theory at high precision.

The baryon sector adds one new parameter, the self-coupling $g$ controlling sound speed; everything else is structurally identical to the dark-matter sector with $m_β = m_α$. The electromagnetic sector adds $c_γ$ and $e$ with their own dimensions. None of these introduce hidden combinations the way $ℓ \, / \, m_α$ does.

## 6 · What the field abstraction in Morphis should look like

The current Morphis field type carries a grid and components and supports the basic spectral operations. Now that we have used it for a full pipeline cycle, several patterns have become clear about what the abstraction should look like to support the kind of work the field-theory branch will do.

The field type should carry a type-level grade tag. The grade is not decorative — it constrains at compile time what operations are valid, and it drives the dispatch on operations like the geometric product whose output grade depends on inputs. Adding a scalar field to a bivector field should be a compile error, not a runtime panic.

The field type should be explicit about whether its components are stored in real space or in Fourier space. The cleanest expression is two distinct types — `Field<G, Real>` and `Field<G, Spectral>` — with the FFT pair as the only operations that move between them. This makes every FFT visible in the type system, prevents the accidental application of a real-space pointwise product to spectral coefficients (which is a real bug I have seen happen), and lets the cost of an algorithm be read off from the type signatures.

The operations that feel like they should be primitive on the field type: the pointwise geometric product with grade-correct output, the grade-reverse involution, the FFT pair, the spectral gradient, divergence, and Laplacian, the spectral Poisson solver (a stateless object holding the precomputed inverse-Laplacian Fourier multiplier), and the $L^2$ inner product $\int dV \, \bar{α} β$ that pairs an even-subalgebra field with another and returns a scalar.

The Madelung decomposition deserves to be a first-class operation specifically for the even-subalgebra field type, with an inverse that takes density and velocity *potential* (not velocity directly) and produces the wavefunction. The asymmetric signature is the whole point — forward Madelung gives $(ρ, S)$; inverse Madelung takes $(ρ, φ_v)$ and reconstructs $α$. The inverse Madelung from $(ρ, \mathbf{v})$ is then a higher-level operation that first solves $\mathbf{v} = ∇ φ_v$ for irrotational $\mathbf{v}$ via the spectral relation $φ_v = ∇^{-2} (∇ \cdot \mathbf{v})$, then calls the underlying inverse Madelung on $(ρ, φ_v)$.

The Schrödinger propagator should be a stateful object holding precomputed kinetic propagator coefficients $\exp\!\big( -\mathbb{1} \, ν \, k^2 \, Δt \, / \, (2 \, a^2) \big)$, recomputed only when $Δt$, $ν$, or $a$ changes. Methods `kinetic_step(α, Δt)`, `potential_step(α, Φ, Δt)`, and `strang_step(α, Φ, Δt)` cover the integration vocabulary. The full cosmological step composes a Strang-split Schrödinger update with a Poisson solve at the right point in the splitting.

A possibly contrarian opinion: do not make $(ρ, S)$ a first-class field type that stores state. The dynamics is genuinely cleaner in $α$-form, phase unwrapping in three dimensions is a nightmare we should not depend on for state recovery, and the Madelung pair is best treated as an output for diagnostics rather than a stored representation.

The code-review pass surfaced a concrete set of methods that Hermes currently has to construct manually — decomposing the even-subalgebra field into scalar and pseudoscalar arrays, applying spectral operations component-by-component, and reassembling. These should be primitives on the Morphis field type:

- `EvenField::grad` returns the pair of grade-1 fields (gradient of scalar part, gradient of pseudoscalar part).
- `EvenField::laplacian` applies the spectral Laplacian componentwise.
- `EvenField::integrate` returns $\int dV \, |α|^2$ and $\int dV \, m_α \, |α|^2$ for the norm and mass diagnostics that every conservation test needs.
- `EvenField::kinetic_energy(ν)` returns $(ν \, / \, 2) \, \int dV \, |∇α|^2$, the kinetic-energy functional in the simplified-units form.

Inside Hermes itself, the FFT helper functions `fft_3d` and `ifft_3d` are currently duplicated between `src/dynamics/schrodinger_dynamics.rs` and `src/scenes/cosmic_web_field/init.rs`. They should live in a shared module — `physics::spectral` is the natural home — with no behavioral change. The duplication is the kind of small inconsistency that compounds: the two copies will drift, and when they do the bug will live wherever the next person did not look.

## 7 · Tests that should exist before the run is trusted

The existing test suite has 166 tests, comprehensive for the particle-mesh path, and zero for the field-theoretic sector. This is the single largest blocker between the current state of the code and a run whose results can be trusted. None of the cheap, definitive tests below is currently present in the suite — and most of them have pass-fail criteria that take only a few lines to set up. The list that follows is roughly priority-ordered, with the cheapest and most diagnostic at the top.

**Spectral derivative correctness.** Build a one-dimensional test: $α(x) = \sin(2 π m \, x \, / \, L)$ for $m = 1, 2, \ldots, N \, / \, 2 - 1$, take $∂_x$ and $∇^2$ spectrally, check against the analytic answers $(2 π m \, / \, L) \cos(\ldots)$ and $-(2 π m \, / \, L)^2 \sin(\ldots)$. Floating-point precision is the green light. Also assert that the Nyquist mode is exactly zero in the real-space output of an odd-order spectral derivative applied to a real field.

**Laplacian self-adjointness and sign.** For two random complex fields $α, β$ on the grid, check $⟨α \,|\, ∇^2 β⟩ = ⟨∇^2 α \,|\, β⟩$ and $-⟨α \,|\, ∇^2 α⟩ \geq 0$. Costs almost nothing, catches sign errors instantly.

**Poisson round-trip.** Build a known density of low-mode sinusoids, compute $Φ = ∇^{-2} (ρ - ⟨ρ⟩)$, then $∇^2 Φ$, and compare to $ρ - ⟨ρ⟩$. Floating-point precision.

**Schrödinger plane wave.** Initialize $α(\mathbf{x}, 0) = \exp(\mathbb{1} \, \mathbf{k}_0 \cdot \mathbf{x})$, evolve under no gravity, check $α(\mathbf{x}, t) = \exp(\mathbb{1} \, \mathbf{k}_0 \cdot \mathbf{x} - \mathbb{1} \, ω \, t)$ with $ω = ν \, k_0^2 \, / \, 2$. Tests sign and the factor of $1 \, / \, 2$ simultaneously.

**Schrödinger Gaussian wave packet.** Initialize a Gaussian envelope with a uniform momentum phase, evolve under no gravity, verify the centroid moves at $\mathbf{v}_g = ν \, \mathbf{k}_0$ and the width grows as $σ(t) = σ_0 \, \sqrt{1 + (ν \, t \, / \, σ_0^2)^2}$. Tests the kinetic step end-to-end.

**Norm and energy conservation.** Run any of the above for many timesteps, assert $\| α(t) \|^2 - \| α(0) \|^2$ stays at floating-point round-off (split-step is exactly unitary), and assert the energy $⟨H⟩$ drifts only as $\mathcal{O}(Δt^2)$ with Strang splitting.

**Madelung round-trip.** Take a known $(ρ, \mathbf{v})$ with $\mathbf{v}$ curl-free, build $α$ via inverse Madelung, decompose, recover $(ρ, \mathbf{v})$ to floating-point precision modulo a global phase constant.

**Zel'dovich consistency.** Generate $δ_{\text{lin}}$ from a known power spectrum, build $α$ via the velocity-potential route described in section 3, check that $m_α \, |α|^2 = \bar{ρ}_α \, (1 + δ_{\text{lin}})$ and $ν \, ∇ \arg(α) = \mathbf{v}_{\text{lin}} = \dot{D}_+ \, \mathbf{χ}$ to floating-point precision.

**Linear growth test.** Start with small $δ_{\text{lin}}$, evolve for a small $Δt$, compare growth of the modes of $|α|^2$ to the linear-theory growth factor $D_+(t) \, / \, D_+(t_0)$. This is the integration-level sanity check that gravity is sourced and applied with the right magnitude.

**Particle-mesh / field cross-check.** Same $δ_{\text{lin}}$, run particle-mesh and field side-by-side, compare power spectra and density fields at small $t$. They should be statistically indistinguishable in the linear regime — any divergence at large scales is a calibration error, not physics.

**Spherical symmetry.** A spherically symmetric initial overdensity should remain spherically symmetric under self-gravity, up to grid discretization at the cell scale. Catches asymmetries introduced by FFT-shift or wavenumber-array mistakes.

**Schrödinger-Newton soliton.** This is the integration-level test for the full gravity-coupled wavefunction sector. There is a known stationary solution to $\mathbb{1} \, ℓ \, ∂_t α = -ℓ^2 \, / \, (2 \, m_α \, a^2) \, ∇^2 α + m_α \, Φ \, α$ with $∇^2 Φ = 4 π G \, m_α \, |α|^2$ — the boson-star or Schrödinger-Newton ground state. Initialize close to the published ground-state profile, evolve, verify the profile is stationary and the central density does not drift more than a percent over many dynamical times. If this passes, the gravity coupling is solid.

## Findings from the code-review pass

The summary below pulls together what the code-review pass through Hermes and Morphis verified, what it caught, and what was already in better shape than expected. The findings are not new conclusions; they are concrete grounding for the analysis in the preceding sections.

What is verified correct: the signs and factors in the split-step integrator match the cosmological Schrödinger-Poisson equation. Kinetic phase $-ℓ \, k^2 \, Δt \, / \, (2 \, m_α \, a^2)$, potential phase $-m_α \, Φ \, Δt \, / \, ℓ$, Poisson source $4 π G \, \bar{ρ} \, a^2 \, δ$ with the inverse Laplacian dividing by $-k^2$ and the zero mode projected out, and right-multiplication by $\cos θ + \mathbb{1} \, \sin θ$ for even-subalgebra phase rotation are all correct. The configuration-to-runtime mapping $ℓ = (ℓ \, / \, m) \cdot m$ at `src/scenes/cosmic_web_field/mod.rs:49` is correct. FFT conventions are internally consistent across both repos: Hermes uses R2C with no normalization on forward and $1 \, / \, N^D$ on inverse; Morphis uses C2C with the same normalization split. The discrete sin² Green's function in the particle-mesh Poisson solver and the continuous $-k^2$ kernel in the spectral solver are the right division of responsibilities — the particle-mesh chain matches its CIC deposition, the spectral chain is exact.

What needs surgical fixes: two bugs in `src/scenes/cosmic_web_field/init.rs`. The Zel'dovich wavefunction phase at lines 161–162 uses $\mathbf{v}(\mathbf{x}) \cdot \mathbf{x}$ where it should use the velocity potential $φ_v(\mathbf{x})$ — the failure mode predicted in section 3, and one the code's own sibling routine warns against in a comment. The per-mode amplitude at lines 68–69 has $σ = \sqrt{P \, / \, V_{\text{box}}}$ where it should be $σ = \sqrt{P \cdot V_{\text{box}}}$ — the cosmic-web particle-mesh init at `src/scenes/cosmic_web_pm/init.rs:167-168` has the right form. Both are local, surgical, and testable.

What is structurally tight at current parameters: the de Broglie wavelength at typical flow speeds is below the lattice Nyquist scale, as worked out in section 2. Either $ν$ goes up, the grid refines, or the box shrinks — but the current configuration cannot resolve the wave content the formulation is built around, and this is the most likely diagnosis behind the field developing grid-scale noise.

What is a workaround that should not survive: the $50\,000\times$ perturbation boost at high redshift, discussed in section 3, makes the linear initial state visible at the cost of breaking its physical meaning. A lower starting redshift or 2LPT replaces it cleanly.

What needs a discipline pass: the equations as written in project notes 1 and 7 omit the $a^{-2}$ in the kinetic term that the cosmological form requires and that the Hermes implementation already includes. The notes should be updated to match. The R2C random-amplitude generator in the field initialization should enforce strict Hermitian symmetry at the special modes. The duplicated FFT helpers in Hermes should consolidate into `physics::spectral`. None of these blocks the firming-up work; they are housekeeping that becomes harder to do later.

What is the largest gap: the test suite has 166 tests for the particle-mesh path and zero for the field-theoretic sector. Section 7's list is the path to closing it.

## Sequencing the branch

The point of having tests at every level is that the work can be sequenced so each layer is verified before the next is built on top. The sequence I would follow, in order:

Zeroth, before anything structural, fix the two surgical bugs in `src/scenes/cosmic_web_field/init.rs` — the velocity-potential phase form and the inverted per-mode amplitude. These are local, the corrected forms are already documented in sections 3 and 4 of this note, and the corresponding tests in section 7 (Madelung round-trip and Zel'dovich consistency) will exercise the fix. This unblocks the rest of the work by making the existing field-theory path produce sensible initial states again.

First, get the spectral kit in Morphis to floating-point precision on the cheap tests in section 7 — derivative correctness, Laplacian self-adjointness, Poisson round-trip. Do not touch dynamics until these pass.

Second, build the split-step Schrödinger propagator in Hermes consuming Morphis fields, and verify with the plane-wave and Gaussian-packet tests. No gravity yet.

Third, add the Poisson coupling and verify with the Schrödinger-Newton soliton.

Fourth, replace the existing initialization with the velocity-potential Zel'dovich route, share $δ_{\text{lin}}$ between the particle-mesh and field engines, drop the $50\,000\times$ perturbation boost in favor of a lower starting redshift, and verify the linear-growth and cross-check tests.

Only then, fifth, run the full cosmological pipeline and look at the cosmic web.

Each step has a definitive pass-fail criterion. The pattern that produces past pain is running the full pipeline before each layer is verified, then trying to debug "the field doesn't behave" with three or four layers of physics simultaneously active and no idea which is at fault. The discipline of verifying every layer before adding the next is what turns the simulation from a thing that works most of the time into a thing whose results can be reported and trusted.
