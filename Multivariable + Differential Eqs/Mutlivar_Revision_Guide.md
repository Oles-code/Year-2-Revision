# MATH50004 Multivariable Calculus — Exam Revision Guide (Q1–Q3)

---

## How to Use This Guide

This guide is built from a reading of the 2022–2025 exam papers and focuses on what comes up **repeatedly and under pressure** in Questions 1–3. The exam structure is predictable:

- **Q1** tests subscript notation, vector-calculus identities, and applications of the integral theorems. It is usually the most technically fiddly question — marks are lost on index manipulation, not ideas.
- **Q2** is almost always a *surface-integral workshop*: you're given a field and a surface, and asked to evaluate the same integral by **two or three different methods** (divergence theorem vs. Stokes vs. projection vs. direct parameterisation). The skill being tested is *choosing the right method* and not getting lost in the algebra.
- **Q3** is calculus of variations: set up a Lagrangian, apply the Euler–Lagrange equation (usually a short form), and frequently add a constraint (Lagrange multiplier). Geodesic problems on curvilinear surfaces are a recurring theme.

Sections are ordered for learning, not exam order. A marker at the top of each section notes which recent papers leaned on it.

---

## 1. Subscript Notation (Index Notation)

**Exam relevance:** *Every* Q1 since at least 2022 has required fluent subscript-notation manipulation. This is where the most marks are lost per minute if you're shaky.

### 1.1 The Three Tools

**Einstein summation convention.** A repeated index in a product is implicitly summed from 1 to 3. So `a_i b_i` means `Σᵢ aᵢbᵢ = a·b`. A repeated (dummy) index can always be renamed: `a_i b_i = a_j b_j`. An index appearing only once (free index) must appear on both sides of an equation.

**Kronecker delta `δᵢⱼ`.** Equals 1 if `i = j`, 0 otherwise. Its action is a *contraction*: `δᵢⱼ aⱼ = aᵢ`. Think of it as a machine that identifies its two indices.

**Permutation (Levi-Civita) symbol `εᵢⱼₖ`.** Equals +1 for even permutations of (1,2,3), −1 for odd, 0 if any two indices match. Totally antisymmetric: `εᵢⱼₖ = εⱼₖᵢ = εₖᵢⱼ = −εⱼᵢₖ`.

**Core expressions:**
- `a·b = aᵢbᵢ`
- `(a × b)ᵢ = εᵢⱼₖ aⱼ bₖ`
- `(∇φ)ᵢ = ∂φ/∂xᵢ` (write `∂ᵢ` for short)
- `∇·A = ∂Aᵢ/∂xᵢ`
- `(∇ × A)ᵢ = εᵢⱼₖ ∂Aₖ/∂xⱼ`
- `(∇²A)ᵢ = ∂²Aᵢ/∂xⱼ∂xⱼ` (Laplacian of a *vector* acts component-wise in Cartesians)

### 1.2 The One Identity You Must Know Cold

$$\varepsilon_{ijk}\,\varepsilon_{ilm} = \delta_{jl}\delta_{km} - \delta_{jm}\delta_{kl}$$

This collapses double cross products into dot products. Note that the *first* index is shared. If you see a different shared index (e.g. `εᵢⱼₖ εₖₗₘ`), cycle the symbols to get the shared index into first position: `εᵢⱼₖ = εₖᵢⱼ`, so

$$\varepsilon_{ijk}\varepsilon_{klm} = \varepsilon_{kij}\varepsilon_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}.$$

This cycling is where most people slip up — write it out slowly the first time.

### 1.3 How to Prove a Vector Identity — The Recipe

1. Take the *i*-th component of the identity you want.
2. Expand every cross product using `ε`, every divergence/curl using `ε` or `∂ᵢ`.
3. Wherever `εε` appears with a shared index, apply the ε–δ identity.
4. Use `δ` to rename indices until every dummy index is paired and every free index matches on both sides.
5. Reassemble into vectors.

**Worked example: `div(A × B) = B · curl A − A · curl B`**

$$\nabla\cdot(A\times B) = \partial_i(\varepsilon_{ijk}A_jB_k) = \varepsilon_{ijk}\,(\partial_i A_j)\,B_k + \varepsilon_{ijk}A_j(\partial_i B_k).$$

In the first term, `εᵢⱼₖ ∂ᵢAⱼ = (curl A)ₖ`, so it reads `Bₖ(curl A)ₖ = B · curl A`. ✓

In the second term, `εᵢⱼₖ ∂ᵢBₖ = −εⱼᵢₖ ∂ᵢBₖ = −(curl B)ⱼ`, so the term is `−Aⱼ(curl B)ⱼ = −A · curl B`. ✓

**Worked example: `2(curl u) × u + ∇(u · u) = 2(u·∇)u`** (2023 Q1a, marked hard)

The *i*-th component is

$$2\varepsilon_{ijk}(\text{curl}\,u)_j u_k + \partial_i(u_j u_j).$$

Expand the curl: `(curl u)ⱼ = εⱼₗₘ ∂ₗuₘ`. Then the first term becomes `2εᵢⱼₖ εⱼₗₘ (∂ₗuₘ)uₖ`. Cycle `εᵢⱼₖ = εⱼₖᵢ` so the shared `j` sits first in both:
`2εⱼₖᵢ εⱼₗₘ (∂ₗuₘ)uₖ = 2(δₖₗ δᵢₘ − δₖₘ δᵢₗ)(∂ₗuₘ)uₖ = 2(∂ₖuᵢ)uₖ − 2(∂ᵢuₖ)uₖ`.

The product-rule term `∂ᵢ(uⱼuⱼ) = 2uⱼ ∂ᵢuⱼ`. After renaming `j → k` this cancels the `−2(∂ᵢuₖ)uₖ` piece.

What's left: `2uₖ ∂ₖuᵢ = 2[(u·∇)u]ᵢ`. ✓

### 1.4 The `curl curl = grad div − Laplacian` Identity

$$\nabla \times (\nabla \times A) = \nabla(\nabla\cdot A) - \nabla^2 A.$$

**Proof (in coordinates).** The *i*-th component of `curl curl A` is

$$\varepsilon_{ijk}\partial_j(\varepsilon_{klm}\partial_l A_m).$$

Cycle the first epsilon: `εᵢⱼₖ = εₖᵢⱼ`, shared index `k` first in both:
`εₖᵢⱼ εₖₗₘ ∂ⱼ∂ₗAₘ = (δᵢₗ δⱼₘ − δᵢₘ δⱼₗ)∂ⱼ∂ₗAₘ = ∂ⱼ∂ᵢAⱼ − ∂ⱼ∂ⱼAᵢ`.

First piece: `∂ᵢ(∂ⱼAⱼ) = [∇(div A)]ᵢ`. Second piece: `−∇²Aᵢ`. ✓

**Why this matters for Q1 2022 and 2025.** If `div A = 0` (A is *solenoidal*), the identity simplifies to `curl curl A = −∇²A`. This immediately tells you: for a divergence-free field, "double curl" and "Laplacian" differ only by a sign. Maxwell-type derivations in 2025 Q1 used exactly this to get a wave equation.

### 1.5 Common Subscript-Notation Slip-Ups (from examiners' reports)

- Confusing the Laplacian of a *scalar* with the Laplacian of a *vector*. The scalar Laplacian is `∂ⱼ∂ⱼφ`. The vector Laplacian has a *free* index: `[∇²A]ᵢ = ∂ⱼ∂ⱼAᵢ`. Write the index.
- Using a dummy index a third time. Each dummy must appear exactly twice in a single product.
- Writing `∇(u·u)` as `∂ᵢ(uⱼuⱼ)` and forgetting to differentiate both copies: it's `2uⱼ∂ᵢuⱼ`, not `uⱼ∂ᵢuⱼ`.

---

## 2. The Operators ∇, div, curl, ∇² in Cartesians

**Exam relevance:** Every Q1 assumes absolute fluency here. Not the place to lose marks.

### 2.1 Definitions (Cartesian)

For scalar field `φ(x,y,z)` and vector field `A = A₁i + A₂j + A₃k`:

$$\nabla\varphi = \left(\frac{\partial\varphi}{\partial x}, \frac{\partial\varphi}{\partial y}, \frac{\partial\varphi}{\partial z}\right), \qquad \nabla\cdot A = \frac{\partial A_1}{\partial x} + \frac{\partial A_2}{\partial y} + \frac{\partial A_3}{\partial z},$$

$$\nabla \times A = \begin{vmatrix} i & j & k \\ \partial_x & \partial_y & \partial_z \\ A_1 & A_2 & A_3 \end{vmatrix}, \qquad \nabla^2\varphi = \frac{\partial^2\varphi}{\partial x^2} + \frac{\partial^2\varphi}{\partial y^2} + \frac{\partial^2\varphi}{\partial z^2}.$$

### 2.2 The Two Always-Zero Identities

$$\text{curl}(\nabla\varphi) = 0, \qquad \text{div}(\nabla\times A) = 0.$$

Geometrically: gradient fields have no "circulation" (the potential function makes loop integrals zero), and curl fields have no "sources" (they just rotate).

In subscript notation:
- `[curl ∇φ]ᵢ = εᵢⱼₖ ∂ⱼ∂ₖφ = 0` because `εᵢⱼₖ` is antisymmetric in `j,k` but `∂ⱼ∂ₖ` is symmetric — the sum is over a tensor times its antisymmetric twin, which is zero.
- Same argument for `div curl A = εᵢⱼₖ ∂ᵢ∂ⱼAₖ = 0`.

### 2.3 Product Rules — Quick Reference

Key formulae you should be able to write without thinking:

- `∇(φψ) = φ∇ψ + ψ∇φ`
- `∇·(φA) = φ (∇·A) + ∇φ · A`
- `∇×(φA) = φ (∇×A) + ∇φ × A`
- `∇·(A×B) = B · curl A − A · curl B`
- `∇×(A×B) = (B·∇)A − B(∇·A) − (A·∇)B + A(∇·B)`
- `∇(A·B) = (B·∇)A + (A·∇)B + B×curl A + A×curl B`

### 2.4 Terminology

- `∇·A = 0` everywhere: A is **solenoidal** ("source-free", "divergence-free", "incompressible" in fluids).
- `∇×A = 0` everywhere: A is **irrotational**. In a simply-connected region this is equivalent to `A = ∇φ` for some potential φ (see §4).

### 2.5 The Special Field `r` and Its Powers

With `r = xi + yj + zk` and `r = |r| = (x²+y²+z²)^{1/2}`:

- `∇r = r/r = r̂` (the unit radial vector)
- `∇(rⁿ) = n rⁿ⁻² r`
- `∇·r = 3`, `∇×r = 0`
- `∇(1/r) = −r/r³`
- `∇²(1/r) = 0` for `r ≠ 0` (appeared in 2025 Q1b)

**Deriving `∇²(1/r)`:** `∇(1/r) = −r/r³`. Then `∇²(1/r) = ∇·(−r/r³) = −(1/r³)(∇·r) − r·∇(1/r³)`. The first piece is `−3/r³`. For the second, `∇(1/r³) = −3r/r⁵`, so `−r·(−3r/r⁵) = 3r²/r⁵ = 3/r³`. Sum: `−3/r³ + 3/r³ = 0`. ✓

---

## 3. Conservative Fields and Potentials

**Exam relevance:** 2022 Q3 relied on recognising extremals are stationary via a gradient argument; 2023 Q1b(ii) explicitly tested "circulation is zero because the field is irrotational". Problem sheets 1, 2, 5 drill this.

### 3.1 Three Equivalent Conditions (in a simply-connected region)

For a continuously differentiable vector field `F`:

1. `F = ∇φ` for some scalar `φ` (a **potential**).
2. `∮_γ F · dr = 0` for every closed curve `γ`.
3. `curl F = 0`.

(1) ⟹ (2) is immediate: along any path from A to B, `∫ F·dr = ∫ ∇φ·dr = φ(B) − φ(A)`, so closed paths give zero.

(1) ⟹ (3): curl of a gradient is always zero.

(2) ⟹ (1): fix A; define `φ(P) = ∫_A^P F·dr`. Path-independence makes φ well-defined. A small displacement `δx i` gives `φ(x+δx,y,z) − φ(x,y,z) = ∫ F₁ dx`, so `∂φ/∂x = F₁`. Similarly for the other components.

(3) ⟹ (2): Stokes' theorem. `∮_γ F·dr = ∫_S curl F · n̂ dS = 0`. (This is why simple-connectedness matters — you need to span γ by a surface sitting inside the region.)

### 3.2 Finding a Potential — The Method

Given `F = (F₁, F₂, F₃)` with `curl F = 0`:

1. Integrate `F₁` in `x`: `φ = ∫F₁ dx + g(y,z)`.
2. Differentiate the result in `y` and match with `F₂`. This pins down `∂g/∂y` and hence `g(y,z) = h(y,z) + k(z)`.
3. Differentiate in `z` and match with `F₃` to pin down `k(z)`.

Sheet 2 Q1 is the textbook example: `v = (2xy+z²)i + (2yz+x²)j + (2xz+y²)k`. Potential: `φ = x²y + y²z + z²x`.

### 3.3 Why "Circulation = 0 When Irrotational" Matters

In 2023 Q1b(ii) you compute a big circulation integral and then are asked to *explain* why it's zero when α = 4 (the irrotational value). The intended one-sentence answer is: "Because `curl u = 0` in a simply-connected region, a potential exists, so `∮ = φ(end) − φ(start) = 0` on a closed curve." Examiners reward this succinctly.

---

## 4. Path and Surface Integrals — Evaluation Technique

### 4.1 Path Integrals — The Parametrisation Recipe

To evaluate `∫_γ F · dr`:

1. Parametrise `γ` by some `t`: `r(t) = (x(t), y(t), z(t))`, with `t` running from `t₀` to `t₁`.
2. Compute `dr/dt = (x', y', z')`.
3. Substitute into `F` so everything becomes a function of `t`.
4. The integral becomes `∫_{t₀}^{t₁} F(t) · r'(t) dt`.

If the field is conservative, use `∫ F·dr = φ(end) − φ(start)` — way faster.

### 4.2 Surface Integrals — Four Methods

When you need `∫_S f dS` or `∫_S A · n̂ dS`, you have a choice:

**(a) Direct parametrisation.** If `S` is given by `r(u,v)`, use
$$dS = \left|\frac{\partial r}{\partial u} \times \frac{\partial r}{\partial v}\right|\, du\, dv, \qquad \hat n\, dS = \pm\frac{\partial r}{\partial u} \times \frac{\partial r}{\partial v}\, du\, dv.$$
Sign of `n̂` is fixed by the orientation condition (e.g., "n̂·k > 0", or "outward").

**(b) Projection theorem.** If `S` isn't parallel to the z-axis at any point and `n̂` is its unit normal, project onto the x–y plane:
$$\int_S f \,dS = \int_\Sigma \frac{f(x,y,z(x,y))}{|\hat n \cdot k|} \, dx\, dy.$$
For a surface `z = φ(x,y)`, `n̂ = (−φₓ, −φᵧ, 1)/√(1+|∇φ|²)`, and `|n̂·k| = 1/√(1+|∇φ|²)`, so
$$\int_S f\, dS = \int_\Sigma f(x,y,\varphi(x,y))\,\sqrt{1+\varphi_x^2+\varphi_y^2}\, dx\, dy.$$
Use when `S` is given explicitly as `z = φ(x,y)`.

**(c) Convert to a volume integral** via the divergence theorem (if `S` is closed or can be closed up).

**(d) Convert to a path integral** via Stokes' theorem (for integrands involving `curl A · n̂`).

### 4.3 How to Pick a Method

| You see... | Try first... |
|---|---|
| Integrand has `curl A · n̂` and the surface is open | Stokes (convert to `∮ A · dr`) |
| Integrand has `A · n̂` and the surface is closed | Divergence theorem |
| Open surface, non-closed-able easily, simple explicit equation `z = f(x,y)` | Projection |
| Surface given parametrically with weird boundary | Direct parametrisation |

**The Q2 exam trick:** examiners frequently ask for the *same* integral by multiple methods. The point is to force you to recognise that whichever way you go, you get the same answer — this is the content of Stokes' theorem / the divergence theorem.

---

## 5. Green's Theorem (2D)

**Exam relevance:** Tested directly in 2024 Q1. Also the ancestor of Stokes and Divergence — understanding it well makes the 3D versions obvious.

### 5.1 Statement

For a plane region `R` bounded by a simple closed curve `C` (traversed so that `R` is on the left, i.e. anticlockwise for a convex `R`), and functions `L(x,y), M(x,y)` with continuous partials on `R`:

$$\oint_C (L\,dx + M\,dy) = \iint_R \left(\frac{\partial M}{\partial x} - \frac{\partial L}{\partial y}\right) dx\, dy.$$

### 5.2 Three Ways to Think of It

**As 2D Stokes.** Write `F = Li + Mj`, so `curl F = (∂M/∂x − ∂L/∂y)k`, and `F·dr = L dx + M dy`. The theorem reads `∮ F·dr = ∫∫ (curl F)·k dA`.

**As 2D divergence theorem.** Put `F = Mi − Lj`. Then `div F = ∂M/∂x − ∂L/∂y`. On `C`, `n̂ ds = dy i − dx j` (outward normal to the anticlockwise curve). Then `F · n̂ ds = M dy − (−L)(−dx) = M dy + L dx`? Careful: compute `(Mi − Lj) · (dy i − dx j) = M dy + L dx`. So Green's reads `∮ F · n̂ ds = ∫∫ div F dA`. ✓

**As an area formula.** Pick `L = −y/2, M = x/2`. Then `∂M/∂x − ∂L/∂y = 1`, so
$$\text{Area}(R) = \tfrac{1}{2}\oint_C (x\,dy - y\,dx).$$
Use this to compute the area of an ellipse, cycloid arch (Sheet 2 Q9), etc. Parametrise the boundary, grind through the integral.

### 5.3 Verification Strategy (2024 Q1)

Exams love "verify Green's theorem for this field on this region". The approach:

1. Compute `∂M/∂x − ∂L/∂y`, integrate over `R`.
2. Parametrise each segment of `C`, compute `∮ L dx + M dy` on each, add.
3. Check the two agree.

**Classic trap on a diamond/quadrilateral (2024 Q1b):** the quadrilateral with vertices `(±a, 0), (0, ±b)` is a **diamond**, NOT a rectangle. The four edges are lines like `x/a + y/b = 1`, not coordinate lines. Examiners flagged this as the most common mistake of 2024.

---

## 6. The Divergence Theorem

**Exam relevance:** Tested every year. The workhorse of Q1 and Q2.

### 6.1 Statement

For a closed surface `S` with outward unit normal `n̂`, enclosing volume `τ`, and a vector field `A` with continuous first derivatives throughout `τ`:

$$\oint_S A \cdot \hat n\, dS = \int_\tau \nabla\cdot A\, d\tau.$$

### 6.2 Proof Idea (Core Case)

Assume `τ` is convex and `S` is split by projection onto `z = 0` into an upper part `S₁: z = f₁(x,y)` and lower part `S₂: z = f₂(x,y)`, both projecting onto the region `Σ`. Consider just the `∂A₃/∂z` piece:

$$\int_\tau \frac{\partial A_3}{\partial z}\,d\tau = \iint_\Sigma \left[\int_{f_2}^{f_1} \frac{\partial A_3}{\partial z} dz\right] dx\, dy = \iint_\Sigma \bigl(A_3(x,y,f_1) - A_3(x,y,f_2)\bigr)\, dx\, dy.$$

On `S₁`, the outward normal has `n̂·k > 0`, and using the projection theorem `dx dy = (n̂·k) dS`. On `S₂`, `n̂·k < 0`, so `dx dy = −(n̂·k) dS`. Substituting:

$$\int_\tau \frac{\partial A_3}{\partial z}\,d\tau = \int_{S_1} A_3 (\hat n \cdot k) dS + \int_{S_2} A_3 (\hat n \cdot k) dS = \oint_S A_3 (\hat n \cdot k) dS.$$

Repeating with projections onto `x = 0` and `y = 0` gives the matching statements for `∂A₁/∂x` and `∂A₂/∂y`. Adding the three gives the theorem. ✓

For non-convex surfaces, cut into convex pieces; the surface integrals over shared internal boundaries cancel.

### 6.3 Gauss' Flux Theorem (a very examinable application)

For a closed surface `S`:

$$\oint_S \frac{\hat n \cdot r}{r^3}\, dS = \begin{cases} 0 & \text{if } O \notin \text{interior of } S \\ 4\pi & \text{if } O \in \text{interior of } S. \end{cases}$$

**Proof.** If `O` is exterior, then `r ≠ 0` throughout `τ`; compute `div(r/r³) = 0` (exercise, similar to `∇²(1/r) = 0`) so divergence theorem gives zero.

If `O` is interior, surround `O` with a small sphere `S_ε` of radius `ε` and apply divergence theorem to the region between `S_ε` and `S`. Since `div(r/r³) = 0` there, we get `∮_{S + S_ε} (n̂·r/r³) dS = 0`. On `S_ε`, the outward normal (relative to the region between, i.e. pointing *towards* `O`) is `−r̂`, and `r = ε`, so
$$\oint_{S_\varepsilon} \frac{(-\hat r)\cdot r}{r^3}\,dS = -\oint_{S_\varepsilon} \frac{\varepsilon}{\varepsilon^3}\,dS = -\frac{1}{\varepsilon^2}\cdot 4\pi\varepsilon^2 = -4\pi.$$
Hence `∮_S (n̂·r/r³) dS = 0 − (−4π) = 4π`. ✓

This appeared as a sub-question in 2025 Q1c.

### 6.4 Useful Corollaries (Sheet 3)

Applying the divergence theorem with clever choices:

- `A = r` gives `∮_S r · n̂ dS = 3V`.
- `A = φc` (c constant vector) gives `∮_S φ n̂ dS = ∫_τ ∇φ dτ`.
- `A = c × G` gives `∮_S n̂ × G dS = ∫_τ curl G dτ`.

These are *not* just trivia — Sheet 3 Q4 derives the last two, and 2025 Q1c(i) asks directly for `∫_V div F dV` when `F = n̂ on S`, which the first corollary gives you (answer: surface area of `S`).

### 6.5 Common Exam Setups

**Spherical region:** `div A` often simplifies to something like `x² + y² + z² = r²`. Then `dV = r² sin θ dr dθ dφ`, and the integral factorises cleanly. (2023 Q1c: `div A = r²`, integrate `∫₀ᵃ r⁴ dr · ∫₀^π sin θ dθ · ∫₀^{2π} dφ = (a⁵/5) · 2 · 2π = 4πa⁵/5`.)

**Closed cylinder with flat/tilted caps** (2024 Q2): break the surface integral into curved side + top + bottom, evaluate each, verify the total matches `∫_V div A dV`.

**Open surface, closed off artificially:** if you're asked for `∫_S A · n̂ dS` on an *open* surface but divergence theorem feels natural, close `S` with a flat cap `S_cap`, apply the theorem, then *subtract* the cap contribution.

---

## 7. Stokes' Theorem

**Exam relevance:** Q2 in 2022, 2023, 2025 hinged on choosing between Stokes and divergence theorem.

### 7.1 Statement

For an open surface `S` bounded by a simple closed curve `γ`, with unit normal `n̂` on `S` chosen so that `n̂` and the orientation of `γ` are related by the right-hand rule, and a vector field `A` with continuous partials:

$$\oint_\gamma A \cdot dr = \int_S (\nabla \times A) \cdot \hat n\, dS.$$

### 7.2 The Most Useful Consequence

**The surface doesn't matter, only the boundary.** If `S₁` and `S₂` are two different open surfaces both bounded by `γ`, then `∫_{S₁} curl A · n̂ dS = ∫_{S₂} curl A · n̂ dS`. Combine them into a closed surface `S₁ − S₂` (with compatible normals); divergence theorem gives zero since `div curl A = 0`.

**Practical use:** if you're asked for `∫_S curl A · n̂ dS` on a horribly curved surface `S`, replace `S` with a simpler surface (usually a disc) that shares the same boundary. The integral is the same.

**Example (Sheet 3 Q7):** `∫_S (∇×A) · n̂ dS` on the upper half of an ellipsoid. Just replace the ellipsoid cap with a flat disc in the `z = 0` plane bounded by the same ellipse. Either surface has the same boundary curve, so either gives the same answer. The flat disc is much easier.

### 7.3 Proof Idea

Write `A = A₁i + A₂j + A₃k` and show `∫_S curl(A₁i) · n̂ dS = ∮_γ A₁ dx` separately for each component. For the first: the surface is `z = f(x,y)`, so `n̂ = (−f_x, −f_y, 1)/√{...}`. Compute
$$\text{curl}(A_1 i) = \left(0,\; \partial_z A_1,\; -\partial_y A_1\right).$$
Then `curl(A₁i) · n̂ dS = [−f_y ∂_z A₁ + (−∂_y A_1)] dx dy` (after the projection theorem cancellation). Setting `Ã₁(x,y) = A₁(x,y,f(x,y))`, the chain rule gives `∂ỹ Ã₁ = ∂_y A_1 + f_y ∂_z A_1`, so the integrand equals `−∂ỹ Ã₁ dx dy`. Green's theorem (2D) then converts this to `∮_C Ã₁ dx`, which on `γ` equals `∮_γ A₁ dx`. Summing over components gives the full theorem.

### 7.4 Verification Strategy (2022 Q2, 2023 Q2)

"Verify Stokes for this field on this surface" workflow:

1. **RHS:** Compute `curl A`. Find `n̂` on the surface. Evaluate `∫_S curl A · n̂ dS` (usually by projection or direct parametrisation).
2. **LHS:** Parametrise the boundary curve `γ` in pieces. For each piece, compute `A · dr` and integrate.
3. Check RHS = LHS.

**Careful with orientation.** After you've fixed `n̂`, the right-hand rule forces the direction of `γ`. Get this wrong and you pick up a sign error — examiners penalise.

**Typical trap:** when `γ` consists of straight-line segments along axes and a curved piece (e.g. the boundary of a half-cylinder has two straight bits and two semicircles), the integrals on some straight segments vanish because `A · dr` happens to be zero there (e.g. if the segment is along the x-axis and `A` has no i-component). Spot these to save time.

---

## 8. Curvilinear Coordinates

**Exam relevance:** 2022, 2023, 2025 all leaned on curvilinear coordinates in Q1 or Q3 (or both). Scale factors for non-standard surfaces are a Q3 staple.

### 8.1 Setup

A curvilinear coordinate system `(u₁, u₂, u₃)` replaces `(x, y, z)`. The position vector `r(u₁, u₂, u₃)` defines the system; the **scale factors** are

$$h_i = \left|\frac{\partial r}{\partial u_i}\right|,$$

and the **unit vectors** are `êᵢ = (1/hᵢ)(∂r/∂uᵢ)`. If the `êᵢ` are mutually perpendicular everywhere, the system is **orthogonal** — all the coordinate systems you'll see in exams (cylindrical, spherical, parabolic, bipolar) are orthogonal.

### 8.2 Standard Systems

**Cylindrical `(r, θ, z)`:** `x = r cos θ, y = r sin θ, z = z`. Scale factors: `h₁ = 1, h₂ = r, h₃ = 1`.

**Spherical `(r, θ, φ)`:** `x = r sin θ cos φ, y = r sin θ sin φ, z = r cos θ`. Scale factors: `h₁ = 1, h₂ = r, h₃ = r sin θ`. (This was asked explicitly in 2023 Q3a.)

### 8.3 Elements of Length, Area, Volume

- `(ds)² = h₁²(du₁)² + h₂²(du₂)² + h₃²(du₃)²` (Pythagoras locally)
- `dτ = h₁ h₂ h₃ du₁ du₂ du₃` (volume of the infinitesimal "curvilinear box")
- `dS` on `uᵢ = const` surface: product of the other two sides, e.g. `dS = h₂ h₃ du₂ du₃` on `u₁ = const`.

**Sphere of radius `a`:** the surface `r = a` has `dS = a² sin θ dθ dφ`. Surface area = `∫₀^{2π} ∫₀^π a² sin θ dθ dφ = 4πa²`.

### 8.4 Gradient, Divergence, Curl, Laplacian

$$\nabla\varphi = \sum_i \frac{1}{h_i}\frac{\partial\varphi}{\partial u_i}\,\hat e_i$$

$$\nabla\cdot A = \frac{1}{h_1 h_2 h_3}\left[\frac{\partial}{\partial u_1}(h_2 h_3 A_1) + \frac{\partial}{\partial u_2}(h_1 h_3 A_2) + \frac{\partial}{\partial u_3}(h_1 h_2 A_3)\right]$$

$$\nabla \times A = \frac{1}{h_1 h_2 h_3}\begin{vmatrix} h_1\hat e_1 & h_2\hat e_2 & h_3\hat e_3 \\ \partial/\partial u_1 & \partial/\partial u_2 & \partial/\partial u_3 \\ h_1 A_1 & h_2 A_2 & h_3 A_3 \end{vmatrix}$$

$$\nabla^2\varphi = \frac{1}{h_1 h_2 h_3}\left[\frac{\partial}{\partial u_1}\!\left(\frac{h_2 h_3}{h_1}\frac{\partial\varphi}{\partial u_1}\right) + \frac{\partial}{\partial u_2}\!\left(\frac{h_1 h_3}{h_2}\frac{\partial\varphi}{\partial u_2}\right) + \frac{\partial}{\partial u_3}\!\left(\frac{h_1 h_2}{h_3}\frac{\partial\varphi}{\partial u_3}\right)\right]$$

**Cylindrical Laplacian:**
$$\nabla^2\varphi = \frac{1}{r}\frac{\partial}{\partial r}\!\left(r\frac{\partial\varphi}{\partial r}\right) + \frac{1}{r^2}\frac{\partial^2\varphi}{\partial\theta^2} + \frac{\partial^2\varphi}{\partial z^2}.$$

**Spherical Laplacian:**
$$\nabla^2\varphi = \frac{1}{r^2}\frac{\partial}{\partial r}\!\left(r^2\frac{\partial\varphi}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\!\left(\sin\theta\frac{\partial\varphi}{\partial\theta}\right) + \frac{1}{r^2\sin^2\theta}\frac{\partial^2\varphi}{\partial\varphi^2}.$$

### 8.5 Laplacian of a Vector in Curvilinear (Watch Out!)

In Cartesians, `∇²A` acts component-wise: `(∇²A)ᵢ = ∇²(Aᵢ)`. **This is NOT true in curvilinear coordinates**, because the unit vectors `êᵢ` themselves vary with position.

The safe route: use `∇²A = ∇(∇·A) − ∇×(∇×A)`. In 2022 Q1 you had to compute `∇²A` in cylindrical polars for `A = v(r,θ)r̂ + w(r,θ)θ̂`, and the only way through was via this identity and the curvilinear curl/div formulas. Students who tried to Laplacian-each-component got the wrong answer.

### 8.6 Parametric Surfaces — Scale Factors in 2D

For a surface given by `r(u,v)` (two parameters, not three), the surface element is

$$dS = \left|\frac{\partial r}{\partial u} \times \frac{\partial r}{\partial v}\right|\, du\, dv.$$

For a path on the surface parameterised as `u = u(t), v = v(t)`,
$$ds^2 = E\,du^2 + 2F\,du\,dv + G\,dv^2, \quad \text{where } E = r_u\cdot r_u,\; F = r_u\cdot r_v,\; G = r_v\cdot r_v.$$

On an orthogonal surface (`F = 0`), this simplifies to `ds² = E du² + G dv²`, and if you write `v = v(u)`, `ds = √{E + G(v')²} du`. This is exactly the setup for geodesic problems in 2023 Q3 and 2025 Q3.

**Worked example (2025 Q3-style):** surface `z = r²` parametrised by `(r, θ)` with `x = r cos θ, y = r sin θ, z = r²`:
- `r_r = (cos θ, sin θ, 2r)`, so `r_r · r_r = 1 + 4r²`.
- `r_θ = (−r sin θ, r cos θ, 0)`, so `r_θ · r_θ = r²`.
- `r_r · r_θ = 0` (orthogonal).
- `|r_r × r_θ| = √{(1+4r²)·r² − 0} = r√{1+4r²}`, so `dS = r√{1+4r²} dr dθ`.
- A path `r = r(θ)` has `ds = √{(1+4r²)(r')² + r²} dθ`.

---

## 9. Calculus of Variations — The Euler-Lagrange Equation

**Exam relevance:** Q3 every year. Worth ≈ 1/3 of the total MVC marks.

### 9.1 The Setup

You're given a "cost functional"
$$I[y] = \int_{x_1}^{x_2} L(x, y, y')\, dx$$
with `y(x₁) = y₁, y(x₂) = y₂` fixed. You want the specific curve `y(x)` making `I` stationary (usually a minimum).

### 9.2 The Vanishing Lemma (required to derive E-L)

**Statement.** If `g` is continuous on `[x₁, x₂]` and `∫_{x₁}^{x₂} g(x) η(x) dx = 0` for *every* smooth function `η` vanishing at the endpoints, then `g ≡ 0`.

**Proof sketch.** If not, `g(x₀) ≠ 0` at some `x₀`. WLOG `g(x₀) > 0`. By continuity, `g > c > 0` on a small neighbourhood `NH` of `x₀`. Now choose a specific `η` — a *bump function* — that's smooth, zero outside `NH`, positive and with positive integral on `NH`. (Sheet 5 Q1 constructs one: `q(x)q(1−x)` where `q(x) = e^{-1/x}` for `x>0`, zero otherwise.) Then `∫ g η > 0`, contradicting the hypothesis. ✓

### 9.3 Derivation of the Euler-Lagrange Equation

Suppose `Y(x)` is the extremising curve. Any nearby curve can be written `y(x) = Y(x) + ε η(x)` with `η(x₁) = η(x₂) = 0`. Then `I(ε) := I[Y + εη]` is a function of `ε`, and it's stationary at `ε = 0`. So `I'(0) = 0`.

Differentiating under the integral,
$$I'(\varepsilon) = \int_{x_1}^{x_2}\!\left(\eta\frac{\partial L}{\partial y} + \eta'\frac{\partial L}{\partial y'}\right) dx.$$
Integrate the second term by parts:
$$\int_{x_1}^{x_2} \eta'\frac{\partial L}{\partial y'} dx = \left[\eta \frac{\partial L}{\partial y'}\right]_{x_1}^{x_2} - \int_{x_1}^{x_2} \eta\frac{d}{dx}\!\left(\frac{\partial L}{\partial y'}\right) dx.$$
The boundary term vanishes since `η(x₁) = η(x₂) = 0`. Setting `ε = 0` and combining,
$$\int_{x_1}^{x_2} \eta(x)\left[\frac{\partial L}{\partial y} - \frac{d}{dx}\!\left(\frac{\partial L}{\partial y'}\right)\right] dx = 0.$$
Since this holds for *every* admissible `η`, the Vanishing Lemma forces
$$\boxed{\;\frac{\partial L}{\partial y} - \frac{d}{dx}\!\left(\frac{\partial L}{\partial y'}\right) = 0.\;}$$

This is the **Euler-Lagrange equation**. Its solutions are called **extremals**.

**A caveat worth knowing.** E-L is necessary but not sufficient — a solution might be a saddle point rather than a min or max. To confirm a minimum strictly you'd check `I''(0) > 0`, which is rarely asked but sometimes appears (2022 Q3a).

### 9.4 Three Short Forms (by which variable is missing)

**Case 1: `L` independent of `y`** (i.e. `∂L/∂y = 0`). Then `d/dx (∂L/∂y') = 0`, so
$$\frac{\partial L}{\partial y'} = \text{const.}$$

**Case 2: `L` independent of `y'`.** Then `d/dx (∂L/∂y') = 0` trivially, so E-L reduces to `∂L/∂y = 0` — an *algebraic* equation, not differential.

**Case 3: `L` independent of `x`** (explicitly — i.e. `∂L/∂x = 0`, though `L` still depends on `x` through `y(x)` and `y'(x)`). In this case
$$L - y'\frac{\partial L}{\partial y'} = \text{const.}$$

**Proof of Case 3.** By the chain rule,
$$\frac{dL}{dx} = \frac{\partial L}{\partial x} + y'\frac{\partial L}{\partial y} + y''\frac{\partial L}{\partial y'} = y'\frac{\partial L}{\partial y} + y''\frac{\partial L}{\partial y'}$$
(using `∂L/∂x = 0`). By E-L, `∂L/∂y = d/dx(∂L/∂y')`. Substituting and combining with the product rule,
$$\frac{dL}{dx} = y'\frac{d}{dx}\!\left(\frac{\partial L}{\partial y'}\right) + y''\frac{\partial L}{\partial y'} = \frac{d}{dx}\!\left(y'\frac{\partial L}{\partial y'}\right).$$
Integrating gives `L − y' ∂L/∂y' = const`. ✓

**Which case applies when.** Check whether `L` contains `x` explicitly, `y` explicitly, or `y'` explicitly. "Explicitly" means if you wrote `L` out, would the symbol literally appear? The implicit `x`-dependence through `y(x)` never counts.

### 9.5 The Three Canonical Examples

**(1) Shortest path between two points:** `L = √{1 + (y')²}`, Case 1 applies (no `y`). `∂L/∂y' = y'/√{1+(y')²} = const`, hence `y' = const`, hence a straight line.

**(2) Brachistochrone (curve of quickest descent):** `L = √{(1+(y')²)/y}`, Case 3 applies (no explicit `x`). Using the short form and the substitution `y = α² sin²θ`:

$$x = \alpha^2(\theta - \tfrac{1}{2}\sin 2\theta), \quad y = \tfrac{1}{2}\alpha^2(1 - \cos 2\theta).$$

This is a **cycloid** — the path traced by a point on a rolling wheel.

**(3) Catenoid (minimal surface of revolution):** minimise `I = ∫ x√{1+(y')²} dx` (a curve `y(x)` rotated around the y-axis). `L = x√{1+(y')²}` is independent of `y`, so Case 1: `∂L/∂y' = xy'/√{1+(y')²} = const = β`. Solving:
$$x = \beta \cosh\!\left(\frac{y - \gamma}{\beta}\right),$$
the **catenary**. The surface of revolution is a **catenoid**. (For some boundary conditions no continuous catenoid exists — the problem bifurcates; relevant in soap film experiments.)

### 9.6 Verifying a Curve is an Extremal (2022 Q3a)

A common exam part: "A family `y = Y(x) + ε η(x)` is given; show `dI/dε = 0` at `ε = 0`." Strategy:

1. Substitute into `I` to get `I(ε)`.
2. Differentiate with respect to `ε` (this is just ordinary calculus of a function of `ε`).
3. Evaluate at `ε = 0`.
4. Show it equals zero — often by integration by parts and orthogonality.

Don't forget the last part: "what is the stationary value of `I`?" Students lost marks in 2022 for skipping that.

---

## 10. Constrained Variational Problems (Lagrange Multipliers)

**Exam relevance:** Appeared in 2022 Q3 (integral constraint), 2024 Q3 (fixed arclength).

### 10.1 Setup

Find `y(x)` stationarising `I = ∫ L dx` subject to an integral constraint `J = ∫ g dx = J₀`. Boundary conditions `y(x₁) = y₁, y(x₂) = y₂` are still imposed.

### 10.2 The Multiplier Rule

There exists a constant `λ` (a **Lagrange multiplier**) such that the constrained extremal satisfies the E-L equation for the new Lagrangian `L + λg`:

$$\frac{\partial}{\partial y}(L + \lambda g) - \frac{d}{dx}\frac{\partial}{\partial y'}(L + \lambda g) = 0.$$

### 10.3 The Derivation Idea

Consider two-parameter perturbations `y = Y + ε₁η₁ + ε₂η₂`. The first-order change in `I` is zero, but we also need the first-order change in `J` to be zero (to respect the constraint). So you need `(η∂L/∂Y − d/dx(η∂L/∂Y'))` orthogonal to `η` *whenever* the corresponding expression for `g` is. Sheet 5 Q7 proves: if `∫ ηf dx = 0` whenever `∫ ηg dx = 0`, then `f = λg` for some constant `λ`. That constant is the multiplier.

### 10.4 Solution Procedure

1. Form `L̃ = L + λg`.
2. Solve the E-L equation for `L̃`. The solution will involve `λ` and two constants of integration.
3. Apply the two boundary conditions to fix the constants. The answer still has `λ` in it.
4. Substitute into the constraint `J = J₀` to solve for `λ`.

**Worked example (2022 Q3c):** minimise `I = ∫₀¹ ((y')² − 2xy + y²) dx` subject to `∫₀¹ y dx = 1`, with `y(0) = 0, y(1) = 1`.

- `L̃ = (y')² − 2xy + y² + λy`.
- E-L: `∂L̃/∂y − d/dx (2y') = (−2x + 2y + λ) − 2y'' = 0`, so `y'' − y = λ/2 − x`.
- General solution: `y = A eˣ + B e⁻ˣ + x − λ/2`.
- Boundary conditions pin down `A, B` in terms of `λ`.
- Finally, `∫₀¹ y dx = 1` fixes `λ`.

### 10.5 Fixed-Arclength Constraint (2024 Q3d,e)

Minimise a functional subject to `∫ √{1 + (y')²} dx = L_fixed`. The constrained E-L becomes "E-L for `L + λ√{1+(y')²}`." For the minimal surface of revolution with fixed arclength, this gives `y = b cosh(x/b) − λ` — a shifted catenary. The shift comes from the `λ`-term acting like a constant subtracted from `L`.

---

## 11. The Euler-Lagrange Equation in Higher Dimensions

**Exam relevance:** Used implicitly in any problem with multiple unknown functions (e.g. parametrised curves `(x(t), y(t))` in 2022 Q3 implicit, 2024 isoperimetric problem). The PDE version showed up in the 2024 Q1 solution notes.

### 11.1 Multiple Dependent Variables

For `I = ∫ L(t, x₁,...,xₙ, x'₁,...,x'ₙ) dt`, you get *n* simultaneous E-L equations, one per `xᵢ`:

$$\frac{\partial L}{\partial x_i} - \frac{d}{dt}\frac{\partial L}{\partial x'_i} = 0, \qquad i = 1, \dots, n.$$

The derivation mirrors the scalar case, but you use independent perturbations `ηᵢ` and apply the vanishing lemma component-wise.

### 11.2 The Isoperimetric Inequality (a beautiful application)

Maximise the area `A = ½∮ (x dy − y dx)` enclosed by a simple closed curve of *fixed perimeter* `l`. Taking `L = ½(xy' − yx') + λ√{(x')² + (y')²}` and applying the two E-L equations gives (after algebra) a circle of radius `λ`. Since the perimeter is `l`, `λ = l/(2π)`, and the enclosed area is `l²/(4π)`. This gives the **isoperimetric inequality**

$$4\pi A \leq l^2,$$

with equality iff the curve is a circle.

### 11.3 Higher-Dimensional (PDE) Euler-Lagrange

For an integral `I = ∫_R L(r, f(r), ∇f(r)) dx dy` (or higher dim) with `f` prescribed on the boundary `C`, the extremal `f` satisfies

$$\frac{\partial L}{\partial f} - \nabla \cdot (\nabla_{\!\nabla f} L) = 0,$$

where `∇_{∇f} L` means the vector whose i-th component is `∂L/∂(∂f/∂xᵢ)`.

**Soap-film / minimal surface equation.** For `L = √{1 + |∇f|²}` (area of the graph `z = f(x,y)`), the E-L equation becomes

$$\nabla \cdot \left(\frac{\nabla f}{\sqrt{1 + |\nabla f|^2}}\right) = 0,$$

which expands to the nonlinear PDE

$$(1 + f_y^2)\,f_{xx} + (1 + f_x^2)\,f_{yy} - 2 f_x f_y f_{xy} = 0.$$

Sheet 5 Q11 verifies the plane and Scherk's surface `f = log(cos x / cos y)` are solutions.

---

## 12. Geodesics — The Recurring Q3 Theme

**Exam relevance:** 2023 Q3 (geodesic on a sphere), 2025 Q3 (geodesic on a paraboloid).

### 12.1 The Template

A geodesic is the shortest path between two points on a surface. The recipe:

1. Parametrise the surface, compute `ds²` in terms of two parameters `(u, v)`.
2. Write `v = v(u)` (or `u = u(v)`) and express arclength as `L = ∫ √{ds²/du²} du`.
3. Identify which short form of E-L applies. Surface geodesics almost always have explicit `u`-independence (because most surfaces have symmetry), so **Case 3** applies: `L − v' ∂L/∂v' = const`.
4. Solve the resulting first-order ODE.

### 12.2 Geodesic on a Sphere (2023 Q3)

On a unit sphere, `ds² = dθ² + sin²θ dφ²`. Writing `φ = φ(θ)` with `φ' = dφ/dθ`, arclength is
$$L = \int \sqrt{1 + (\varphi')^2 \sin^2\theta}\, d\theta.$$
`L` is independent of `φ` (Case 1): `∂L/∂φ' = φ' sin²θ/√{1 + (φ')² sin²θ} = K` (constant). Solving for `φ'`:
$$\frac{d\varphi}{d\theta} = \frac{K\,\csc\theta}{\sqrt{\sin^2\theta - K^2}}.$$
This integrates (not obviously) to `sin(α − φ) = β cot θ`, where `β = K/√{1−K²}`. This is a **great circle** — a plane through the origin intersects the sphere in such a curve, which you can verify from the equation `(sin α)y − (cos α)x = β z`.

### 12.3 Geodesic on `z = r²` (2025 Q3)

`ds² = (1 + 4r²) dr² + r² dθ²`. Writing `r = r(θ)`:
$$L = \int \sqrt{(1 + 4r^2)(r')^2 + r^2}\, d\theta.$$
Case 3 (no explicit `θ`): `L − r' ∂L/∂r' = C`. Grinding through:
$$\frac{d\theta}{dr} = \frac{C \sqrt{1 + 4r^2}}{r\sqrt{r^2 - C^2}}, \quad\text{so}\quad \theta - \theta_0 = C\int_{r_0}^r \frac{\sqrt{1+4r^2}}{r\sqrt{r^2 - C^2}}\,dr.$$
`C` is fixed by the condition that the curve passes through the second point.

### 12.4 Key Takeaway

Once you know `ds²` on a surface, the E-L machinery is automatic. The skill tested is:

- Computing `ds²` from a parametrisation (scale factors, §8.6).
- Recognising the right short form.
- Not panicking when the resulting integral looks ugly — the question usually provides a substitution or asks only for the integral form of the answer.

---

## 13. Key Problem Sheet Questions to Review

The problems below are the ones closest in spirit to exam questions. If short on time, hit these first.

### Sheet 1 — Subscript notation and gradient identities

- **Q7** — The three "big" vector identities via subscript notation: `(a×b)·(a×b) = a²b² − (a·b)²`; `(a×b)·(c×d)` and `(a×b)×(c×d)`. Exactly the skill tested in Q1 every year.
- **Q9** — `curl(φA)`, `div(A×B)`, `A × curl A = ½∇(|A|²) − (A·∇)A`. Prove these in subscript notation; the last one is basically 2023 Q1a backwards.
- **Q8** — Kronecker-delta simplifications. Good warm-up; don't skip.

### Sheet 2 — Path, surface and 2D integrals

- **Q1** — Conservative field, find potential, use it to evaluate a line integral. Bread-and-butter method for Q1 circulation sub-parts.
- **Q3** — Three paths between the same two points, three evaluations of `∫ F · dr`. (Values differ because `F` is *not* conservative — check `curl F`!) Good practice distinguishing conservative vs. non-conservative.
- **Q5, Q6** — Surface integrals by projection. Q5 uses a parabolic cylinder (projection onto `x = 0`), Q6 uses a hemisphere with polar coordinates on the projected disc. Directly mirrors 2023 Q2 and 2025 Q2.
- **Q8** — Verifying Green's theorem on a rectangle. Template for 2024 Q1.
- **Q9** — Area formula from Green's theorem, applied to a cycloid. Technique reused in the 2022 revisit-your-examples section.

### Sheet 3 — Divergence and Stokes theorems

- **Q1** — `∫_V φ div A dV = −∫_V A · ∇φ dV` when `φ` vanishes on the boundary. This is the 3D analogue of integration by parts. Used in the PDE Euler-Lagrange derivation (§11.3).
- **Q2, Q3** — Quick applications: `∮ r · n̂ dS = 3V`; `∮ (n̂·r)/r² dS = ∫ dV/r²`. The first appeared in 2025 Q1c.
- **Q4** — The divergence theorem variants for `∇φ` and `curl A`: `∮ n̂ φ dS = ∫ ∇φ dτ` and `∮ n̂ × A dS = ∫ curl A dτ`. The first underlies the alternative definition of `∇φ`; both reappear in 2024 Q1a(iii).
- **Q6** — Verifying divergence theorem on a closed cone. Practice for the closed-cylinder verification of 2024 Q2.
- **Q7** — Replacing a messy surface (ellipsoid cap) with a flat disc sharing the same boundary, using Stokes. Exactly the trick for 2025 Q2.
- **Q8, Q9, Q10** — Stokes verification on spheres/cones with ring or annular boundaries.

### Sheet 4 — Curvilinear coordinates

- **Q3** — Bipolar coordinates: find scale factors for a non-standard orthogonal system. Drills the `|∂r/∂uᵢ|` machinery.
- **Q4, Q5** — Custom curvilinear system `x = ½(u² − v²), y = uv`: scale factors, unit vectors, divergence and curl computed two ways (curvilinear formula vs. converting to Cartesian). Exceptional practice for 2022 Q1d and 2023 Q3a.
- **Q7–Q9** — Changes of variable in double integrals using Jacobians. Important for projection-theorem integrals when the projected region is awkward.
- **Q10, Q11** — Parametric surface areas (helicoid, torus). Practice for 2025 Q3c-style surface area questions.

### Sheet 5 — Calculus of variations

- **Q1** — Construction of a bump function. Underpins the vanishing lemma proof; exam could ask for a sketch or one property.
- **Q2** — Verifying explicitly that `dI/dε = 0` along a family of perturbed paths. Exact template for 2022 Q3a.
- **Q3** — Basic 1D E-L with given boundary conditions. Answer is `y = sin x`.
- **Q5, Q6** — Geodesics on non-flat surfaces (including the great-circle problem, which is 2023 Q3 essentially pre-solved).
- **Q7** — The multiplier-existence lemma. This is the missing step in the Lagrange-multiplier derivation; worth understanding the argument.
- **Q8, Q9, Q10** — Constrained variational problems. Q8 is the simplest (quadratic `L`, linear constraint). Q9 and Q10 are closer to exam-difficulty (fourth-order-like equations, multiple constraints).
- **Q11** — Expanding the minimal-surface equation and verifying known solutions. Direct practice for the PDE E-L.

---

## 14. Quick-Reference Cheat Sheet

### Key identities (memorise)
- `εᵢⱼₖ εᵢₗₘ = δⱼₗδₖₘ − δⱼₘδₖₗ`
- `curl(curl A) = ∇(div A) − ∇²A`
- `div(A × B) = B · curl A − A · curl B`
- `div(φA) = φ div A + ∇φ · A`
- `curl(∇φ) = 0`, `div(curl A) = 0`

### Integral theorems (know the statements perfectly)
- **Green's:** `∮_C (L dx + M dy) = ∫∫_R (∂M/∂x − ∂L/∂y) dA`
- **Divergence:** `∮_S A · n̂ dS = ∫_τ div A dτ` (S closed, bounds τ)
- **Stokes:** `∮_γ A · dr = ∫_S curl A · n̂ dS` (γ bounds S, right-hand rule)
- **Gauss flux:** `∮_S (n̂·r)/r³ dS = 0` (O outside) or `4π` (O inside)

### Curvilinear (cylindrical, spherical)
- Cylindrical: `h₁=1, h₂=r, h₃=1`, `dV = r dr dθ dz`
- Spherical: `h₁=1, h₂=r, h₃=r sin θ`, `dV = r² sin θ dr dθ dφ`

### Euler-Lagrange (all three forms)
- Full: `∂L/∂y − d/dx(∂L/∂y') = 0`
- No `y`: `∂L/∂y' = const`
- No `x`: `L − y' ∂L/∂y' = const`

### Constrained
- Replace `L` by `L + λg`, apply E-L, then fix `λ` via the constraint.

---

*Good luck.*