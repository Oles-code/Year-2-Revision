# Analysis II — Term 2: Complex Analysis

## Comprehensive Exam Revision Guide

**MATH50001 — Imperial College London**

---

**How to use this guide.** This document distils the key examinable content from the Term 2 (Complex Analysis) lecture notes and problem sheets. It is structured around topics that appear most frequently and in the most challenging contexts across the 2022–2025 final exams. Proofs are included where they are short enough to plausibly appear on an exam. At the end, specific problem sheet questions are referenced that practise the most exam-relevant skills.

**Exam structure.** The paper has 6 questions of equal weight (≈20 marks each). Questions 4–6 cover Term 2 (this material). Typically: Q4 tests complex differentiability and Cauchy–Riemann equations, harmonic functions, and complex logarithm/powers; Q5 tests Laurent series, singularities, and residues; Q6 tests applications — computing real integrals via residues, Rouché's theorem, maximum modulus principle, and Möbius transformations.

---

## Part 1: The Complex Plane — Foundations

### 1.1 Complex Numbers

> **Definition.** A *complex number* is z = x + iy where x, y ∈ ℝ and i² = −1. The *real part* is Re(z) = x and the *imaginary part* is Im(z) = y.

> **Definition.** The *complex conjugate* of z = x + iy is z̄ = x − iy. The *modulus* is |z| = √(x² + y²) = √(zz̄).

**Key identities (used constantly):**

- |z₁z₂| = |z₁||z₂|, |z₁/z₂| = |z₁|/|z₂|
- Re(z) = (z + z̄)/2, Im(z) = (z − z̄)/(2i)
- |z|² = zz̄
- Triangle inequality: |z₁ + z₂| ≤ |z₁| + |z₂|
- Reverse triangle inequality: ||z₁| − |z₂|| ≤ |z₁ − z₂|

### 1.2 Polar Form and De Moivre's Formula

For z = x + iy ≠ 0, let r = |z| and θ = arg(z). Then z = r(cos θ + i sin θ) = re^{iθ}.

> **Definition.** The *principal value of the argument*, denoted Arg(z), is the unique θ satisfying −π < θ ≤ π.

**De Moivre's formula:** For n ∈ ℕ, z^n = r^n(cos nθ + i sin nθ).

⚠️ **Exam warning:** Arg(z₁) + Arg(z₂) ≠ Arg(z₁z₂) in general — the sum may fall outside (−π, π]. This subtlety is tested repeatedly in the context of the logarithm.

### 1.3 Topology of ℂ

Since ℂ ≅ ℝ², we can use the metric d(z₁, z₂) = |z₁ − z₂|. Key topological notions carry over directly from Term 1.

> **Definition.** The *open disc* of radius r centred at z₀ is D_r(z₀) = {z ∈ ℂ : |z − z₀| < r}.

> **Definition.** Ω ⊆ ℂ is *open* if every point of Ω is an interior point. It is *compact* if and only if it is closed and bounded (Heine–Borel applies here since ℂ ≅ ℝ²).

> **Definition.** An open set Ω ⊆ ℂ is *connected* if any two points in Ω can be joined by a curve lying entirely in Ω.

---

## Part 2: Holomorphic Functions and the Cauchy–Riemann Equations

### 2.1 Complex Differentiability

> **Definition.** Let Ω₁, Ω₂ ⊆ ℂ be open sets and f : Ω₁ → Ω₂. We say f is *differentiable* (or *holomorphic*) at z₀ ∈ Ω₁ if the limit
>
> f'(z₀) = lim_{h→0} [f(z₀ + h) − f(z₀)] / h
>
> exists, where h ∈ ℂ and z₀ + h ∈ Ω₁.

**Crucial remark.** Unlike real differentiation, the increment h may approach 0 from *any direction* in the complex plane. This makes complex differentiability a much stronger condition than real differentiability.

> **Definition.** f is *holomorphic on an open set* Ω if it is holomorphic at every point of Ω. f is *holomorphic on a closed set* if it is holomorphic on an open set containing it. f is *entire* if it is holomorphic on all of ℂ.

**Standard examples:**

- f(z) = z^n is entire with f'(z) = nz^{n−1}
- Any polynomial p(z) = a₀ + a₁z + ⋯ + aₙzⁿ is entire
- f(z) = 1/z is holomorphic on ℂ \ {0}
- **f(z) = z̄ is NOT holomorphic** — the limit h̄/h has no limit as h → 0 (it depends on direction)

> **Proposition 2.1.** f is holomorphic at z₀ if and only if there exists a ∈ ℂ such that f(z₀ + h) − f(z₀) − ah = hψ(h), where ψ(h) → 0 as h → 0. Moreover, a = f'(z₀).

> **Corollary.** Every holomorphic function is continuous.

**Algebraic properties:** If f, g are holomorphic on Ω, then so are f + g, fg, f/g (where g ≠ 0), and g ∘ f (chain rule applies):

- (f + g)' = f' + g'
- (fg)' = f'g + fg'
- (f/g)' = (f'g − fg')/g²
- (g ∘ f)'(z) = g'(f(z)) · f'(z)

### 2.2 The Cauchy–Riemann Equations ⭐

This is the single most important topic in the course — tested every year without exception.

Write f(x + iy) = u(x, y) + iv(x, y), where u and v are real-valued functions of two real variables.

> **Theorem (Cauchy–Riemann Equations).** If f is holomorphic at z₀, then the partial derivatives of u and v exist at z₀ and satisfy:
>
> ∂u/∂x = ∂v/∂y,    ∂u/∂y = −∂v/∂x.

> **Theorem (Converse / Sufficient condition).** Suppose f = u + iv is defined on an open set Ω ⊆ ℂ. If u and v are *continuously differentiable* on Ω and satisfy the Cauchy–Riemann equations on Ω, then f is holomorphic on Ω and
>
> f'(z) = ∂f/∂z(z) = uₓ − iuᵧ = uₓ + ivₓ.

**Proof of the sufficient condition (examinable — short).**

Let h = h₁ + ih₂. Since u is differentiable as a real function of two variables:

u(x + h₁, y + h₂) − u(x, y) = uₓh₁ + uᵧh₂ + |h|φ₁(h)

and similarly for v. Combining:

f(z + h) − f(z) = (uₓ + ivₓ)h₁ + (uᵧ + ivᵧ)h₂ + |h|φ(h).

Using the C–R equations vₓ = −uᵧ and uₓ = vᵧ:

f(z + h) − f(z) = (uₓ − iuᵧ)(h₁ + ih₂) + |h|φ(h) = (uₓ − iuᵧ)h + |h|φ(h).

Dividing by h and letting h → 0 gives f'(z) = uₓ − iuᵧ. ∎

**The Wirtinger operators.** Define:

∂/∂z = (1/2)(∂/∂x + (1/i)∂/∂y),    ∂/∂z̄ = (1/2)(∂/∂x − (1/i)∂/∂y).

Then f is holomorphic at z₀ if and only if (∂f/∂z̄)(z₀) = 0, and in that case f'(z₀) = (∂f/∂z)(z₀).

### 2.3 Cauchy–Riemann Equations in Polar Coordinates ⭐

> **Theorem.** If f = u + iv is holomorphic, then in polar coordinates (r, θ):
>
> ∂u/∂r = (1/r)(∂v/∂θ),    ∂v/∂r = −(1/r)(∂u/∂θ).

**How to verify (exam technique):** Given f(z), write z = re^{iθ}, compute u(r,θ) and v(r,θ), then check these two equations.

**Example (exam-style).** For f(z) = z² = r²e^{2iθ}:

u(r,θ) = r² cos 2θ, v(r,θ) = r² sin 2θ.

Then uᵣ = 2r cos 2θ, (1/r)vθ = (1/r)(2r² cos 2θ) = 2r cos 2θ. ✓

And vᵣ = 2r sin 2θ, −(1/r)uθ = −(1/r)(−2r² sin 2θ) = 2r sin 2θ. ✓

---

## Part 3: Elementary Functions

### 3.1 The Complex Exponential

> **Definition.** For z = x + iy, e^z = eˣ cos y + ieˣ sin y.

**Key properties:**

1. Extends the real exponential: e^{x+0i} = eˣ
2. Multiplicative: e^{z₁+z₂} = e^{z₁}e^{z₂}
3. 2πi-periodic: e^{z+2πi} = e^z
4. |e^z| = eˣ (depends only on the real part)
5. arg(e^z) = y + 2πk, k ∈ ℤ
6. Entire, and (d/dz)e^z = e^z
7. For any holomorphic g: (d/dz)e^{g(z)} = e^{g(z)}g'(z)

⚠️ **Exam warning (2025 markers):** cos z is NOT bounded for complex z — because cos z = (e^{iz} + e^{−iz})/2 and the real part of iz can be arbitrarily large.

### 3.2 Complex Trigonometric Functions

> **Definition.** For z ∈ ℂ:
>
> cos z = (e^{iz} + e^{−iz})/2,    sin z = (e^{iz} − e^{−iz})/(2i).

**Properties:**

- sin z and cos z are entire functions
- (d/dz) sin z = cos z, (d/dz) cos z = −sin z
- sin²z + cos²z = 1
- **The zeros of sin z are exactly z = nπ, n ∈ ℤ (all real)**
- **The zeros of cos z are exactly z = (n + 1/2)π, n ∈ ℤ (all real)**

**Proof that sin z has only real zeros (exam-style):** Write sin z = sin(x + iy) = sin x cosh y + i cos x sinh y. For sin z = 0 we need both sin x cosh y = 0 and cos x sinh y = 0. Since cosh y ≥ 1 > 0, the first equation gives sin x = 0, so x = nπ. Then cos(nπ) = (−1)ⁿ ≠ 0, so the second equation gives sinh y = 0, hence y = 0.

### 3.3 The Complex Logarithm ⭐

> **Definition (multi-valued).** For z = re^{iθ} ≠ 0:
>
> log z = {ln r + i(θ + 2πk) : k ∈ ℤ}.

Each value is called a *branch* of the logarithm, and each branch satisfies e^{log z} = z.

> **Definition (principal value).** Log z = ln|z| + i Arg z, where −π < Arg z ≤ π.

**Properties:**

- log(z₁z₂) = log z₁ + log z₂ (as sets, i.e. multi-valued)
- **Log(z₁z₂) ≠ Log(z₁) + Log(z₂) in general** — this fails when the arguments sum to something outside (−π, π]
- Log z is holomorphic on ℂ \ (−∞, 0] (the *branch cut*)
- (d/dz) Log z = 1/z on ℂ \ (−∞, 0]

**Common exam computations:**

- Log(−1) = iπ
- Log(2i) = ln 2 + iπ/2
- Log(1 − i) = ln√2 − iπ/4

⚠️ **Exam warning (2024 markers):** Many students were "not familiar with the notion of the principal value of the log-function." Make sure you know the definition and the branch cut.

### 3.4 Complex Powers ⭐

> **Definition (single-valued).** For z ≠ 0 and α ∈ ℂ:
>
> z^α = e^{α Log z} = e^{α(ln|z| + i Arg z)}.

As a multi-valued function: z^α = e^{α log z} takes (in general) infinitely many values.

**Key exam example:** i^α as multi-valued:

i^α = e^{α log i} = e^{α(iπ/2 + 2πki)} = e^{αi(π/2 + 2πk)}, k ∈ ℤ.

In particular, i^i = e^{i·Log(i)} = e^{i·iπ/2} = e^{−π/2} (principal value), and all values |i^α| are equal when α is purely imaginary.

---

## Part 4: Integration in the Complex Plane

### 4.1 Contour Integrals

> **Definition.** A *parametrised curve* γ is a function z : [a,b] → ℂ. It is *smooth* if z'(t) exists, is continuous, and never zero. It is *piecewise smooth* if it is smooth on finitely many subintervals.

> **Definition.** For f continuous on a smooth curve γ parametrised by z(t):
>
> ∫_γ f(z) dz = ∫_a^b f(z(t)) z'(t) dt.

> **Definition.** The *length* of γ is length(γ) = ∫_a^b |z'(t)| dt.

**Properties:**

- **Linearity:** ∫_γ (αf + βg) dz = α∫_γ f dz + β∫_γ g dz
- **Reverse orientation:** ∫_{γ⁻} f(z) dz = −∫_γ f(z) dz
- **ML-inequality:** |∫_γ f(z) dz| ≤ sup_{z∈γ} |f(z)| · length(γ)

The ML-inequality is used constantly in estimates — make sure you are comfortable applying it.

**Standard parametrisation:** The circle C_r(z₀) centred at z₀ of radius r has parametrisation z(t) = z₀ + re^{it}, t ∈ [0, 2π], with z'(t) = ire^{it} and length = 2πr.

**Fundamental computation:** ∮_{|z−z₀|=r} (z − z₀)^n dz = {2πi if n = −1; 0 otherwise} for any integer n.

### 4.2 Primitives

> **Definition.** A *primitive* of f on Ω is a holomorphic function F such that F'(z) = f(z) for all z ∈ Ω.

> **Theorem (Fundamental Theorem of Contour Integration).** If f is continuous on Ω and F is a primitive of f, and γ is a piecewise-smooth curve in Ω from w₁ to w₂, then ∫_γ f(z) dz = F(w₂) − F(w₁).

> **Corollary.** If f has a primitive in Ω and γ is a closed curve in Ω, then ∮_γ f(z) dz = 0.

**Key example:** f(z) = 1/z has NO primitive in ℂ \ {0}, because ∮_{|z|=1} (1/z) dz = 2πi ≠ 0.

### 4.3 Cauchy–Goursat Theorem ⭐

> **Theorem (Cauchy–Goursat for a triangle).** Let Ω ⊆ ℂ be open, T a triangle whose interior lies in Ω. If f is holomorphic in Ω, then ∮_T f(z) dz = 0.

The proof uses a nested triangle subdivision argument — it is long (≈1 page) but has appeared in exam-adjacent contexts.

> **Corollary (for a disc).** If f is holomorphic in a disc, then ∮_γ f(z) dz = 0 for any closed curve γ in that disc.

> **Theorem.** A holomorphic function in an open disc has a primitive in that disc.

### 4.4 Simply Connected Domains and the General Cauchy Theorem

> **Definition.** An open set Ω ⊆ ℂ is *simply connected* if any two curves in Ω with the same endpoints are homotopic in Ω.

**Examples:** Any disc is simply connected. Any convex open set is simply connected. ℂ \ (−∞, 0] is simply connected. The punctured plane ℂ \ {0} is **not** simply connected.

> **Theorem (Cauchy–Goursat for simply connected domains).** If f is holomorphic in a simply connected open set Ω, then ∮_γ f(z) dz = 0 for any closed, piecewise-smooth curve γ ⊆ Ω.

> **Theorem.** Any holomorphic function in a simply connected domain has a primitive.

> **Theorem (Deformation Theorem).** If γ₁ and γ₂ are simple, closed, piecewise-smooth curves with γ₂ inside γ₁, and f is holomorphic in a domain containing the region between them, then ∮_{γ₁} f(z) dz = ∮_{γ₂} f(z) dz.

This is the workhorse for computing integrals: deform the contour to a small circle around each singularity.

---

## Part 5: Cauchy's Integral Formula and Applications

### 5.1 Cauchy's Integral Formula ⭐

> **Theorem.** Let f be holomorphic inside and on a simple, closed, piecewise-smooth curve γ. For any z₀ interior to γ:
>
> f(z₀) = (1/2πi) ∮_γ f(z)/(z − z₀) dz.

**Proof (examinable — short).** By the deformation theorem, replace γ by a small circle γ_r = {z : |z − z₀| = r} inside γ. Split f(z) = f(z₀) + (f(z) − f(z₀)):

(1/2πi) ∮_{γ_r} f(z)/(z − z₀) dz = f(z₀) + (1/2πi) ∮_{γ_r} [f(z) − f(z₀)]/(z − z₀) dz.

The second integral vanishes by the ML-inequality: |∮_{γ_r} [f(z)−f(z₀)]/(z−z₀) dz| ≤ sup|f(z)−f(z₀)|/r · 2πr → 0 as r → 0 by continuity. ∎

### 5.2 Generalised Cauchy Integral Formula ⭐

> **Theorem.** Under the same hypotheses, f has infinitely many complex derivatives in Ω, and for any z inside γ and any n ≥ 0:
>
> f^{(n)}(z) = n!/(2πi) ∮_γ f(ξ)/(ξ − z)^{n+1} dξ.

**Exam technique:** To evaluate ∮_γ g(z)/(z − z₀)^{n+1} dz where g is holomorphic inside γ, use:

∮_γ g(z)/(z − z₀)^{n+1} dz = 2πi · g^{(n)}(z₀)/n!

### 5.3 Key Consequences

> **Corollary (Infinite differentiability).** If f is holomorphic in Ω, then f', f'', f''', … all exist and are holomorphic in Ω.

> **Corollary (Liouville's Theorem).** If f is entire and bounded, then f is constant.

**Proof (examinable — very short).** Let |f(z)| ≤ M for all z. Fix z₀ and integrate over a circle of radius r centred at z₀:

|f'(z₀)| = |1/(2πi) ∮_{γ_r} f(z)/(z−z₀)² dz| ≤ (1/2π) · M/r² · 2πr = M/r → 0 as r → ∞.

So f'(z₀) = 0 for all z₀, hence f is constant. ∎

> **Corollary (Fundamental Theorem of Algebra).** Every polynomial of positive degree with complex coefficients has at least one zero in ℂ.

> **Theorem (Morera's Theorem — converse of Cauchy).** If f is continuous in an open disc D and ∮_T f(z) dz = 0 for every triangle T ⊂ D, then f is holomorphic in D.

---

## Part 6: Power Series and Taylor Series

### 6.1 Power Series

> **Definition.** A *power series* is ∑_{n=0}^∞ aₙzⁿ, aₙ ∈ ℂ.

> **Theorem.** There exists 0 ≤ R ≤ ∞ (the *radius of convergence*) such that the series converges absolutely for |z| < R and diverges for |z| > R. Moreover:
>
> 1/R = lim sup_{n→∞} |aₙ|^{1/n}.

> **Theorem.** A power series with radius of convergence R > 0 is holomorphic on |z| < R, and
>
> f'(z) = ∑_{n=1}^∞ naₙz^{n−1},
>
> with the same radius of convergence.

> **Corollary.** A power series is infinitely complex differentiable in its disc of convergence, and all higher derivatives are obtained by term-by-term differentiation.

### 6.2 Taylor's Theorem ⭐

> **Theorem (Taylor Expansion).** Let f be holomorphic in an open set Ω, z₀ ∈ Ω, and suppose the closed disc D̄_{|z−z₀|}(z₀) ⊂ Ω. Then:
>
> f(z) = ∑_{n=0}^∞ f^{(n)}(z₀)/n! · (z − z₀)ⁿ.
>
> The series converges for |z − z₀| < R, where R = dist(z₀, ∂Ω).

**Key Taylor series (memorise):**

- eᶻ = ∑ zⁿ/n!, R = ∞
- 1/(1−z) = ∑ zⁿ, R = 1
- sin z = ∑ (−1)ⁿz^{2n+1}/(2n+1)!, R = ∞
- cos z = ∑ (−1)ⁿz^{2n}/(2n)!, R = ∞
- Log(1+z) = ∑ (−1)^{n+1}zⁿ/n, R = 1

### 6.3 Uniform Convergence and Holomorphic Limits

> **Theorem.** If {fₙ} is a sequence of holomorphic functions converging uniformly on compact subsets of Ω, then the limit function is holomorphic in Ω.

⚠️ This is FALSE for real differentiable functions — the uniform limit of differentiable functions need not be differentiable.

---

## Part 7: Laurent Series and Singularities ⭐

### 7.1 Laurent Series

> **Theorem (Laurent Expansion).** Let f be holomorphic in the annulus A = {z : r < |z − z₀| < R}. Then:
>
> f(z) = ∑_{n=−∞}^∞ aₙ(z − z₀)ⁿ,
>
> where aₙ = (1/2πi) ∮_γ f(ξ)/(ξ − z₀)^{n+1} dξ for any circle γ in A centred at z₀.

The sum ∑_{n=0}^∞ aₙ(z − z₀)ⁿ is the *analytic part* (or regular part). The sum ∑_{n=1}^∞ a_{−n}(z − z₀)^{−n} is the *principal part*.

**Exam technique for finding Laurent series:** Typically, you decompose f using partial fractions and then expand each piece as a geometric series in the correct region.

**Example.** Find the Laurent series of f(z) = 9/[(z−4)(z+5)] for 4 < |z| < 5.

Partial fractions: f(z) = 1/(z−4) − 1/(z+5).

For 1/(z−4) with |z| > 4: write 1/(z−4) = (1/z) · 1/(1 − 4/z) = (1/z)∑(4/z)ⁿ = ∑ 4ⁿ/z^{n+1}.

For −1/(z+5) with |z| < 5: write −1/(z+5) = (−1/5) · 1/(1 + z/5) = (−1/5)∑(−z/5)ⁿ = ∑ (−1)^{n+1}zⁿ/5^{n+1}.

Combine for the full Laurent series in 4 < |z| < 5.

### 7.2 Isolated Singularities ⭐

> **Definition.** z₀ is an *isolated singularity* of f if f is holomorphic in some punctured disc 0 < |z − z₀| < R but not at z₀ itself.

**Classification (from the Laurent series ∑ aₙ(z − z₀)ⁿ):**

| Type | Condition | Principal part |
|---|---|---|
| **Removable singularity** | aₙ = 0 for all n < 0 | Empty |
| **Pole of order m** | a_{−m} ≠ 0, aₙ = 0 for n < −m | Finite (m terms) |
| **Essential singularity** | Infinitely many aₙ ≠ 0 for n < 0 | Infinite |

⚠️ **Exam warning (2024 markers):** The definition of "pole of order m" requires BOTH that a_{−m} ≠ 0 AND that a_{−m−1} = a_{−m−2} = ⋯ = 0. Many students forget the second part.

**Equivalent characterisations:**

- **Removable:** lim_{z→z₀} f(z) exists (and is finite)
- **Pole of order m:** lim_{z→z₀} (z − z₀)^m f(z) exists and is non-zero
- **Essential:** lim_{z→z₀} f(z) does not exist (even as ∞)

> **Proposition.** If f is holomorphic at z₀ and has a zero of order m at z₀, then 1/f has a pole of order m at z₀.

### 7.3 L'Hôpital's Rule (Complex Version)

> **Theorem (L'Hôpital).** Suppose f, g are holomorphic in Ω, and z₀ ∈ Ω with f(z₀) = g(z₀) = 0 and g'(z₀) ≠ 0. Then lim_{z→z₀} f(z)/g(z) = f'(z₀)/g'(z₀).

More generally, if f and g both have zeros of order ≥ m at z₀, then lim_{z→z₀} f(z)/g(z) = lim_{z→z₀} f^{(m)}(z)/g^{(m)}(z), provided the denominator's derivative is non-zero.

This is useful for quickly identifying the type of singularity: if f(z₀) = g(z₀) = 0 and g'(z₀) ≠ 0, then f/g has a removable singularity at z₀.

### 7.4 The Identity Theorem (Isolation of Zeros)

> **Theorem.** If f is holomorphic and not identically zero on a connected open set Ω, then the zeros of f are isolated. That is, for each zero z₀, there exists r > 0 such that f(z) ≠ 0 for 0 < |z − z₀| < r.

This has the important consequence: if two holomorphic functions agree on a set with a limit point in Ω, then they agree everywhere on Ω.

---

## Part 8: Residues and the Residue Theorem ⭐⭐

### 8.1 Definition of Residue

> **Definition.** If f has an isolated singularity at z₀ with Laurent series ∑ aₙ(z − z₀)ⁿ, the *residue* of f at z₀ is
>
> Res[f, z₀] = a_{−1}.

The significance: ∮_{|z−z₀|=r} f(z) dz = 2πi · Res[f, z₀] (for small enough r).

### 8.2 Computing Residues ⭐

**Simple pole (order 1):** If f has a simple pole at z₀:

Res[f, z₀] = lim_{z→z₀} (z − z₀)f(z).

**Special case:** If f(z) = g(z)/h(z) where g(z₀) ≠ 0, h(z₀) = 0, h'(z₀) ≠ 0 (simple zero of denominator):

Res[f, z₀] = g(z₀)/h'(z₀).

**Pole of order m:** If f has a pole of order m at z₀:

Res[f, z₀] = (1/(m−1)!) · lim_{z→z₀} (d/dz)^{m−1} [(z − z₀)^m f(z)].

**By Laurent expansion:** Sometimes the easiest approach is to directly expand f in a Laurent series and read off the coefficient of (z − z₀)^{−1}.

### 8.3 The Residue Theorem ⭐

> **Theorem.** Let f be holomorphic in an open set Ω except for isolated singularities z₁, …, zₖ. Let γ ⊂ Ω be a simple, closed, piecewise-smooth curve enclosing exactly these singularities. Then:
>
> ∮_γ f(z) dz = 2πi ∑_{j=1}^k Res[f, zⱼ].

This is the central computational tool of complex analysis.

---

## Part 9: Computing Real Integrals via Residues ⭐⭐

This topic appears on every exam. There are two main types.

### 9.1 Type I: Trigonometric Integrals ∫₀^{2π} R(cos θ, sin θ) dθ

**Method:** Substitute z = e^{iθ}, so dz = ie^{iθ} dθ = iz dθ, giving dθ = dz/(iz). Then:

cos θ = (z + z⁻¹)/2,    sin θ = (z − z⁻¹)/(2i).

The integral becomes a contour integral around |z| = 1.

**Example.** Compute ∫₀^{2π} sin²θ/(2 + cos θ) dθ.

Substituting: sin²θ = −(z − 1/z)²/4 = −(z² − 1)²/(4z²), and 2 + cos θ = 2 + (z² + 1)/(2z) = (4z + z² + 1)/(2z).

The integral becomes (i/2)∮_{|z|=1} (z² − 1)² / [z²(z² + 4z + 1)] dz.

Find poles inside |z| = 1: z = 0 (order 2) and z = −2 + √3 (order 1). Compute residues and apply the residue theorem.

### 9.2 Type II: Improper Integrals ∫_{−∞}^∞ f(x) dx

**Method:** Use a semicircular contour γ = [−R, R] ∪ C_R, where C_R is the semicircle of radius R in the upper (or lower) half-plane. If f(z) → 0 fast enough on C_R as R → ∞, the contribution from C_R vanishes.

**Standard setup:** If f(z) = P(z)/Q(z) is rational with deg Q ≥ deg P + 2 and Q has no real zeros, then:

∫_{−∞}^∞ f(x) dx = 2πi · (sum of residues of f in the upper half-plane).

**For integrals involving e^{iξx}:** Close in the upper half-plane if ξ > 0 (so that |e^{iξz}| = e^{−ξy} → 0 for y > 0), and in the lower half-plane if ξ < 0.

**Example.** Compute ∫_{−∞}^∞ e^{−iξx}/(1 + x²) dx for ξ < 0.

Close in the upper half-plane. The only pole inside is z = i (simple pole).

Res[e^{−iξz}/(1 + z²), i] = e^{−iξ·i}/(2i) = e^ξ/(2i).

So ∫_{−∞}^∞ e^{−iξx}/(1 + x²) dx = 2πi · e^ξ/(2i) = πe^ξ = πe^{−|ξ|}.

⚠️ **Exam warning (2023 markers):** "Students struggled with choosing the contour, with many choosing it to pass right through a pole." Always check that your contour avoids ALL singularities.

### 9.3 Worked Example: ∫₀^{2π} 1/(2 − cos θ) dθ (from notes)

Substitute z = e^{iθ}: cos θ = (z + z⁻¹)/2, dθ = dz/(iz).

∫₀^{2π} 1/(2 − cos θ) dθ = ∮_{|z|=1} 1/[2 − (z+1/z)/2] · dz/(iz) = ∮_{|z|=1} 1/[(4z − z² − 1)/(2z)] · dz/(iz)

= ∮_{|z|=1} 2/(−iz² + 4iz − i) dz = ∮_{|z|=1} 2/[−i(z² − 4z + 1)] dz.

Roots of z² − 4z + 1 = 0: z = 2 ± √3. Only z = 2 − √3 ≈ 0.27 is inside |z| = 1.

Res at z = 2−√3: 2/[−i · 2(z − 2)] evaluated at z = 2−√3 gives 2/[−i(−2√3)] = 1/(i√3).

Result: 2πi · 1/(i√3) = 2π/√3.

---

## Part 10: Zeros, Rouché's Theorem, and the Argument Principle ⭐⭐

### 10.1 Zeros of Holomorphic Functions

> **Definition.** If f is holomorphic at z₀ and f(z₀) = 0, the *order* of the zero is the smallest n ≥ 1 such that f^{(n)}(z₀) ≠ 0.

Equivalently, f(z) = (z − z₀)ⁿ g(z) where g is holomorphic and g(z₀) ≠ 0.

> **Theorem (Isolation of zeros).** The zeros of a non-constant holomorphic function are isolated.

### 10.2 The Argument Principle

> **Theorem.** If f is holomorphic inside and on γ (simple, closed) except for poles, then:
>
> (1/2πi) ∮_γ f'(z)/f(z) dz = N − P,
>
> where N = number of zeros (counted with multiplicity) and P = number of poles (counted with order) of f inside γ.

### 10.3 Rouché's Theorem ⭐

> **Theorem (Rouché).** Let f and g be holomorphic in an open set Ω, and let γ ⊂ Ω be a simple, closed, piecewise-smooth curve containing only points of Ω in its interior. If
>
> |g(z)| < |f(z)| for all z ∈ γ,
>
> then f and f + g have the same number of zeros (counted with multiplicity) inside γ.

**Exam technique:** To count zeros of a function h(z) inside a contour γ:

1. Split h = f + g where one piece "dominates" on γ
2. Verify |g(z)| < |f(z)| on γ (usually using the triangle inequality)
3. Count zeros of f (the simpler function) inside γ

**Example (exam-style).** Show that z³ + 5z + 1 has no zeros with |z| ≤ 1.

Let f(z) = 5z and g(z) = z³ + 1. On |z| = 1: |f(z)| = 5, and |g(z)| = |z³ + 1| ≤ |z|³ + 1 = 2 < 5 = |f(z)|.

By Rouché, z³ + 5z + 1 has the same number of zeros as 5z inside |z| = 1, which is exactly 1 zero (at z = 0, but wait — we need to check: z³ + 5z + 1 at z = 0 gives 1 ≠ 0, so the zero is perturbed to somewhere inside the disc).

⚠️ **Exam warning (2025 markers):** "For |z| > 1, Rouché's theorem cannot be directly applied. One has to show there are no roots exactly on |z| = 1, then subtract the number inside from the total."

**Counting zeros for |z| > 1:** Apply Rouché on a large circle |z| = R to find the total number of zeros, then subtract those inside |z| ≤ 1.

**Example (roots in an annulus — from notes).** Show all roots of w(z) = z⁷ − 2z² + 8 lie in 1 < |z| < 2.

*On |z| = 2:* Let f(z) = z⁷, g(z) = −2z² + 8. Then |f(z)| = 128 and |g(z)| ≤ 8 + 8 = 16 < 128. By Rouché, w has 7 zeros inside |z| = 2. Since w has degree 7, ALL roots are inside |z| ≤ 2.

*On |z| = 1:* Let f(z) = 8, g(z) = z⁷ − 2z². Then |g(z)| ≤ 1 + 2 = 3 < 8 = |f(z)|. By Rouché, w has the same number of zeros as 8 inside |z| = 1, which is 0.

Therefore all 7 roots lie in the annulus 1 < |z| < 2.

---

## Part 11: Maximum Modulus Principle ⭐

### 11.1 Open Mapping Theorem

> **Theorem.** If f is a non-constant holomorphic function on an open set Ω, then f(Ω) is open. That is, f maps open sets to open sets.

**Proof (examinable — uses Rouché).** Let w₀ = f(z₀) ∈ f(Ω). Choose δ > 0 so that the closed disc D̄_δ(z₀) ⊂ Ω and f(z) ≠ w₀ on the circle |z − z₀| = δ (possible since zeros are isolated). Let ε = min_{|z−z₀|=δ} |f(z) − w₀| > 0. We claim D_ε(w₀) ⊂ f(Ω).

For any w₁ with |w₁ − w₀| < ε, set F(z) = f(z) − w₀ and G(z) = w₀ − w₁. On |z − z₀| = δ: |F(z)| ≥ ε > |w₀ − w₁| = |G(z)|. By Rouché, F + G = f(z) − w₁ has the same number of zeros as F(z) = f(z) − w₀ inside the disc, which has at least one (at z₀). So w₁ ∈ f(Ω). ∎

### 11.2 Maximum Modulus Principle

> **Theorem (Maximum Modulus Principle).** If f is a non-constant holomorphic function on an open set Ω ⊂ ℂ, then |f| cannot attain a maximum in Ω.

**Proof (examinable — short, uses Open Mapping Theorem).** Suppose |f| attains a maximum at z₀ ∈ Ω, i.e. |f(z)| ≤ |f(z₀)| for all z in a neighbourhood. Let D ⊂ Ω be a small open disc around z₀. Since f is an open mapping, f(D) is open and contains f(z₀). So f(D) contains a disc around f(z₀), which includes points with |w| > |f(z₀)| — contradiction. ∎

> **Corollary.** If f is holomorphic on a bounded open set Ω and continuous on its closure Ω̄, then max_{z∈Ω̄} |f(z)| is attained on the boundary ∂Ω.

### 11.3 Schwarz Lemma ⭐

> **Lemma (Schwarz).** If f is holomorphic in D = {z : |z| < 1}, f(0) = 0, and |f(z)| ≤ 1 for all z ∈ D, then |f(z)| ≤ |z| for all z ∈ D.

**Proof (examinable — from PS7).** Consider g(z) = f(z)/z. Since f(0) = 0, g is holomorphic in D. On |z| = ρ < 1: |g(z)| = |f(z)|/ρ ≤ 1/ρ. By the maximum modulus principle, |g(z)| ≤ 1/ρ for all z with |z| ≤ ρ. Fixing z and letting ρ → 1 gives |g(z)| ≤ 1, i.e. |f(z)| ≤ |z|. ∎

⚠️ **Exam warning (2022 markers):** "A couple of students simply tried to invoke Schwarz Lemma, but that was what the question was asking you to prove." You need to be able to reproduce the proof, not just state the result.

---

## Part 12: Harmonic Functions ⭐

### 12.1 Definition and Connection to Holomorphic Functions

> **Definition.** A twice continuously differentiable function u : Ω → ℝ (where Ω ⊆ ℝ²) is *harmonic* if Δu = ∂²u/∂x² + ∂²u/∂y² = 0.

> **Theorem.** If f = u + iv is holomorphic, then both u and v are harmonic.

**Proof (short).** From C–R: uₓ = vᵧ and uᵧ = −vₓ. Differentiating: uₓₓ = vᵧₓ and uᵧᵧ = −vₓᵧ. By equality of mixed partials: uₓₓ + uᵧᵧ = vᵧₓ − vₓᵧ = 0. ∎

### 12.2 Harmonic Conjugates ⭐

> **Definition.** Given a harmonic function u, a function v such that f = u + iv is holomorphic is called a *harmonic conjugate* of u.

> **Theorem (Existence of Harmonic Conjugates).** If u is harmonic in an open disc D (or more generally, in a simply connected domain), then there exists a harmonic conjugate v, defined by the line integral:
>
> v(x, y) = ∫_γ (−∂u/∂y dx + ∂u/∂x dy),
>
> where γ is any path from a fixed base point (x₀, y₀) to (x, y) in D. The integral is path-independent because u is harmonic and D is simply connected.

**Exam technique for finding v given u (the practical method):**

1. Compute uₓ and uᵧ
2. Use C–R: vₓ = −uᵧ, vᵧ = uₓ
3. Integrate vᵧ = uₓ with respect to y to get v(x,y) = ∫ uₓ dy + C(x)
4. Differentiate with respect to x and use vₓ = −uᵧ to find C'(x)
5. Integrate to get C(x) (up to a constant)

**Example.** Given u(x,y) = x² − y² (harmonic: uₓₓ = 2, uᵧᵧ = −2, sum = 0).

Step 1: uₓ = 2x, uᵧ = −2y.

Step 2: vᵧ = uₓ = 2x, so v = 2xy + C(x).

Step 3: vₓ = 2y + C'(x) = −uᵧ = 2y, so C'(x) = 0, hence C(x) = const.

Result: v = 2xy + c, and f(z) = (x² − y²) + i(2xy) = z².

**Example 2 (exam-style, from notes).** Given u(x, y) = x³ − 3xy² + y, with f(1) = 1 + i.

Step 1: Check harmonic: uₓₓ = 6x, uᵧᵧ = −6x, so Δu = 0. ✓

Step 2: uₓ = 3x² − 3y², uᵧ = −6xy + 1.

Step 3: vᵧ = uₓ = 3x² − 3y², so v = 3x²y − y³ + C(x).

Step 4: vₓ = 6xy + C'(x) = −uᵧ = 6xy − 1, so C'(x) = −1, giving C(x) = −x + c.

Step 5: v(x,y) = 3x²y − y³ − x + c. Using f(1) = 1 + i: u(1,0) + iv(1,0) = 1 + i(−1+c) = 1+i gives c = 2.

So f(z) = (x³ − 3xy² + y) + i(3x²y − y³ − x + 2) = z³ + iz − i + 2i = z³ + iz + i.

**Key identity involving the Laplacian and Wirtinger operators:**

4(∂/∂z)(∂/∂z̄) = Δ.

If f is holomorphic then Δ|f(z)|² = 4|f'(z)|² (from PS8). This follows because u and v are harmonic, so Δ(u² + v²) reduces to 2(uₓ² + uᵧ² + vₓ² + vᵧ²) = 4|f'|² by C–R.

⚠️ **Exam warning (2025 markers):** "Schwarz reflection is not applicable because f is not necessarily real-valued on the real line." Be careful about when the Schwarz reflection principle applies.

---

## Part 13: Conformal Mappings and Möbius Transformations

### 13.1 Conformal Mappings

> **Theorem (Angle Preservation).** Let f be holomorphic in Ω, and let γ₁, γ₂ be curves in Ω intersecting at z₀ with z₁'(0), z₂'(0), f'(z₀) all non-zero. Then the angle between f(γ₁) and f(γ₂) at f(z₀) equals the angle between γ₁ and γ₂ at z₀.

The proof follows from the chain rule: (f ∘ zⱼ)'(0) = f'(z₀)zⱼ'(0), so arg(f ∘ z₂)'(0) − arg(f ∘ z₁)'(0) = arg z₂'(0) − arg z₁'(0).

> **Definition.** f is *conformal* in Ω if it is holomorphic and f'(z) ≠ 0 for all z ∈ Ω.

> **Theorem.** If f : Ω → ℂ is a holomorphic local injection, then f'(z) ≠ 0 for all z ∈ Ω. Moreover, the inverse f⁻¹ is also holomorphic with (f⁻¹)'(w₀) = 1/f'(f⁻¹(w₀)).

### 13.2 Möbius Transformations ⭐

> **Definition.** A *Möbius transformation* is T(z) = (az + b)/(cz + d) where a, b, c, d ∈ ℂ and ad − bc ≠ 0.

**Matrix representation.** T corresponds to the matrix (a b; c d). Composition of Möbius transformations corresponds to matrix multiplication. The inverse is T⁻¹(w) = (dw − b)/(−cw + a).

**Three basic types (every Möbius transformation is a composition of these):**

- **Scaling/rotation:** z ↦ az (rotation by arg(a), dilation by |a|)
- **Translation:** z ↦ z + b
- **Inversion:** z ↦ 1/z

> **Corollary.** Möbius transformations map circles and lines to circles and lines (where lines are considered circles of infinite radius).

**Key properties:**

- A Möbius transformation is uniquely determined by the images of three points
- They are conformal on their domain
- They preserve the cross-ratio

### 13.3 The Cross-Ratio ⭐

> **Theorem (Cross-Ratio Preservation).** If T is a Möbius transformation mapping distinct points (z₁, z₂, z₃) to (w₁, w₂, w₃), then for all z:
>
> (z − z₁)(z₂ − z₃) / [(z − z₃)(z₂ − z₁)] = (w − w₁)(w₂ − w₃) / [(w − w₃)(w₂ − w₁)].

**Exam technique:** To find the Möbius transformation mapping z₁ ↦ w₁, z₂ ↦ w₂, z₃ ↦ w₃, set the cross-ratios equal and solve for w in terms of z. If one of the points is ∞, take a limit (e.g. the factor (z₂ − z₃)/(z₂ − z₁) becomes (1 − z₃/z₂)/(1 − z₁/z₂) → 1 as z₂ → ∞, leaving (z − z₁)/(z − z₃)).

**Example (2023 exam style).** Find T mapping 1 ↦ −1, i ↦ 0, −1 ↦ 1.

Cross-ratio: (z − 1)(i + 1)/[(z + 1)(i − 1)] = (w + 1)(0 − 1)/[(w − 1)(0 + 1)].

Simplifying the left side: (z − 1)/(z + 1) · (i+1)/(i−1) = (z−1)/(z+1) · (−1) = −(z−1)/(z+1) = (w+1)/(w−1).

Solving: w = (z − i)/(zi − 1).

### 13.4 The Half-Plane ↔ Unit Disc Map ⭐

The upper half-plane H = {z : Im(z) > 0} is conformally equivalent to the unit disc D = {w : |w| < 1}.

**H → D:** f(z) = (i − z)/(i + z).

**D → H:** g(w) = i(1 − w)/(1 + w).

**Verification that f maps H to D:** For z = x + iy with y > 0:

|f(z)|² = |i − z|²/|i + z|² = [x² + (y−1)²]/[x² + (y+1)²] < 1 (since y > 0).

**Verification that g maps D to H:** For w = u + iv with |w| < 1:

Im(g(w)) = (1 − u² − v²)/[(1+u)² + v²] > 0 (since |w|² = u² + v² < 1).

The boundary ∂H = ℝ maps to ∂D = {|w| = 1}, and ±∞ both map to −1.

**For real a, b, c, d with ad − bc > 0:** T(z) = (az+b)/(cz+d) maps the upper half-plane to itself, since Im(T(z)) = (ad − bc)Im(z)/|cz + d|² > 0 when Im(z) > 0.

---

## Part 14: Schwarz Reflection Principle

> **Definition.** An open set Ω ⊂ ℂ is *symmetric with respect to the real axis* if z ∈ Ω ⟺ z̄ ∈ Ω. Define Ω⁺ = {z ∈ Ω : Im z > 0}, Ω⁻ = {z ∈ Ω : Im z < 0}, I = {z ∈ Ω : Im z = 0}.

> **Theorem (Symmetry Principle).** Let f₊ and f₋ be holomorphic in Ω⁺ and Ω⁻ respectively, both extending continuously to I with f₊(x) = f₋(x) for all x ∈ I. Then the function f defined by f₊ on Ω⁺ ∪ I and f₋ on Ω⁻ is holomorphic on all of Ω.

The proof uses Morera's theorem: for any triangle T in Ω, if T doesn't cross I then ∮_T f dz = 0 trivially. If T does cross I, approximate by triangles lying in one half-plane and use continuity.

> **Theorem (Schwarz Reflection Principle).** Suppose f is holomorphic in Ω⁺, extends continuously to I, and **is real-valued on I**. Then there exists a holomorphic function F on Ω such that F|_{Ω⁺} = f. The extension is given by F(z) = f̄(z̄) for z ∈ Ω⁻.

**Proof sketch.** Define f₋(z) = f̄(z̄) for z ∈ Ω⁻. Since f is holomorphic in Ω⁺ with Taylor series f(z) = ∑ aₙ(z − z₀)ⁿ, we get f₋(z) = ∑ āₙ(z − z̄₀)ⁿ, which is holomorphic in Ω⁻. Since f is real-valued on I, f₊(x) = f(x) = f̄(x) = f₋(x) for x ∈ I. The Symmetry Principle then gives holomorphicity on all of Ω. ∎

**Example.** The square root f(z) = √z = r^{1/2}e^{iθ/2} (for z = re^{iθ}, θ ∈ (0,π)) is holomorphic in the upper half-plane, continuous and real-valued on [0,∞). By Schwarz Reflection, it extends holomorphically to the lower half-plane via F(z) = f̄(z̄) = r^{1/2}e^{iθ/2} for θ ∈ (−π, 0).

⚠️ **Exam warning (2025):** The hypothesis that f is real-valued on the real axis is essential. Don't apply this theorem when it's not met.

---

## Summary of Exam Topic Frequency (2022–2025)

| Topic | 2022 | 2023 | 2024 | 2025 | Priority |
|-------|------|------|------|------|----------|
| Cauchy–Riemann equations | ✓ | ✓ | ✓ | ✓ | **Essential** |
| Laurent series & singularity classification | ✓ | ✓ | ✓ | ✓ | **Essential** |
| Residues & Residue Theorem | ✓ | ✓ | ✓ | ✓ | **Essential** |
| Rouché's theorem (counting zeros) | ✓ | ✓ | ✓ | ✓ | **Essential** |
| Computing real integrals via residues | ✓ | ✓ | ✓ | ✓ | **Essential** |
| Complex logarithm / powers | — | ✓ | ✓ | ✓ | **High** |
| Harmonic functions / conjugates | — | ✓ | ✓ | ✓ | **High** |
| Maximum Modulus Principle | ✓ | ✓ | ✓ | — | **High** |
| Schwarz Lemma | ✓ | — | ✓ | — | **Medium** |
| Möbius transformations | ✓ | ✓ | — | — | **Medium** |
| Cauchy's Integral Formula | ✓ | ✓ | ✓ | ✓ | **High** (implicit) |
| Liouville's Theorem | ✓ | — | — | — | **Medium** |

---

## Key Examiner Warnings (from markers' comments)

1. **cos z is not bounded for complex z** (2025) — don't assume real-variable bounds carry over
2. **Rouché cannot be directly applied for |z| > R** (2025) — count total zeros on a large circle, subtract those inside
3. **Schwarz reflection requires f real on ℝ** (2025) — verify the hypothesis before applying
4. **Pole of order m requires both a_{−m} ≠ 0 and all earlier coefficients zero** (2024) — state the full definition
5. **|y| ≤ R does NOT imply 1/|y| ≤ 1/R** (2024) — careful with inequality directions
6. **Know the principal value of log** (2024) — many students were unfamiliar
7. **Don't use Liouville when f is not entire** (2023) — it only works for bounded entire functions
8. **Don't choose contours passing through poles** (2023)
9. **The Schwarz Lemma must be proved, not just stated** (2022) — learn the short proof

---

## Recommended Problem Sheet Questions

**PS1 (Complex numbers, basics):** Q3 (modulus and argument), Q5 (geometric descriptions — builds intuition for the complex plane).

**PS2 (Logarithm, powers, CR equations):** ⭐ Q1b (computing i^i, 2^i — exam staple), Q1c (Log computations), Q2 (CR equations for specific functions — drill this), Q1d (Bernoulli's paradox — tests multi-valuedness understanding).

**PS3 (Contour integrals, simply connected domains):** ⭐ **Q1** (∮ 1/(z²−4) dz — partial fractions + deformation), ⭐ **Q2** (∮ 1/(z³−1) dz — identifying which poles are inside the contour), Q3 (simply connected domains), **Q8** (∮ dz/[z(z−1)⋯(z−k)] — systematic use of CIF with multiple poles).

**PS4 (Cauchy integral formula, higher derivatives):** ⭐ **Q1** (∮ z²/(z−1)ⁿ dz — different answers for n=1,2,3,>3), **Q2b** (∮ sin z/(z+2)³ dz — CIF for second derivative), Q3 (ML-inequality: polynomials can't approximate 1/z on unit circle), **Q5** (|f(z)| ≤ C(1+|z|)ⁿ ⟹ f is polynomial of degree ≤ n), Q6 (entire + bounded Re ⟹ constant, via e^f and Liouville).

**PS5 (Laurent series, singularities):** ⭐ **Q1** (Laurent series about specific points), **Q2** (Laurent in a specific annulus), ⭐ **Q3** (Laurent in three different regions — exam classic), Q5a (zero of order m ⟹ pole of order m for 1/f).

**PS6 (Real integrals via residues):** ⭐ **Q1** (∫₀^{2π} sin²θ/(2+cosθ) dθ — the z=e^{iθ} method appears every exam), Q2 (holomorphic continuation — tests analytic continuation).

**PS7 (Maximum modulus, Schwarz lemma, improper integrals):** ⭐ **Q3** (prove Schwarz Lemma — tested 2022, learn this proof!), ⭐ **Q4** (∫ e^{−iξx}/(1+x²) dx — semicircular contour), Q1 (max modulus applied to polynomial on unit circle), Q2 (no holomorphic f with |f(z)| = e^{|z|}  — uses max mod).

**PS8 (Harmonic functions, Laplacian):** **Q1** (Δ = 4∂_z∂_{z̄} and Δ|f|² = 4|f'|²), Q3–Q4 (Poisson integral formula — appeared in 2024 exam).

---

*End of revision guide. Good luck with the exam!*
