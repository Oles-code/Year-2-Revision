# Analysis II — Term 1: Comprehensive Exam Revision Guide

**MATH50001 — Imperial College London**

---

**How to use this guide.** This document distills the key examinable content from the Term 1 lecture notes. It is structured around the topics that appear most frequently and in the most challenging contexts across the 2022–2025 final exams. Proofs are included where they are short enough to plausibly appear on an exam. At the end, specific problem sheet questions are referenced that practise the most exam-relevant skills.

**Exam structure.** The paper has 6 questions of equal weight (20 marks each). Questions 1–3 cover Term 1 (this material). Typically: Q1 tests differentiation in Rⁿ (continuity, differentiability, IFT/Inverse FT); Q2 tests metric spaces (metric verification, completeness, compactness); Q3 tests topology of function spaces (C[0,1] with d∞, open/closed, connected, compact sets).

---

# Part 1: Differentiation in Higher Dimensions

## 1.1 Euclidean Spaces, Norms, and Inner Products

The Euclidean space Rⁿ consists of ordered n-tuples x = (x¹, x², ..., xⁿ). It is a vector space with component-wise addition and scalar multiplication.

> **Definition (Inner Product on Rⁿ).**
> For x, y ∈ Rⁿ, the inner product is ⟨x, y⟩ = Σᵢ xⁱyⁱ. This satisfies: ⟨x,y⟩ = ⟨y,x⟩, linearity in each argument, and ⟨x,x⟩ ≥ 0 with equality iff x = 0.

> **Definition (Norm on Rⁿ).**
> ‖x‖ = √⟨x,x⟩ = (Σᵢ (xⁱ)²)^(1/2). This satisfies:
> - (i) ‖x‖ ≥ 0 with equality iff x = 0
> - (ii) ‖λx‖ = |λ|‖x‖
> - (iii) ‖x + y‖ ≤ ‖x‖ + ‖y‖ (triangle inequality)

> **Cauchy–Schwarz Inequality.**
> |⟨x, y⟩| ≤ ‖x‖ · ‖y‖ for all x, y ∈ Rⁿ.

> **Reverse Triangle Inequality.**
> |‖x‖ − ‖y‖| ≤ ‖x − y‖ for all x, y ∈ Rⁿ.

⚠ **Exam Tip:** The reverse triangle inequality is used frequently to prove continuity of the norm function g(x) = ‖x‖ — this appeared on the 2025 exam. The key step: |g(z) − g(w)| = |‖z‖ − ‖w‖| ≤ ‖z − w‖ < ε.

---

## 1.2 Open Sets in Rⁿ

> **Definition (Open Set in Rⁿ).**
> A set U ⊆ Rⁿ is open if for every u ∈ U, there exists δ > 0 such that Bδ(u) ⊆ U, where Bδ(u) = {x ∈ Rⁿ : ‖x − u‖ < δ}.

To prove a set Ω = {(x,y,z) : x > 0, y > 0, z > 0} is open (a common exam question), fix (x,y,z) ∈ Ω and let δ = min{x, y, z} > 0. Then for any (x', y', z') with ‖(x,y,z) − (x',y',z')‖ < δ, each coordinate satisfies |xⁱ − x'ⁱ| < δ, so x'ⁱ > 0.

⚠ **Exam Tip:** Identifying open/closed/neither for given sets appears on almost every exam. Remember: the preimage of an open set under a continuous map is open. For a set defined by a strict inequality like {y ≠ x sin(y)}, think of it as the preimage of R\{0} under the continuous map (x,y) ↦ y − x sin(y).

---

## 1.3 Continuity of Maps on Rⁿ

> **Definition (Continuity at a Point).**
> A map f: Ω → Rᵐ (where Ω ⊆ Rⁿ is open) is continuous at p ∈ Ω if for every ε > 0 there exists δ > 0 such that ‖x − p‖ < δ implies ‖f(x) − f(p)‖ < ε.

**Common technique for showing continuity at the origin:** Use the squeeze theorem. If f(x,y,z) = xyz/(x²+y²+z²), bound using |x| ≤ ‖(x,y,z)‖ to get |f| ≤ ‖(x,y,z)‖ → 0. This is exactly the technique used in the 2023 Q1(a) solution.

**Common technique for showing discontinuity:** Approach the point along two different paths and get different limits. E.g., along y = 0 versus y = x.

---

## 1.4 Differentiability

> **Definition (Differentiability of f: Ω → Rᵐ at p ∈ Ω).**
> f is differentiable at p if there exists a linear map Λ ∈ L(Rⁿ; Rᵐ) such that
>
> lim\_{h→0} ‖f(p+h) − f(p) − Λ[h]‖ / ‖h‖ = 0.
>
> We write Df(p) = Λ (the derivative / Jacobian of f at p).

**Key Fact:** Differentiability implies continuity (but not vice versa).

> **Uniqueness of the Derivative.**
> If it exists, the derivative Df(p) is unique.

### 1.4.1 Partial Derivatives

> **Definition (Partial Derivative).**
> Dⱼf(p) = lim\_{t→0} [f(p + teⱼ) − f(p)] / t, where eⱼ is the j-th standard basis vector.

⚠ **Exam Tip:** On the exam, when asked to find partial derivatives at the origin for a piecewise-defined function, you MUST use the limit definition — you cannot just 'plug in' to the formula that applies away from the origin. This was noted as a common error in the 2024 markers' comments.

> **Theorem 1.12 (Continuous Partials ⟹ Differentiable).**
> Let Ω ⊆ Rⁿ be open, and f: Ω → Rᵐ. If all first partial derivatives Dⱼfⁱ exist and are continuous on Ω, then f is differentiable at every point in Ω.

**Why this matters:** This is the standard way to show differentiability. Compute the partials, check they're continuous (usually by quotient/product of continuous functions). Then cite this theorem. If the domain is open, you're done.

> **Jacobian Matrix.**
> When f is differentiable, its derivative Df(p) is represented by the m×n Jacobian matrix whose (i,j)-entry is Dⱼfⁱ(p).

**Proving differentiability from the definition:** Sometimes you must work directly. The strategy is:
1. Guess the linear map Λ (usually from partial derivatives).
2. Show ‖f(p+h) − f(p) − Λ[h]‖ / ‖h‖ → 0. Bound the numerator, often using ‖h‖² or similar.

See the 2024 exam Q1(c) for a clean example.

---

## 1.5 Chain Rule

> **Chain Rule (Theorem 1.8).**
> If g: Ω → Rᵐ is differentiable at p, and f: Ω' → Rˡ is differentiable at g(p), where g(Ω) ⊆ Ω', then f ∘ g is differentiable at p and
>
> D(f ∘ g)(p) = Df(g(p)) ∘ Dg(p).
>
> In matrix form: the Jacobian of the composition is the product of the Jacobians.

---

## 1.6 Schwarz's Theorem (Symmetry of Mixed Partials)

> **Schwarz's Theorem (Theorem 1.13).**
> If f: Ω → R has continuous second partial derivatives on Ω (i.e., f is C²), then DᵢDⱼf = DⱼDᵢf for all i, j.

Problem Sheet 4, Exercise 4.2 gives a beautiful counterexample: a function where the mixed partials exist but are NOT equal at the origin, because the second partials are not continuous there.

---

## 1.7 Inverse Function Theorem

> **Inverse Function Theorem (Theorem 1.15).**
> Let Ω ⊆ Rⁿ be open, f: Ω → Rⁿ be continuously differentiable, and q ∈ Ω such that det Df(q) ≠ 0. Then there exist open sets U ∋ q and V ∋ f(q) such that:
> - (i) f: U → V is a bijection,
> - (ii) f⁻¹: V → U is continuously differentiable,
> - (iii) Df⁻¹(y) = [Df(f⁻¹(y))]⁻¹ for all y ∈ V.

**How to apply it (exam recipe):**
1. Show f is C¹ (usually all partials are continuous).
2. Compute the Jacobian Df at the given point.
3. Check det Df ≠ 0.
4. Conclude f is locally invertible by the IFT. If asked for Df⁻¹, use formula (iii).

⚠ **Exam Tip:** This appears in nearly every exam (2022 Q1 via IFT, 2023 Q1(c), 2025 Q1(c)). Common pitfall from the 2023 markers: "Many students worked hard to find their way around the key idea, and spent more time than required." Know the theorem statement cold.

---

## 1.8 Implicit Function Theorem

> **Implicit Function Theorem — Low-dimensional (Theorem 1.16).**
> Let Ω ⊆ R² be open, F: Ω → R continuously differentiable, and (x', y') ∈ Ω with F(x', y') = 0 and D₂F(x', y') ≠ 0. Then there exist open sets A ∋ x' and B ∋ y', and a C¹ map f: A → B such that:
>
> F(x, y) = 0 for (x, y) ∈ A × B  ⟺  y = f(x) for x ∈ A.

> **General Implicit Function Theorem (Theorem 1.17).**
> Let Ω ⊆ Rⁿ, Ω' ⊆ Rᵐ be open, f: Ω × Ω' → Rᵐ be C¹, and (a,b) ∈ Ω × Ω' with f(a,b) = 0. If the m×m matrix [Dₙ₊ⱼfⁱ(a,b)] is invertible, then there exist open sets A ∋ a, B ∋ b, and a C¹ map g: A → B with f(x, g(x)) = 0 for all x ∈ A.

**Computing derivatives of the implicit function:** Differentiate the relation F(x, y, g(x,y)) = 0 with respect to x using the chain rule. For example, in the 2022 exam Q1(d), this gives D₁g(1,1) = −D₁h/D₃h evaluated at the known point.

⚠ **Exam Tip:** The IFT appeared on the 2022 exam explicitly. The key step students miss: verifying that the partial derivative with respect to the "solved" variable is non-zero. Always check D₂F ≠ 0 (or det M ≠ 0 in higher dimensions).

---
---

# Part 2: Metric Spaces and Topology

## 2.1 Metrics — Definitions and Verification

> **Definition (Metric).**
> A metric on a set X is a function d: X × X → R satisfying:
> - **(M1)** d(x,y) ≥ 0 for all x, y, and d(x,y) = 0 ⟺ x = y  (positivity + separation)
> - **(M2)** d(x,y) = d(y,x) for all x, y  (symmetry)
> - **(M3)** d(x,y) ≤ d(x,z) + d(z,y) for all x, y, z  (triangle inequality)

**This is THE most commonly tested topic.** Every single exam from 2022–2025 asks you to verify that some function is (or isn't) a metric. The integral-based metric d(x,y) = |∫ₓʸ h(t) dt| for some positive function h appears in 2023, 2024, and 2025.

### 2.1.1 Exam recipe: Verifying an integral-based metric

**Given d(x,y) = |∫ₓʸ h(t) dt| where h > 0:**

**(M1)** Since h > 0, the integral |∫ₓʸ h(t) dt| ≥ 0, and equals 0 iff x = y (because h is strictly positive, the integral over a non-degenerate interval is strictly positive).

**(M2)** d(x,y) = |∫ₓʸ h(t) dt| = |−∫ʸₓ h(t) dt| = |∫ʸₓ h(t) dt| = d(y,x).

**(M3)** d(x,z) = |∫ₓᶻ h(t) dt| = |∫ₓʸ h(t) dt + ∫ʸᶻ h(t) dt| ≤ |∫ₓʸ h(t) dt| + |∫ʸᶻ h(t) dt| = d(x,y) + d(y,z). This uses the triangle inequality for absolute values.

### 2.1.2 Showing a function is NOT a metric

The most common failure is the triangle inequality. For f(x,y) = |x−y|², try x=0, y=4, z=2: f(0,4) = 16, but f(0,2) + f(2,4) = 4 + 4 = 8 < 16. This was tested in 2022 Q2(a).

### 2.1.3 Showing √|x−y| IS a metric

For g(x,y) = √|x−y|, the triangle inequality requires: √|x−y| ≤ √|x−z| + √|z−y|. Start from |x−y| ≤ |x−z| + |z−y| ≤ (√|x−z| + √|z−y|)². Take square roots (both sides non-negative).

### 2.1.4 Norm-induced metrics and when d is NOT from a norm

A norm ‖·‖ induces a metric via d(x,y) = ‖x−y‖. A metric comes from a norm only if d(0, 2v) = 2·d(0, v) for all v. The 2025 Q2(b) tests exactly this: showing that d\*(x,y) = |∫₀ˣ e^{−t²} dt − ∫₀ʸ e^{−t²} dt| does NOT satisfy homogeneity because e^{−t²} is not constant.

---

## 2.2 Key Metric Spaces for the Exam

**(C[0,1], d∞):** The space of continuous functions f: [0,1] → R with d∞(f,g) = sup\_{t ∈ [0,1]} |f(t) − g(t)|. This is the single most important metric space for the exam.

**Discrete metric:** d\_disc(x,y) = 0 if x = y, 1 if x ≠ y. Every subset is both open and closed. Only eventually constant sequences converge. Heine-Borel fails: sets can be closed and bounded but not compact.

**Rail metric on R²:** d\_rail(x,y) = ‖x−y‖ if x, y are collinear through the origin, otherwise ‖x‖ + ‖y‖. Tested in 2023 Q2(c).

---

## 2.3 Open and Closed Sets in Metric Spaces

> **Definition (Open Set).**
> U ⊆ X is open in (X,d) if for every u ∈ U, there exists δ > 0 with Bδ(u) ⊆ U.

> **Definition (Closed Set).**
> V ⊆ X is closed in (X,d) if for every sequence (xₙ) in V converging in (X,d), the limit belongs to V.
>
> *Equivalently: V is closed ⟺ X\V is open.*

> **Key Properties of Open/Closed Sets.**
> - Any union of open sets is open. Any finite intersection of open sets is open.
> - Any intersection of closed sets is closed. Any finite union of closed sets is closed.
> - ∅ and X are both open and closed in any metric space.

A set can be: open but not closed, closed but not open, both, or neither. "Open" is NOT the opposite of "closed."

---

## 2.4 Interior, Closure, Boundary, Limit Points

> **Definition (Interior point of V):** x ∈ X such that ∃δ > 0 with Bδ(x) ⊆ V. The interior Int(V) is the set of all interior points.
>
> **Definition (Limit point of V):** x ∈ X such that ∀δ > 0, (Bδ(x) ∩ V)\{x} ≠ ∅.
>
> **Definition (Closure of V):** V̄ = V ∪ {limit points of V}. Equivalently, the smallest closed set containing V.
>
> **Definition (Boundary of V):** ∂V = V̄ \ Int(V). Equivalently, x ∈ ∂V if every ball about x meets both V and X\V.

**Worked example (2022 Q2(b) style).** For Y = R² \ {(x,0) : x ∈ Q}:
- **Closure** = R² (every point in R² is a limit point by density of irrationals).
- **Interior** = {(x,y) : y ≠ 0} (points on the x-axis with irrational x-coordinate are not interior points, because every ball around them contains rational x-axis points).
- **Boundary** = {(x,y) : y = 0} (the entire x-axis).

---

## 2.5 Continuous Maps Between Metric Spaces

> **Definition (Continuity in Metric Spaces).**
> f: (X, dₓ) → (Y, d\_Y) is continuous at x ∈ X if: ∀ε > 0 ∃δ > 0 such that dₓ(x', x) < δ ⟹ d\_Y(f(x'), f(x)) < ε.
>
> *Equivalent: f is continuous iff f⁻¹(U) is open in X for every open set U in Y.*
>
> *Equivalent: f is continuous iff f⁻¹(V) is closed in X for every closed set V in Y.*

⚠ **Exam Tip:** The preimage characterisation is extremely useful on exams. To show a set like V = {g ∈ C[0,1] : g(0)+g(1) ∈ [0,1]} is closed, show it's the preimage of the closed set [0,1] under the continuous map Φ(g) = g(0)+g(1). This was 2022 Q3(a)–(b).

---

## 2.6 Connectedness

> **Definition (Connected Set).**
> T ⊆ X is **disconnected** if there exist open sets U, V in X with: U ∩ V = ∅, T ⊆ U ∪ V, T ∩ U ≠ ∅, T ∩ V ≠ ∅.
>
> T is **connected** if it is not disconnected.

> **Theorem 2.26 (Continuous Image of Connected is Connected).**
> If f: (A₁, d₁) → (A₂, d₂) is continuous and S ⊆ A₁ is connected, then f(S) is connected.

*Proof.* Suppose f(S) is disconnected: ∃ open U, V in A₂ with U ∩ V = ∅, f(S) ⊆ U ∪ V, f(S) ∩ U ≠ ∅, f(S) ∩ V ≠ ∅. Then U' = f⁻¹(U) and V' = f⁻¹(V) are open (f continuous), disjoint, cover S, and both meet S. This contradicts S connected. □

> **Theorem 2.25 (Intervals are Connected).**
> For a < b, [a,b] is connected in (R, d₁). In fact, any interval in R is connected, and every connected subset of R is an interval.

> **Definition (Path-Connected).**
> (X, d) is path-connected if for any a, b ∈ X there exists a continuous map γ: [0,1] → X with γ(0) = a and γ(1) = b.

> **Theorem 2.30 (Path-Connected ⟹ Connected).**
> If (X, d) is path-connected, then it is connected. The converse is false (the topologist's sine curve is connected but not path-connected).

**Exam strategy for showing connectedness:** Show path-connectedness instead — it's usually much easier. For a convex set like D = {f ∈ C[0,1] : d∞(f,0) ≤ 1}, the path γ(t) = (1−t)f + tg works. This appeared in 2023 Q3(c) and 2025 Q3(b). You need to verify:
1. γ(t) ∈ D for all t ∈ [0,1] (use the triangle inequality / convexity of the ball).
2. γ is continuous from [0,1] to (C[0,1], d∞) (bound d∞(γ(t), γ(s)) ≤ 2M|t−s| where M bounds f and g).

**Strategy for showing disconnectedness:** Either find two suitable open sets U, V explicitly, or find a continuous surjection onto {0, 1} (or any disconnected set). For W = {g ∈ C[0,1] : g(0)+g(1) ∈ Q}, compose with a function h that separates Q (2022 Q3(c)).

⚠ **Exam Tip:** The 2024 markers noted a common error: "The most common major error was claiming falsely that a connected set is path-connected." While path-connected ⟹ connected, the converse fails. However, for convex subsets of normed spaces, both notions coincide.

---

## 2.7 Compactness

> **Definition (Compact Set).**
> Y ⊆ X is compact if every open cover of Y has a finite sub-cover.

> **Definition (Sequentially Compact).**
> X is sequentially compact if every sequence in X has a subsequence that converges in X.

> **Theorem (Compact ⟺ Sequentially Compact for metric spaces).**
> In a metric space, compactness and sequential compactness are equivalent.

### Compact ⟹ Sequentially Compact (Theorem 2.39) — Proof

*Proof.* Suppose X is compact but NOT sequentially compact. Then ∃ a sequence (xₙ) with no convergent subsequence. So for every x ∈ X, ∃ εₓ > 0 such that B\_{εₓ}(x) contains only finitely many terms of (xₙ). The collection {B\_{εₓ}(x) : x ∈ X} is an open cover. By compactness, extract a finite subcover. But each set in a finite subcover contains only finitely many terms ⟹ the whole sequence has only finitely many terms, a contradiction. □

> **Compact ⟹ Closed and Bounded (Theorem 2.33 + Lemma 2.36).**
> If Y is compact in (X, d), then Y is closed in X and bounded.

> **Closed Subset of Compact is Compact (Proposition 2.32).**
> If X is compact and Y ⊆ X is closed, then Y is compact.

*Proof.* Let R be an open cover for Y. Since Y is closed, X\Y is open. So R ∪ {X\Y} is an open cover for X. X compact ⟹ finite subcover of R ∪ {X\Y} covering X. Remove X\Y ⟹ finite subcover of R covering Y. □

> **Heine-Borel Theorem (Theorem 2.37).**
> In (Rⁿ, d₂), X is compact ⟺ X is closed and bounded.
>
> ⚠ **WARNING: This is specific to Rⁿ with the Euclidean metric! It fails in general metric spaces.**

⚠ **Exam Tip:** The 2023 markers explicitly noted: "A very common mistake: it is not possible to apply Heine-Borel theorem here because the metric is not Euclidean." The same warning appears in 2024 and 2025. If you're in C[0,1], an arbitrary metric space, or using a non-Euclidean metric, you CANNOT use Heine-Borel.

> **Continuous Image of Compact is Compact (Theorem 2.42).**
> If f: X → Y is continuous and X is compact, then f(X) is compact.

### 2.7.1 Showing D = {f ∈ C[0,1] : d∞(f, 0) ≤ 1} is NOT compact

This is a perennial exam question (2024 Q3(d), 2025 Q3(c)). The strategy: find a sequence in D with no convergent subsequence.

**Method:** Consider fₙ(t) = tⁿ (or fₙ(t) = sin(2πnt)). Each fₙ ∈ D. The pointwise limit of tⁿ is the discontinuous function that is 0 on [0,1) and 1 at t = 1. Any convergent subsequence in d∞ must converge uniformly, but the uniform limit of continuous functions is continuous — contradiction.

⚠ **Exam Tip:** Always use: d∞ convergence = uniform convergence, and the fact that uniform limits of continuous functions are continuous.

---

## 2.8 Completeness

> **Definition (Cauchy Sequence).**
> (xₙ) is Cauchy in (X, d) if ∀ε > 0 ∃N such that ∀n, m ≥ N: d(xₙ, xₘ) < ε.

> **Definition (Complete Metric Space).**
> (X, d) is complete if every Cauchy sequence in X converges to a limit in X.

**Key facts:**
- (Rⁿ, d₂) is complete.
- (Q, d₁) is not complete.
- ((0,1], d₁) is not complete.
- (C[0,1], d∞) IS complete.
- (C[0,1], d₂) is NOT complete.

### 2.8.1 Showing a space is NOT complete

Find a Cauchy sequence that doesn't converge in the space.

For integral-based metrics like d(x,y) = |∫ₓʸ h(t)dt| on (0,1): the sequence aₙ = 1 − 1/n is Cauchy (bound d(aₘ, aₙ)) but converges to 1 ∉ (0,1).

For (R, d\*) with d\*(x,y) = |∫ₓʸ e^{−t²}dt|: the sequence n → ∞ is Cauchy (because ∫\_N^∞ e^{−t²} dt → 0) but doesn't converge in (R, d\*) since d\*(n, a) stays bounded away from 0 for any a. This appeared in 2025 Q2(c).

### 2.8.2 C[0,1] with d∞ is Complete (Theorem 2.49)

Convergence in d∞ is uniform convergence. The uniform limit of continuous functions is continuous (this is the key theorem from Analysis I). So if (fₙ) is Cauchy in d∞, it converges uniformly to some function, which is continuous, hence in C[0,1].

### 2.8.3 C[0,1] with d₂ is NOT Complete

Consider the step-approximation sequence φₙ that converges in d₂ to the step function ψ(t) = −1 for t < 0 and ψ(t) = 1 for t ≥ 0. This is Cauchy in d₂ but the "limit" is not continuous, so the sequence doesn't converge in (C[−1,1], d₂).

---

## 2.9 Arzelà–Ascoli Theorem

> **Definition (Uniformly Bounded and Equicontinuous).**
> A family F ⊆ C[a,b] is **uniformly bounded** if ∃M > 0 such that |f(t)| ≤ M for all f ∈ F and t ∈ [a,b].
>
> F is **equicontinuous** if ∀ε > 0 ∃δ > 0 such that |x − y| < δ ⟹ |f(x) − f(y)| < ε for ALL f ∈ F.

> **Arzelà–Ascoli Theorem.**
> If F ⊆ C[a,b] is uniformly bounded and equicontinuous, then every sequence in F has a subsequence that converges uniformly (i.e., in d∞).

This appeared in 2023 Q3(d). The strategy: show the family is uniformly bounded (e.g., |f(x)| ≤ M for all f and x), then show equicontinuity (often via a Lipschitz condition |f(x) − f(y)| ≤ K|x − y|).

⚠ **Exam Tip:** 2024 markers: "Some students tried to apply Arzelà-Ascoli in part (d), but many sequences in C[0,1] are not equicontinuous." Only use A-A when you can verify BOTH conditions for the specific family.

---

## 2.10 Banach Fixed Point Theorem (Contraction Mapping)

> **Banach Fixed Point Theorem.**
> Let (X, d) be a **complete** metric space, and T: X → X be a **contraction**, i.e., ∃c ∈ [0,1) with d(T(x), T(y)) ≤ c·d(x,y) for all x, y. Then T has a unique fixed point x\* ∈ X, and for any x₀ ∈ X, the iterates Tⁿ(x₀) → x\*.

Both conditions matter:
1. X must be COMPLETE (see Problem Sheet 10, Ex 10.5: f(x) = x² on (0, 1/3) is a contraction with no fixed point because the space is not complete).
2. The Lipschitz constant must be strictly less than 1 (see Ex 10.6: d(f(x), f(y)) ≤ d(x,y) is not enough).

---
---

# Part 3: Key Problem Sheet Questions

**The following exercises are particularly worth revising, either because closely related questions have appeared on exams, or because they develop exam-critical techniques.**

## Differentiation (Sheets 1–5)

**Sheet 1, Exercise 1.1:** Cauchy-Schwarz, triangle inequality, reverse triangle inequality. Foundational — the reverse triangle inequality was directly tested in 2025 Q1(b).

**Sheet 4, Exercise 4.1:** Differentiability of quadratic forms Q(x) = xAxᵀ. This exact question appeared on the 2023 exam Q1(b).

**Sheet 4, Exercise 4.2:** A function where mixed partials exist but D₁D₂f ≠ D₂D₁f. Great for understanding why continuity of second partials is needed in Schwarz's theorem.

**Sheet 4, Exercise 4.5:** Applying the Inverse Function Theorem to f(x,y) = (x+y−xy, x²). Standard exam-style IFT question.

**Sheet 5 (IFT questions):** The exercises on the Implicit Function Theorem — computing derivatives of implicitly defined functions. Compare directly with 2022 Q1(c)–(d).

## Metric Spaces (Sheets 6–8)

**Sheet 6, Exercises on metric verification:** Practice with d(x,y) = |x³−y³|, the discrete metric, and the sup metric. These build the skills needed for the integral-based metrics that appear every year.

**Sheet 7, Open/Closed sets and subspace topology:** Understanding how openness/closedness depends on the ambient space and metric. Critical for the C[0,1] questions.

**Sheet 8, Closure/Interior/Boundary exercises:** The 2022 exam Q2(b) directly tests finding closure, interior, and boundary of R² minus a set.

**Sheet 8 Unseen:** Work through these — they test metric space concepts in less familiar settings, exactly as the exam does.

## Connectedness and Compactness (Sheets 9–10)

**Sheet 9, Connectedness exercises:** Showing sets are connected via path-connectedness (the standard exam technique). Also: disconnectedness proofs using continuous surjections to {0,1}.

**Sheet 9, Exercises on C[0,1]:** Showing D[0,1] is not open, not closed, but is connected — this was 2023 Q3 verbatim.

**Sheet 10, Completeness exercises:** Cauchy sequences, the Banach fixed point theorem. Exercise 10.4 (convergence of an iterative sequence to a root) is a typical fixed-point application.

**Sheet 10, Exercise 10.5 and 10.6:** Understanding why the Banach theorem requires completeness AND strict contraction.

## Compactness and Arzelà–Ascoli (Sheet 11)

**Sheet 11, Compactness exercises:** Showing compactness and non-compactness in various settings. These directly prepare you for Q3 of the exam.

**Sheet 11 Unseen problems:** These are excellent exam practice — they combine compactness, completeness, and function space arguments in challenging ways.

---

*Good luck with your revision!*