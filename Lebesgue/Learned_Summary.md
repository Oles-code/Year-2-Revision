# MATH50006 — Exam Preparation Summary

## Part A: Topics You've Found Difficult / Asked Clarifications On

Based on our conversations across this project, these are the areas where you've asked for extra help, gotten confused, or needed conceptual clarification.

### 1. Radon-Nikodym Notation & Integrating w.r.t. Different Measures

You asked about the meaning of `dμ_k/dλ` notation (the Radon-Nikodym derivative) and were confused about why we write `dλ` vs `dλ(x)`, and how to think about integrating with respect to different measures. The core confusion was around **when you need simple functions vs. when you can just use a density**.

**Key takeaway:** If μ = fλ (i.e. μ has density f w.r.t. Lebesgue), then ∫g dμ = ∫g·f dλ. This is the workhorse for all concrete computations. The simple-function machinery stays in the background and is only needed for theoretical proofs (the 4-stage "indicators → simple → non-negative → general" argument).

### 2. Absolute Continuity of Measures (ν ≪ μ)

You asked about this multiple times, including on the 2024 exam and lecture notes. The sticking points were:
- **Why you can't always use R-N backwards:** If μ is "smaller" (blind to part of the space), no density g can recover λ's mass on the invisible region.
- **The standard trap:** μ = 1_{[0,1]}dλ and ν = ½·1_{[0,2]}dλ. Then μ ≪ ν is TRUE but ν ≪ μ is FALSE — because μ([1,2]) = 0 while ν([1,2]) = ½. **This exact pair appeared in both 2023 Q1c and 2025 Q1d.**

### 3. Convergence in Measure vs. a.e. Convergence

You needed conceptual clarification on:
- The trick |sin(u)| ≤ |u| for showing convergence in measure (used in 2023 Q2a, 2024 Q1c)
- The "subsequence of a subsequence" technique: convergence in measure → extract subsequence converging a.e. → apply DCT → get contradiction. This was central to 2024 Q2b (the Ψ question).
- The exam comment confirms this was a common student error: "many only show pointwise convergence, which is not enough."

### 4. Markov's Inequality Applied to Lp

You were confused about why the Markov bound for Lp has ‖f‖^p_p (with the p-th power) in the numerator rather than ‖f‖_p. The answer: Markov applied to g = |f|^p at threshold t^p gives μ({|f| ≥ t}) ≤ ‖f‖^p_p / t^p.

### 5. The Cantor Set / Devil's Staircase Construction

You needed clarification on several parts of 2024 Q4:
- **Why C is uncountable:** bijection with {0,1}^ℕ (binary choice at each stage)
- **Why |f_{k+1} - f_k| ≤ 2^{-k}:** Both f_k and f_{k+1} agree at interval endpoints and are monotone, with f_k(b) - f_k(a) = 2^{-k}
- **Why f = f_k in a neighbourhood of x ∉ C:** x sits in a removed open interval from some stage k; no subsequent stage touches it, so all later f_j agree there
- **What dμ_k/dλ means:** μ_k is Lebesgue restricted to C_k and renormalised to total mass 1

### 6. Measurability: Checking via Generating Families

You asked about the remark that f is measurable iff f⁻¹(A₀) ∈ A for all A₀ in a generating family E₀. The proof trick: define G = {A₀ ⊆ Y : f⁻¹(A₀) ∈ A}, show G is a σ-algebra containing E₀, hence G ⊇ σ(E₀) = A₀.

### 7. Restriction of a Measure

You were confused why μ|_A(B) = μ(B) rather than μ(B∩A). The resolution: every B in the domain F|_A is already a subset of A (by construction as A∩C for C ∈ F), so B∩A = B automatically.

### 8. Mollifiers and Change of Variables (Jacobians)

On Sheet 7 Q4a, you needed clarification on why ∫φ_ε = 1 (the ε^{-d} prefactor cancels with the Jacobian ε^d from the substitution u = x/ε).

---

## Part B: Key Theorems (Precise Statements from Notes)

### Monotone Convergence Theorem (Thm 2.17(ii))
If f_n : X → [0,∞] measurable with f_n ↑ f (i.e. f_n(x) ≤ f_{n+1}(x) for all x and lim f_n(x) = f(x) for all x), then
$$\lim_{n\to\infty}\int f_n\,d\mu = \int f\,d\mu.$$

### Fatou's Lemma (Thm 2.22)
Let f_n : X → [0,∞] measurable. Then
$$\int \liminf_{n\to\infty} f_n\,d\mu \leq \liminf_{n\to\infty}\int f_n\,d\mu.$$

### Dominated Convergence Theorem (Thm 2.24)
Let g ∈ L¹(μ), and f, f_n measurable with |f_n(x)| ≤ g(x) for all x, and f_n → f μ-a.e. Then
$$\left|\int f_n\,d\mu - \int f\,d\mu\right| \leq \int|f_n - f|\,d\mu \to 0.$$

### Egorov's Theorem (Thm 2.9)
Let λ(Ω) < ∞ and f_k → f a.e. on Ω. Then for all δ > 0, there exists F ⊂ Ω compact with λ(Ω\F) < δ and f_k → f **uniformly** on F.

### Lusin's Theorem (Thm 2.11)
If λ(Ω) < ∞ and f : Ω → ℝ is measurable, then for all δ > 0, there exists F ⊂ Ω compact with λ(Ω\F) < δ and f|_F is **continuous**.

### Vitali's Theorem (Thm 2.36)
Let μ(X) < ∞, f, f_n ∈ L¹(μ). Then the following are equivalent:
- (i) f_n → f in measure AND {f_n} has uniformly absolutely continuous integrals
- (ii) ∫|f_n − f| dμ → 0

### Absolute Continuity of the Integral (Prop 2.32)
Let f ∈ L¹(μ). Then ∀ε > 0, ∃δ > 0 such that μ(A) < δ ⟹ ∫_A |f| dμ < ε.

### Hölder's Inequality (Cor 2.42(i))
For conjugate exponents 1/p + 1/q = 1, f ∈ Lp(μ), g ∈ Lq(μ):
$$\|fg\|_1 \leq \|f\|_p\|g\|_q.$$

### Minkowski's Inequality (Cor 2.42(ii))
For 1 ≤ p ≤ ∞ and f, g ∈ Lp(μ):
$$\|f+g\|_p \leq \|f\|_p + \|g\|_p.$$

### Fubini's Theorem (Thm 3.5)
Let (X_i, A_i, μ_i) be σ-finite, i = 1,2. If f ≥ 0 or f ∈ L¹(μ₁ ⊗ μ₂), then
$$\int f\,d(\mu_1\otimes\mu_2) = \int\!\left(\int f(x_1,x_2)\,d\mu_2\right)d\mu_1 = \int\!\left(\int f(x_1,x_2)\,d\mu_1\right)d\mu_2.$$

### Hahn-Carathéodory Extension (Thm 1.13)
Given a pre-measure μ̃ on an algebra A, the outer measure construction yields a measure μ on (X, Σ) with A ⊂ Σ and μ|_A = μ̃. If μ̃ is σ-finite, the extension is **unique** (Thm 1.14).

### Convergence in Measure: Key Relations (Thm 2.28)
- ∫|f_n − f| dμ → 0 ⟹ f_n →^μ f (when μ(X) < ∞)
- f_n → f a.e. ⟹ f_n →^μ f (when μ(X) < ∞)
- f_n →^μ f ⟹ ∃ subsequence with f_{n_k} → f a.e.
- **Converse of the second is FALSE** (typewriter/sliding-bump counterexample)

### Lebesgue Differentiation Theorem (Thm 4.2)
For f ∈ L¹_loc(ℝⁿ): f(x) = lim_{r↓0} (1/λ(B(x,r))) ∫_{B(x,r)} f(y) dy, for λ-a.e. x.

### Radon-Nikodym Theorem
If ν ≪ μ (both σ-finite), then ν = fμ for some measurable f ≥ 0, i.e. ν(A) = ∫_A f dμ.

---

## Part C: Key Proofs and Repeated Ideas on Exams

### The "Big Three" Proof Techniques

**1. DCT Application (appears in EVERY paper, 12–16 marks)**
The checklist:
- Identify pointwise limit f
- Find dominating function g with |f_n| ≤ g and g ∈ L¹
- Verify g is integrable
- Common dominating functions: 1/(1+x²), √x/(1+x), 1/x², bounded functions on finite domains
- Common trick: split into cases x ≤ n and x > n; use |sin(y)/y| ≤ 1 (or ≤ C on [0,1]) for the first, |sin| ≤ 1 for the second

**2. Fubini Applications (appears 2022, 2023, 2025)**
- For f ≥ 0: Fubini always applies (both integrals may be ∞ but they agree)
- For f changing sign: MUST first show |f| ∈ L¹ (apply Fubini to |f| first, it's ≥ 0)
- The contrapositive: if two iterated integrals disagree, then f ∉ L¹
- Key computation: ∫₀^∞ e^{-xy} dy = 1/x (justify by MCT)

**3. The Layer Cake / Markov Chain (2024 Q3)**
- Markov: μ({|f| ≥ t}) ≤ ‖f‖^p_p / t^p
- Layer cake: ‖f‖^p_p = ∫₀^∞ p t^{p−1} μ(|f| ≥ t) dt (proved via Fubini)
- Tail condition (Tail)_q with q > p ⟹ f ∈ Lp (apply layer cake to truncations, then MCT)

### Proof Ideas That Repeat

**The "subsequence-of-subsequence" contradiction:**  
If you have convergence in measure but need a.e. convergence for DCT, suppose the conclusion fails → extract subsequence where it fails → extract sub-subsequence converging a.e. → DCT applies → contradiction. Used in 2024 Q2b, Vitali's proof.

**The "split into bounded + tail" argument:**  
Write |f| = min(|f|, N) + (|f| − min(|f|, N)). The bounded part is controlled on small sets; the tail is globally small by MCT. Used in: absolute continuity of the integral (Prop 2.32), several DCT applications.

**The "show it's a σ-algebra containing the generators" trick:**  
Define G = {sets with the desired property}, show G is a σ-algebra, show G contains a generating family E₀, conclude G ⊇ σ(E₀). Used in: measurability via generating families, Fubini proof (Lemma 3.6), Dynkin's π-λ theorem applications.

**Bounding with |sin(x)| ≤ |x|:**  
This is used repeatedly for convergence-in-measure problems (2022 Q1, 2023 Q2a, 2024 Q1c, 2025 Q1). The pattern: the set {|f_n| > ε} becomes empty for large n because the bound forces |x| > nε which exceeds the support.

**Lp membership via splitting at singularity/infinity:**  
For f(x) = |x|^{-α}: on (0,1), f ∈ Lp iff αp < 1; on (1,∞), f ∈ Lp iff αp > 1. Consequence: never in Lp(0,∞). Tested in Q1 of every single paper.

---

## Part D: Common and/or Difficult Exam Questions

### Tier 1: Appears in (Almost) Every Paper

**DCT Limit Computations (Category D — hardest marks)**

2022 Q3b: Determine lim_{n→∞} ∫₀^∞ [n sin(x/n)] / [x(1+x²)] dλ.  
*The dominating function requires splitting into x ≤ n and x > n, using that sin(y)/y extends to a bounded continuous function on [0,1]. Answer: π/2.*

2023 Q2b: Determine lim_{n→∞} ∫₀^∞ [√x / (1+nx³)] cos(nx³) dλ.  
*Dominate by √x/(1∨x³), which is L¹. Pointwise limit is 0. Answer: 0.*

2024 Q2a: Determine lim_{n→∞} ∫₀^∞ sin(eˣ)/(1+nx²) dx.  
*Dominate by 1/(1+x²). Answer: 0.*

2025 Q3: Show f_a(x,t) = sin(x)e^{-xt}1_{[0,a]}(x) ∈ L¹, use Fubini to get I(a) = ∫₀^a sin(x)/x dx as an integral involving 1/(1+t²), then take a→∞ using DCT. *The exam comments note common errors: bounding |f_a| with |sin x| ≤ 1 instead of |sin x| ≤ x, and not finding a proper dominating function for part (c).*

**True/False: Convergence Modes (Category A–B, ~16 marks per paper)**

Every paper gives a sequence f_n and asks about: (i) uniform convergence, (ii) convergence in measure, (iii) Lp convergence. The pattern is always:
- Uniform: FALSE if sup|f_n| doesn't → 0 (usually the "spike" doesn't shrink in height)
- Convergence in measure: usually TRUE; show λ({|f_n| > ε}) → 0
- Lp convergence: compute ∫|f_n|^p explicitly

**True/False: Absolute Continuity (appears 2023, 2024, 2025)**

The pair μ = 1_{[0,1]}dλ, ν = ½·1_{[0,2]}dλ appears in 2023 Q1c AND 2025 Q1d. The test: find a set A with μ(A)=0 and check if ν(A)=0.

**True/False: Lp Membership (appears every paper)**

Checking whether f(t) = t^{-α}·(extra factors) is in Lp. Split at problem points, use the αp < 1 / αp > 1 criterion.

### Tier 2: Appears in Multiple Papers (High Marks)

**Fubini with Sign-Changing Functions (2022 Q4a, 2023 Q3, 2025 Q3)**

2023 Q3: f(x,y) = (x−y)/(x+y)³ on E ⊂ ℝ².  
Part (a) E = [1,∞)×[1,2]: bound |f| by 1/(x+y)², show that's L¹ via Fubini, then compute.  
Part (b) E = [1,∞)²: compute both iterated integrals → get +½ and −½ → they disagree → f ∉ L¹.  
*This is the contrapositive trick and it's a very clean exam question.*

2022 Q4a: G(x,y) = e^{-xy}g(x), bound |g(x)| ≤ √x/(1+x). Use Fubini on |G| (non-negative), integrate e^{-xy} in y first to get 1/x, then split ∫|g|/x into [0,1] and [1,∞).

**Cantor Set (2024 Q4 — entire 20 marks)**

Part (a): C uncountable (bijection with {0,1}^ℕ). Part (b): λ(C) = 0. Part (c): Devil's staircase construction.  
*Exam comments say "several people did not recognise the Cantor set or were not able to properly show it is uncountable" and "the main points missing were showing the f_k's are Cauchy and converge uniformly, and correctly showing f is locally constant outside C."*

**Non-Separability of L∞ (2022 Q4b)**

The family {f_t = 1_{[0,t]} : 0 < t < 1} satisfies ‖f_t − f_s‖_{L∞} = 1 for s ≠ t. Hence any ball of radius < ½ contains at most one f_t → no countable dense subset exists.

### Tier 3: Harder / Less Frequent but High-Value

**Hahn-Carathéodory from a Specific Function (2025 Q2 — 20 marks, Category B+D)**

Given f left-continuous, non-decreasing, define μ([a,b)) = f(b)−f(a) and the outer measure μ*. Show μ*([a,b)) = μ([a,b)). Part (a) is easy (choose the trivial cover). Part (b) is the hard direction: use left-continuity to nudge [a,b) to [a,b'], enlarge each [a_j,b_j) to (a_j',b_j), extract finite subcover by compactness, use monotonicity.  
*Comments: "Some common mistakes involved assuming that μ is a measure. A fair number had the beautiful idea of applying Hahn-Carathéodory (not needed)."*

**Vitali + Convergence in Measure (2024 Q2b — 12 marks, Category C)**

Show Ψ(f_n − f) → 0 where Ψ(g) = ∫|g|/(1+|g|) dx, given only convergence in measure.  
*The key: contradiction → subsequence → sub-subsequence converging a.e. → DCT (bounded by 1) → contradiction.*

**Brezis-Lieb Style (2023 Q4 — 16 marks, Category C+D)**

Given f_n → f a.e. with sup ‖f_n‖_p < ∞, show ∫||f_n|^p − |f_n−f|^p − |f|^p| dμ → 0. Uses a given algebraic inequality, the hint to define G^ε_n, DCT for the majorant (1+C(ε))|f|^p, and an ε-first-then-n argument.

**Riemann Sum → Lebesgue Integral (2025 Q4b — 12 marks, Category D)**

Show (1/n)Σf(k/n) → ∫f dλ for bounded measurable f with λ(D_f) = 0. The idea: f_n(x) = Σf(k/n)1_{[(k-1)/n, k/n)} is a simple function; f_n → f a.e. (on continuity points); |f_n| ≤ C (bounded); DCT.  
*Comments: "Not many students noticed that the set C_ε in part 1 is open."*

**Lusin's Theorem Application (2023 Q1d(ii))**

Show there exists closed C' ⊂ (0,1) with h|_{C'} continuous and λ(C') > ¾. Either cite Lusin directly, or construct explicitly by removing small intervals around the rationals.

### Summary of Category D (Hardest 20%) Mark Allocation

| Year | Category D Questions | Marks |
|------|---------------------|-------|
| 2022 | Q3b (DCT computation) | 16 |
| 2023 | Q4b,c (Brezis-Lieb) | 16 |
| 2024 | Q1a(ii) (Cantor counterexample) + Q4c (Devil's staircase) | 16 |
| 2025 | Q2b (outer measure) + Q4b(ii) (Riemann→Lebesgue) | 16 |

The hardest marks consistently test: DCT with non-trivial dominating functions, Fubini on sign-changing integrands, uniform convergence arguments (like the Devil's staircase), and measure construction.