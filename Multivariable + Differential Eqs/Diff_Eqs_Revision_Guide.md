# MATH50004 — Differential Equations: Exam Revision Guide

## How to Use This Document

This guide focuses on the **last three questions of the exam** (Q4, Q5, Q6), which cover differential equations. Looking at papers 2022–2025, the structure is remarkably consistent:

- **Q4** — Existence/uniqueness theory: Picard–Lindelöf (global & local), Picard iterates, maximal solutions, global existence arguments
- **Q5** — Linear systems: matrix exponential $e^{At}$, stability of trivial equilibrium, Lyapunov exponents, exponential bounds
- **Q6** — Nonlinear systems: equilibria, linearised stability, Lyapunov functions, invariant regions, domain of attraction, Poincaré–Bendixson

The material below follows the order of the lecture notes, but prioritises what has actually been tested. I'll highlight which parts are "must know cold" vs "good to know" and show proofs for results that regularly appear as exam questions.

---

## 1. Differential Equations and Initial Value Problems

### 1.1 Basic Setup

A $d$-dimensional ODE is $\dot{x} = f(t,x)$ where $f : D \subset \mathbb{R} \times \mathbb{R}^d \to \mathbb{R}^d$. A **solution** on interval $I$ is a differentiable function $\lambda : I \to \mathbb{R}^d$ with $(t, \lambda(t)) \in D$ and $\dot\lambda(t) = f(t, \lambda(t))$ for all $t \in I$.

The equation is **autonomous** if $f$ doesn't depend on $t$, i.e. $\dot{x} = f(x)$.

**Conceptually:** think of $f(t,x)$ as a "velocity vector" assigned to every point $(t,x)$. A solution is a path $\lambda$ whose velocity at every moment matches the prescribed velocity at its current location.

### 1.2 Translation Invariance (Autonomous Case)

For autonomous systems, if $\lambda(t)$ is a solution, so is $\lambda(t - c)$ for any constant $c$. **Why it matters:** this justifies defining the *flow* $\varphi(t, x)$ as "where you land after time $t$ starting at $x$" without worrying about the starting time.

---

## 2. Picard Iterates and Reformulation as an Integral Equation

**Exam relevance: EVERY exam has an IVP question in Q4 that tests this.** 2025 Q4 asked to compute Picard iterates explicitly; 2022 Q4 asked the same. Understanding this is non-negotiable.

### 2.1 The Integral Reformulation (Proposition 2.1)

$\lambda$ solves the IVP $\dot{x} = f(t,x)$, $x(t_0) = x_0$ on interval $I \ni t_0$ **if and only if** $\lambda$ is continuous and
$$\lambda(t) = x_0 + \int_{t_0}^{t} f(s, \lambda(s)) \, ds \quad \text{for all } t \in I.$$

**Proof sketch:** $(\Rightarrow)$ integrate the ODE from $t_0$ to $t$. $(\Leftarrow)$ differentiate the integral equation using the fundamental theorem of calculus. Both directions use the initial condition $\lambda(t_0) = x_0$.

**Why this matters:** the unknown $\lambda$ appears on both sides of an integral equation, which sets us up for a *fixed-point iteration*.

### 2.2 Picard Iterates (Definition 2.2)

For the IVP $\dot{x} = f(t,x)$, $x(t_0) = x_0$, define:
- $\lambda_0(t) \equiv x_0$ (the constant initial function)
- $\lambda_{n+1}(t) := x_0 + \int_{t_0}^{t} f(s, \lambda_n(s)) \, ds$

**If $\{\lambda_n\}$ converges uniformly** to some $\lambda_\infty$, then $\lambda_\infty$ solves the IVP. Uniform convergence is needed to swap the limit with the integral.

### 2.3 Exam-Standard Computation

For $\dot{x} = ax$, $x(0) = x_0$:
$$\lambda_0 = x_0, \quad \lambda_1 = x_0 + atx_0, \quad \lambda_2 = x_0 + atx_0 + \tfrac{a^2t^2}{2}x_0, \quad \ldots$$
$$\lambda_n(t) = \sum_{k=0}^{n} \frac{(at)^k}{k!} x_0 \to x_0 e^{at}.$$

**This matrix-valued generalisation gives us $e^{At}$ — see §5.**

**2025 Q4(a)(ii) asked to show all Picard iterates are polynomials** for $\dot{x} = t^2 x^3$, $x(0) = 1$. The strategy: show by induction that if $\lambda_n$ is a polynomial, then so is $\lambda_{n+1}$ (because integrating a polynomial in $t$ gives a polynomial, and $t^2 \cdot (\text{polynomial})^3$ is a polynomial).

---

## 3. Lipschitz Continuity

### 3.1 Definitions (Definition 2.5, 2.12)

**Lipschitz continuous:** $\|f(x) - f(\bar{x})\| \leq K \|x - \bar{x}\|$ for all $x, \bar{x}$ (with a *fixed* $K$).

For ODEs $\dot{x} = f(t,x)$ on domain $D$:
- **Globally Lipschitz in $x$:** there exists $K > 0$ with $\|f(t,x) - f(t,y)\| \leq K\|x-y\|$ for *all* $(t,x), (t,y) \in D$.
- **Locally Lipschitz in $x$:** every $(t_0, x_0)$ has a neighbourhood $U$ where a Lipschitz constant exists (which may depend on $U$).

### 3.2 Checking Lipschitz Continuity (Exam Technique)

**Rule of thumb (from the Mean Value Inequality):** if $\partial f / \partial x$ is *bounded* on a region, then $f$ is Lipschitz there with constant equal to the bound on the derivative.

- $f(x) = x^2$ on $\mathbb{R}$: **not** Lipschitz (derivative $2x$ unbounded).
- $f(x) = x^2$ on $[0,1]$: Lipschitz with constant $K = 2$.
- $f(x) = \sqrt{x}$ on $[0,1]$: **not** Lipschitz (derivative blows up at $0$).

**2024 Q4(a)(i) asked:** show $\dot{x} = t^n/x$ is not globally Lipschitz in $x$ on $\mathbb{R} \times (0, \infty)$. The answer: $\partial f/\partial x = -t^n/x^2$, which blows up as $x \to 0^+$, so no global Lipschitz constant can exist.

### 3.3 Continuous Differentiability $\Rightarrow$ Locally Lipschitz (Proposition 2.14)

**Result:** If $f : D \to \mathbb{R}^d$ is $C^1$ (continuously differentiable), then $f$ is locally Lipschitz in $x$.

**Proof (worth knowing — short):** for any $(t_0, x_0) \in D$, pick a compact convex neighbourhood $U \subset D$. Since $f$ is $C^1$ and $U$ is compact, $\|\partial f / \partial x\|$ is bounded by some $K$ on $U$. By the Mean Value Inequality (Theorem 2.8), for any $(t,x), (t,y) \in U$:
$$\|f(t,x) - f(t,y)\| \leq \|\tfrac{\partial f}{\partial x}(t, \xi)\| \cdot \|x-y\| \leq K \|x-y\|.$$
$\square$

**Why this is useful on exams:** in most concrete IVP problems, the right-hand side is a nice polynomial or rational function — which is $C^1$ on its domain — so local Lipschitz continuity comes for free. That's enough for the **local** version of Picard–Lindelöf.

### 3.4 Mean Value Inequality (Theorem 2.8)

For $f : D \subset \mathbb{R}^n \to \mathbb{R}^m$ continuously differentiable, and $[x,y] \subset D$:
$$\|f(x) - f(y)\| \leq \|f'(\xi)\| \cdot \|x - y\|$$
for some $\xi \in [x,y]$. (In higher dimensions it's an *inequality*, not an equality like MVT in 1D.)

---

## 4. Picard–Lindelöf Theorem

**This is the #1 most tested concept in Q4 — every exam tests whether students can distinguish the global vs local versions and apply the correct one.**

### 4.1 Global Version (Theorem 2.11)

If $f : \mathbb{R} \times \mathbb{R}^d \to \mathbb{R}^d$ is continuous and satisfies a **global** Lipschitz condition in $x$:
$$\|f(t,x) - f(t,y)\| \leq K \|x - y\| \quad \text{for all } t \in \mathbb{R}, \, x, y \in \mathbb{R}^d$$

Then every IVP $\dot x = f(t,x)$, $x(t_0) = x_0$ has a **unique** solution on $[t_0 - h, t_0 + h]$ where $h = 1/(2K)$.

**Proof idea (you should know the structure, not every detail):**
1. Reformulate as a fixed-point problem $\lambda = P(\lambda)$ where $P(u)(t) := x_0 + \int_{t_0}^t f(s, u(s))\,ds$ on the Banach space $X = C^0([t_0-h, t_0+h], \mathbb{R}^d)$ with supremum norm.
2. Show $P$ is a contraction: using Lipschitz, $\|P(u_1) - P(u_2)\|_\infty \leq Kh \|u_1 - u_2\|_\infty = \tfrac{1}{2}\|u_1 - u_2\|_\infty$.
3. Apply Banach fixed-point theorem to get a unique fixed point, which solves the IVP.

**Key insight:** the constant $h = 1/(2K)$ makes $Kh = 1/2 < 1$, guaranteeing contraction.

### 4.2 Local Version (Theorem 2.13)

If $f : D \to \mathbb{R}^d$ is continuous and **locally Lipschitz** in $x$, then for every $(t_0, x_0) \in D$, there exists $h > 0$ such that the IVP has a unique solution on $[t_0 - h, t_0 + h]$.

**Quantitative version:** if $W_{\tau, \delta}(t_0, x_0) := [t_0 - \tau, t_0 + \tau] \times \overline{B_\delta(x_0)} \subset D$ with Lipschitz constant $K$ and $\|f\| \leq M$ on this set, then
$$h = \min\{\tau, \tfrac{1}{2K}, \tfrac{\delta}{M}\}.$$

The constraint $\delta/M$ ensures the solution doesn't leave the ball $B_\delta(x_0)$ during time $h$.

### 4.3 When to Use Each Version — The Decision Tree

This is the skill the exam is actually testing:

1. Is the right-hand side $f$ defined on all of $\mathbb{R} \times \mathbb{R}^d$? 
   - If **no** → global version doesn't apply (e.g. 2023 Q4 where $f(t,x) = e^t/x^2$ is undefined at $x=0$). Use local version if $f$ is at least $C^1$ on its domain.
2. Is there a **global** Lipschitz constant $K$ valid for all $x, y$?
   - **Autonomous linear** $\dot{x} = Ax$: yes, $K = \|A\|$ works. Use global version.
   - **Polynomial with degree $\geq 2$** like $\dot{x} = x^2$: no, derivative grows unboundedly. Use local version only.
   - **Bounded derivative** like $\dot{x} = \sin(x)$: yes, $|\cos(x)| \leq 1$.
3. If $f$ is $C^1$ on its domain → automatically locally Lipschitz → local version applies and gives a unique local solution.

**Exam trap:** students often claim the global version applies when only local applies. The 2024 examiner comments explicitly flagged this.

---

## 5. Maximal Solutions

**Exam relevance:** tested in 2023 Q4 and 2024 Q4 (both ask for the maximal solution via separation of variables and verification it cannot be extended).

### 5.1 Definition (Definition 2.16, Theorem 2.17)

The **maximal existence interval** $I_{\max}(t_0, x_0) = (I^-(t_0, x_0), I^+(t_0, x_0))$ is the largest open interval on which the unique solution $\lambda_{\max}$ exists.

**Boundary behaviour theorem (Theorem 2.17):** if $I^+$ is finite, then *either*:
- $\lambda_{\max}$ is unbounded as $t \to I^+$ (blow-up), OR
- $(t, \lambda_{\max}(t))$ approaches the boundary of $D$ as $t \to I^+$.

**Conceptual picture:** the solution can't just "stop" — it has to either run off to infinity, or hit the boundary of the domain.

### 5.2 The Separation-of-Variables Exam Template

This is the standard Q4(a) workflow you should have memorised:

**Step 1 — Apply local Picard–Lindelöf:** check the right-hand side is $C^1$ on its domain, conclude unique local solution exists.

**Step 2 — Separate variables:** from $\dot{x} = g(t) h(x)$, write $\frac{dx}{h(x)} = g(t)\,dt$ and integrate.

**Step 3 — Apply initial condition** to determine the constant of integration.

**Step 4 — Solve explicitly for $x(t)$** — you should get a formula with a $t$-dependent denominator that can blow up.

**Step 5 — Find the blow-up time** by setting the denominator to zero. This gives you $I^+$.

**Step 6 — Verify maximality:** show the solution blows up (i.e. $\|\lambda(t)\| \to \infty$) as $t \to I^+$. By Theorem 2.17, this means the interval cannot be extended.

**Example walk-through (in the style of a typical Q4):** consider $\dot{x} = tx^2$, $x(0) = 1$ on $\mathbb{R} \times \mathbb{R}$.
- $f(t,x) = tx^2$ is $C^1$ → locally Lipschitz → unique local solution exists.
- Separating: $-\tfrac{1}{x} = \tfrac{t^2}{2} + C$, apply $x(0) = 1 \Rightarrow C = -1$.
- Solve: $x(t) = \dfrac{2}{2 - t^2}$.
- Blow-up at $t = \sqrt{2}$, so $I_{\max} = (-\sqrt{2}, \sqrt{2})$.
- As $t \to \sqrt{2}^-$, $x(t) \to +\infty$, confirming maximality.

---

## 6. Global Existence Arguments

**Exam relevance:** 2022 Q4 and 2023 Q4 both had parts asking "prove $I_{\max} = \mathbb{R}$". 2025 Q4 too: "if $f$ is bounded and locally Lipschitz, show solutions exist globally."

### 6.1 The Key Principle

A local solution extends globally (to all of $\mathbb{R}$) **unless it blows up in finite time or approaches the boundary**. So to prove global existence, you need to rule out finite-time blow-up.

### 6.2 Standard Argument Template

Suppose $f$ is locally Lipschitz and we want to show $I^+(t_0, x_0) = \infty$.

1. Assume for contradiction that $I^+ = T < \infty$.
2. By Theorem 2.17, either $\lambda_{\max}(t)$ is unbounded as $t \to T^-$, or $(t, \lambda_{\max}(t))$ approaches $\partial D$.
3. Rule out each possibility using an a priori bound on $\lambda_{\max}$.

**2025 Q4(b)(ii) prototype:** if $f : \mathbb{R} \times \mathbb{R}^d \to \mathbb{R}^d$ is bounded ($\|f\| \leq M$) and locally Lipschitz in $x$, then every solution satisfies
$$\|\lambda(t)\| = \left\| x_0 + \int_{t_0}^t f(s, \lambda(s))\,ds \right\| \leq \|x_0\| + M|t - t_0|,$$
so $\lambda$ grows at most linearly and can never blow up in finite time. Combined with the fact that $D = \mathbb{R} \times \mathbb{R}^d$ has no boundary, we get $I_{\max} = \mathbb{R}$.

**Conceptual note for the 2025 Q4(b)(i) question** ("does bounded imply globally Lipschitz?"): **No**. Counterexample: $f(x) = \sin(x^2)$ is bounded but its derivative $2x\cos(x^2)$ is unbounded, so it's not globally Lipschitz. Boundedness of $f$ and Lipschitz continuity of $f$ are different properties.

---

## 7. General Solutions and Flows

### 7.1 Definitions

The **general solution** $\lambda(t, t_0, x_0)$ is the maximal solution as a function of all three variables.

For **autonomous** systems $\dot{x} = f(x)$, translation invariance means $\lambda(t, t_0, x_0) = \lambda(t - t_0, 0, x_0)$, so we define the **flow** $\varphi(t, x) := \lambda(t, 0, x)$.

### 7.2 Key Properties (Proposition 2.24)

For the flow $\varphi$ of an autonomous system:
- **Initial value:** $\varphi(0, x) = x$
- **Group property:** $\varphi(t, \varphi(s, x)) = \varphi(t+s, x)$ — flowing for time $s$ then $t$ is the same as flowing for time $t+s$
- **Inverse:** $\varphi(-t, \varphi(t, x)) = x$

**Conceptually:** the flow defines a one-parameter family of maps $\{\varphi(t, \cdot)\}_{t \in \mathbb{R}}$ forming a group under composition.

### 7.3 Three Types of Orbits

For any $x \in D$, the orbit $O(x) = \{\varphi(t, x) : t \in J_{\max}(x)\}$ is exactly one of:
1. **A single point (equilibrium):** $f(x) = 0$ so $\varphi(t, x) = x$ for all $t$.
2. **A closed curve (periodic orbit):** $\exists T > 0$ with $\varphi(T, x) = x$ but $f(x) \neq 0$.
3. **A non-closed curve:** $t \mapsto \varphi(t, x)$ is injective.

**Exam fact (Proposition 2.27):** in dimension $d = 1$, there are *no* periodic orbits. Solutions are monotone. (This explains 2024 Q4(b)(ii): $t \mapsto \sin t$ cannot solve an autonomous 1D ODE with $C^1$ RHS because $\sin t$ is periodic and not constant.)

---

## 8. Matrix Exponential

**Exam relevance: EVERY exam Q5 asks to compute $e^{At}$ for a specific $2\times 2$ matrix.** The examiner comments for 2022 flagged that many students produced non-real matrix exponentials, losing lots of marks. You need to be fluent with this.

### 8.1 Definition and Basic Properties

$$e^{At} := \sum_{k=0}^{\infty} \frac{t^k A^k}{k!}$$

This series converges for all $t \in \mathbb{R}$ and $A \in \mathbb{R}^{d \times d}$ (Proposition 3.2 — proof uses sub-multiplicativity $\|AB\| \leq \|A\|\|B\|$ and comparison with $e^{\|A\|}$).

**The flow of $\dot{x} = Ax$ is $\varphi(t, x) = e^{At} x$** (Theorem 3.3).

**Properties (Proposition 3.4) — know these cold:**
1. If $C = T^{-1}BT$, then $e^C = T^{-1} e^B T$. **(Essential for computing $e^{At}$ via Jordan form.)**
2. $e^{-B} = (e^B)^{-1}$
3. **If $BC = CB$:** $e^{B+C} = e^B e^C$. **The commuting condition is essential** — fails in general. Try $B = \begin{pmatrix} 0 & 1 \\ 0 & 0\end{pmatrix}$, $C = \begin{pmatrix} 0 & 0 \\ 1 & 0\end{pmatrix}$.
4. Block-diagonal: $e^{\text{diag}(B_1, \ldots, B_p)} = \text{diag}(e^{B_1}, \ldots, e^{B_p})$.

### 8.2 The 2D Computation Recipe

Given $A \in \mathbb{R}^{2\times 2}$, follow these steps:

**Step 1 — Eigenvalues.** Solve $\det(A - \lambda I) = 0$.

**Step 2 — Identify which of four cases you're in:**
- **(C1)** Two distinct real eigenvalues $a \neq b$: $J = \begin{pmatrix} a & 0 \\ 0 & b \end{pmatrix}$, giving $e^{Jt} = \begin{pmatrix} e^{at} & 0 \\ 0 & e^{bt} \end{pmatrix}$.
- **(C2)** Repeated real eigenvalue $a$, diagonalisable: $J = aI$, giving $e^{Jt} = e^{at} I$.
- **(C3)** Repeated real eigenvalue $a$, non-diagonalisable (defective): $J = \begin{pmatrix} a & 1 \\ 0 & a \end{pmatrix}$, giving $e^{Jt} = e^{at}\begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}$.
- **(C4)** Complex conjugate eigenvalues $a \pm ib$: $J = \begin{pmatrix} a & b \\ -b & a \end{pmatrix}$, giving $e^{Jt} = e^{at}\begin{pmatrix} \cos(bt) & \sin(bt) \\ -\sin(bt) & \cos(bt) \end{pmatrix}$.

**Step 3 — Find transformation matrix $T$** so that $T^{-1}AT = J$:
- Real eigenvalues: $T$ has eigenvectors (and generalised eigenvectors in case C3) as columns.
- Complex eigenvalues: if $u + iv$ is an eigenvector for $a + ib$, use $T = (u | v)$ (real and imaginary parts).

**Step 4 — Assemble:** $e^{At} = T e^{Jt} T^{-1}$.

### 8.3 Phase Portraits for 2D Linear Systems

| Case | Eigenvalues | Portrait |
|---|---|---|
| C1, $a < b < 0$ | two negative | Stable node ("knot with two tangents") |
| C1, $a < 0 < b$ | opposite signs | **Saddle** |
| C1, $0 < a < b$ | two positive | Unstable node |
| C2, $a < 0$ | double negative, diagonalisable | Stable star |
| C3, $a < 0$ | double negative, defective | Stable improper node (one tangent) |
| C4, $a < 0$ | complex with Re < 0 | **Stable focus/spiral** |
| C4, $a = 0$ | pure imaginary | **Centre** (closed orbits) |
| C4, $a > 0$ | complex with Re > 0 | Unstable focus |

```
Saddle (a < 0 < b)         Centre (a = 0, b ≠ 0)      Stable focus (a < 0)
    ↖       ↗                 ↶   ↷                    ⤸  ⤹
      ↖   ↗                 ↶       ↷                  spirals
  ────●────→               ●       ●                    into 0
      ↙   ↘                 ↷       ↶
    ↙       ↘                 ↷   ↶
```

**Sign of $b$ in case C4 determines rotation direction:** $b > 0$ → clockwise; $b < 0$ → anticlockwise.

### 8.4 The Big Picture (Higher Dimensions)

For a general $A \in \mathbb{R}^{d \times d}$ in Jordan normal form, entries of $e^{Jt}$ are of the form $p(t) e^{\rho t}$ where $\rho$ is a real part of an eigenvalue and $p(t)$ is a polynomial of degree at most (block size $- 1$).

**Key exponential estimate (Proposition 3.9):** if $\gamma > \max\{\text{Re}\,\rho : \rho \text{ eigenvalue of } A\}$, then there exists $K > 0$ with
$$\|e^{At}\| \leq K e^{\gamma t} \quad \text{for all } t \geq 0.$$

**Strict inequality caveat:** you can use $\gamma = \max\text{Re}\,\rho$ (equality) **only if** all eigenvalues on that maximum are *semi-simple* (i.e. algebraic multiplicity = geometric multiplicity, or equivalently, their Jordan blocks have size 1 for real, size $2 \times 2$ for complex). Otherwise, polynomial factors $t^n$ grow, and you need strict $\gamma >$.

---

## 9. Lyapunov Exponents and Exponential Bounds

**Exam relevance:** 2022 Q5 ($\|e^{At}\| \leq Ke^{-t}$ for which $a$?), 2024 Q5 (Lyapunov exponents of specific solutions, bound on $\|e^{At}x\|/t^{d-1}$).

### 9.1 Lyapunov Exponent

For a nonzero solution $\lambda(t)$:
$$\sigma_{\text{Lyap}}(\lambda) := \lim_{t \to \infty} \frac{\ln \|\lambda(t)\|}{t}$$
(when the limit exists). This measures the exponential growth/decay rate.

**For $\dot{x} = Ax$:** the Lyapunov exponent of a solution starting in $x \in E_j$ (the generalised eigenspace for eigenvalues with real part $s_j$) is exactly $s_j$. The space decomposes as
$$\mathbb{R}^d = E_1 \oplus \cdots \oplus E_q$$
where $s_1 < s_2 < \cdots < s_q$ are the distinct real parts.

### 9.2 How to Answer Bound Questions (2022 Q5 style)

**Question type:** "For which $a$ does there exist $K > 0$ with $\|e^{At}\| \leq Ke^{-t}$ for all $t \geq 0$?"

**Strategy:** this is asking for the maximum real part of eigenvalues to be $\leq -1$:
$$\|e^{At}\| \leq Ke^{\gamma t} \text{ needs } \gamma > \max\text{Re}\,\rho \text{ (or equality if all are semi-simple)}.$$
So for $\gamma = -1$: need $\max\text{Re}\,\rho < -1$ (strict), or $\leq -1$ if the ones achieving $-1$ are semi-simple.

For "$t \leq 0$": the condition flips — equivalent to asking about $e^{-At}$, i.e. $-A$ has max real part $\leq -1$, i.e. $A$ has *min* real part $\geq 1$.

### 9.3 The $t^{d-1}$ Bound (2024 Q5b)

**Claim:** if all eigenvalues of $A \in \mathbb{R}^{d\times d}$ have zero real part, then
$$\limsup_{t\to\infty} \frac{\|e^{At}x\|}{t^{d-1}} < \infty \quad \text{for all } x \in \mathbb{R}^d.$$

**Why:** entries of $e^{At}$ are of the form $p(t) \cdot (\cos(bt), \sin(bt))$ where $p(t)$ is a polynomial of degree $\leq d-1$ (since the largest Jordan block of a $d\times d$ matrix has size $\leq d$). Since $|\cos|, |\sin| \leq 1$, the whole thing grows at most like $t^{d-1}$.

---

## 10. Stability Notions

**Exam relevance:** Q5 and Q6 both routinely ask for definitions and which stability notion holds in a given example.

### 10.1 The Five Notions (Definition 4.1)

Let $x^*$ be an equilibrium ($f(x^*) = 0$):

| Term | Definition | Intuition |
|---|---|---|
| **Stable** | $\forall \varepsilon > 0, \exists \delta > 0$: $x \in B_\delta(x^*) \Rightarrow \varphi(t,x) \in B_\varepsilon(x^*)$ $\forall t \geq 0$ | "Start close, stay close" |
| **Attractive** | $\exists \delta > 0$: $\lim_{t\to\infty} \varphi(t,x) = x^*$ for all $x \in B_\delta(x^*)$ | "Nearby points eventually converge" |
| **Asymptotically stable** | Stable AND attractive | "Stay close AND converge" |
| **Exponentially stable** | $\exists \delta > 0$, $K \geq 1$, $\gamma < 0$: $\|\varphi(t,x) - x^*\| \leq K e^{\gamma t} \|x - x^*\|$ | "Converge at exponential rate" |
| **Repulsive** | $\exists \delta > 0$: $\lim_{t\to-\infty} \varphi(t,x) = x^*$ for all $x \in B_\delta(x^*)$ | "Time-reversed attractive" |

**Crucial fact: stability and attractivity are independent** (Example 4.3). You can have attractive-but-unstable equilibria (the classic example involves a limit cycle with a stagnation point).

Exponential stability $\Rightarrow$ asymptotic stability (stable + attractive), but the converse is false. Example: $\dot{x} = -x^3$ is asymptotically stable at $0$ but not exponentially (decay is polynomial, not exponential).

### 10.2 Stability of Linear Systems (Theorem 4.5)

For $\dot{x} = Ax$, the trivial equilibrium $0$ is:

**(i) Stable** iff both:
- $\text{Re}\,\rho \leq 0$ for all eigenvalues $\rho$ of $A$, AND
- Every eigenvalue $\rho$ with $\text{Re}\,\rho = 0$ is **semi-simple**.

**(ii) Exponentially stable** iff $\text{Re}\,\rho < 0$ for *all* eigenvalues.

**Why semi-simplicity matters:** a Jordan block $\begin{pmatrix} 0 & 1 \\ 0 & 0\end{pmatrix}$ gives $e^{Jt} = \begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}$, which has an entry that grows to infinity, so we'd lose stability despite having $\text{Re}\,\rho = 0$.

**Proof of $(\Leftarrow)$ direction of (i) — worth knowing:** if (a) and (b) hold, Proposition 3.9 gives $\|e^{At}\| \leq K$ for all $t \geq 0$. Then for $x \in B_\delta(0)$ with $\delta := \varepsilon/K$: $\|\varphi(t,x)\| = \|e^{At}x\| \leq K\|x\| < K\delta = \varepsilon$.

**2024 Q5(a) application:** for $A = \begin{pmatrix} 0 & -1 \\ -2 & -1 \end{pmatrix}$ with eigenvalues $-2$ and $1$ — since there's an eigenvalue with positive real part, the origin is **unstable** and not attractive.

### 10.3 The Attractivity Detail (2023 Q5c)

**Claim:** if $x^* = 0$ is attractive for $\dot{x} = Ax$, then $\lim_{t\to\infty} \varphi(t, x) = 0$ for **all** $x \in \mathbb{R}^d$ (not just nearby ones).

**Why:** attractivity gives $\lim_{t\to\infty}\varphi(t,x) = 0$ for $x \in B_\delta(0)$. For arbitrary $x$, write $x = \frac{\delta}{2\|x\|} \cdot \frac{2\|x\|}{\delta} x = c \cdot y$ where $y \in B_\delta(0)$. By **linearity** $\varphi(t, cy) = c \varphi(t, y) \to 0$.

This argument uses the linearity of $\dot{x} = Ax$ crucially and **would not work for nonlinear systems**.

---

## 11. Hyperbolicity and Linearised Stability

### 11.1 Hyperbolic Equilibrium (Definition 4.7)

An equilibrium $x^*$ of $\dot{x} = f(x)$ is **hyperbolic** if the matrix $f'(x^*)$ has no eigenvalue with zero real part.

**Why this matters:** the behaviour near hyperbolic equilibria is determined by the linearisation (Hartman–Grobman). Near non-hyperbolic equilibria, nonlinear terms can completely change the picture (see Example 4.8: $\dot{x} = \pm x^2, \pm x^3$ all linearise to $\dot{x} = 0$ but have very different stability).

### 11.2 Linearised Stability Theorem (Theorem 4.10)

**Statement:** if $x^*$ is an equilibrium of $\dot{x} = f(x)$ with $f \in C^1$, and all eigenvalues of $f'(x^*)$ have **negative** real parts, then $x^*$ is exponentially stable.

**By time reversal:** if all eigenvalues have **positive** real parts, $x^*$ is exponentially *unstable* (repulsive).

**Why this is so useful:** to check stability of a nonlinear equilibrium, just compute eigenvalues of the Jacobian. No solving the ODE required.

**Proof sketch (too long for full exam but key idea is testable):**
1. Write $\dot{x} = Ax + r(x)$ where $A = f'(0)$ and $r(x) = o(\|x\|)$.
2. Use variation of constants: $\varphi(t, x_0) = e^{At}x_0 + \int_0^t e^{A(t-s)} r(\varphi(s, x_0))\,ds$.
3. Estimate: $\|e^{At}\| \leq Ke^{\gamma t}$ with $\gamma < 0$, and $\|r(x)\| \leq M\|x\|$ for $\|x\|$ small.
4. Get an implicit inequality on $u(t) := e^{-\gamma t}\|\varphi(t, x_0)\|$: $u(t) \leq K\|x_0\| + KM\int_0^t u(s)\,ds$.
5. Apply **Gronwall's Lemma** to get $u(t) \leq K\|x_0\| e^{KMt}$, so $\|\varphi(t, x_0)\| \leq K\|x_0\|e^{(KM + \gamma)t}$.
6. Choose $M < -\gamma/K$ so $KM + \gamma < 0$, giving exponential decay.

### 11.3 Gronwall's Lemma (Lemma 4.9) — Must Know

**Statement:** if $u : [a,b] \to \mathbb{R}$ is continuous with
$$0 \leq u(t) \leq c + d \int_a^t u(s)\,ds,$$
then $u(t) \leq c e^{d(t-a)}$.

**Intuition:** turns an implicit exponential-looking bound (with $u$ on both sides) into an explicit one.

**Why it's useful:** comes up everywhere in ODE theory — linearised stability proofs, comparison arguments, estimates of perturbations. This is one of the most frequently applied tools in this course.

---

## 12. Invariant Sets and Invariance Techniques

**Exam relevance:** *every* Q6 involves checking positive invariance of a region, either to bound the solution or to apply Poincaré–Bendixson. 2022 Q6(c), 2024 Q6(b), 2025 Q6(d) all test this.

### 12.1 Definitions (Definition 4.15)

- **Positively invariant:** $x \in M \Rightarrow O^+(x) \subset M$ (orbits can't escape forward in time).
- **Negatively invariant:** orbits can't escape backward in time.
- **Invariant:** both.

### 12.2 The Nullcline Method (Exam Standard)

To check whether a region $M$ bounded by curves is positively invariant:

1. **Identify the boundary of $M$.** Usually pieces of nullclines or coordinate axes.
2. **On each boundary piece, check the direction of the vector field** $(\dot x, \dot y)$:
   - Does it point into $M$, out of $M$, or along $M$?
   - For positive invariance, it must point inward (or along the boundary).
3. **Special case — axes and nullclines:** if $\dot{x} = 0$ on some curve, the vector field is vertical there. You only need to check the sign of $\dot{y}$ to determine whether trajectories move along or cross.
4. **Invariant axes:** if $\dot{x}$ has factor $x$, then $\{x = 0\}$ is invariant (positively and negatively). Trajectories can't cross.

**Nullcline picture:**
```
     ẏ = 0 curve
        |     
        |  → vector field direction
   ─────●───── ẋ = 0 curve
        |  ↑   equilibrium at intersection
        |
```

### 12.3 Finding Equilibria

Set $\dot{x} = 0$ AND $\dot{y} = 0$ simultaneously. Typical factored form:
- $\dot{x} = x \cdot g(x,y) = 0$ gives $x = 0$ OR $g = 0$.
- Solve each branch's intersection with $\dot{y} = 0$.

**Trap:** don't forget branches. The 2022 Q6 required showing *exactly one* equilibrium — requires ruling out other candidates.

### 12.4 Linearisation at an Equilibrium

The Jacobian at $(x^*, y^*)$:
$$f'(x^*, y^*) = \begin{pmatrix} \partial \dot{x}/\partial x & \partial \dot{x}/\partial y \\ \partial \dot{y}/\partial x & \partial \dot{y}/\partial y \end{pmatrix} \bigg|_{(x^*, y^*)}.$$

**Exam technique:** you don't always need eigenvalues individually:
- If $\text{tr} < 0$ and $\det > 0$: both eigenvalues have negative real part → **exponentially stable**.
- If $\det < 0$: eigenvalues have opposite signs → **saddle** (hyperbolic but unstable, has 1D stable and unstable manifolds).
- If $\text{tr} > 0$ and $\det > 0$: both have positive real part → **repulsive/unstable**.

---

## 13. Lyapunov Functions

**Exam relevance:** featured in 2023 Q6, 2024 Q6, and implicitly in 2025 Q6(d). Must know the direct methods.

### 13.1 Orbital Derivative (Definition 4.22)

For a $C^1$ function $V : D \to \mathbb{R}$, the **orbital derivative** is
$$\dot{V}(x) := V'(x) \cdot f(x) = \sum_{i=1}^d \frac{\partial V}{\partial x_i}(x) f_i(x).$$

This is the time-derivative of $V$ along solutions: $\frac{d}{dt} V(\varphi(t, x)) = \dot{V}(\varphi(t, x))$ (chain rule).

**Key mental model:** if $V$ is "energy," then $\dot V$ is the rate of energy change along trajectories.

### 13.2 Lyapunov Function (Definition 4.24)

$V : D \to \mathbb{R}$ is a **Lyapunov function** if $\dot{V}(x) \leq 0$ for all $x \in D$.

**Consequence:** $V(\varphi(t, x)) \leq V(x)$ for all $t \geq 0$ — $V$ is non-increasing along orbits. So **sublevel sets $\{V \leq c\}$ are positively invariant** (Proposition 4.25).

### 13.3 Direct Method for Stability (Theorem 4.26)

If there's a Lyapunov function $V$ with:
- $V(x^*) = 0$
- $V(x) > 0$ for $x \neq x^*$

then $x^*$ is **stable**.

**Geometric intuition:** if orbits can't increase $V$, and $V$ has a strict local minimum at $x^*$, then orbits near $x^*$ are trapped in progressively smaller sublevel sets.

### 13.4 Direct Method for Asymptotic Stability (Theorem 4.31)

Strengthen: if additionally $\dot{V}(x) < 0$ for all $x \neq x^*$ (in a neighbourhood), then $x^*$ is **asymptotically stable**.

### 13.5 La Salle's Invariance Principle (Theorem 4.28, Corollary 4.30)

**Statement:** for any $x$, the omega limit set satisfies
$$\omega(x) \subseteq \{y \in D : \dot{V}(y) = 0\}.$$

**Stronger version:** $\omega(x) \subseteq \text{largest invariant subset of } \{\dot{V} = 0\}$.

**Why useful:** when the strict Lyapunov condition fails ($\dot{V} < 0$ only on some set but $= 0$ elsewhere), La Salle still tells you where trajectories asymptotically go.

**Typical exam trick (2024 Q6):** the Lyapunov function gives $\dot V \leq 0$ (non-strict), so you can't directly conclude asymptotic stability. But combined with the observation that "$\dot V = 0$ only on a set that cannot contain a full orbit," you conclude that the trajectory must eventually leave that set, and thus $V$ must strictly decrease, converging to the equilibrium.

### 13.6 Domain of Attraction via Sublevel Sets (Corollary 4.33)

If $V$ is a Lyapunov function with strict minimum at $x^*$, and $S_c = \{V \leq c\}$ is a compact subset of $D$, then $S_c \subseteq W^s(x^*)$ (domain of attraction).

**Used in 2025 Q6(d):** to show a specific point $(1/2, 1)$ is in the domain of attraction, find a positively invariant region containing it that's bounded away from other equilibria.

### 13.7 How to Find a Lyapunov Function (Exam Strategy)

- **Try $V(x,y) = x^2 + y^2$** or weighted quadratics $V = ax^2 + by^2$ first — works for many systems.
- **For systems with logarithmic features:** try $V(x,y) = x - \ln x + y - \ln y$ or similar (2024 Q6b hint explicitly uses this for a Lotka–Volterra-style system).
- **Compute $\dot V$** and try to show it has the sign you want.
- **For instability (Exercise 37 in PS8):** find a function with $\dot V > 0$ near the equilibrium instead — reverses the arguments.

### 13.8 Ruling Out Periodic Orbits via Lyapunov Functions (2024 Q6b)

**Key trick:** if $V$ is strictly decreasing along orbits in some region ($\dot V < 0$ whenever $x \neq x^*$), **there can be no periodic orbits** in that region. Why? On a periodic orbit, $V$ would have to return to its starting value after one period, contradicting strict monotonicity.

---

## 14. Omega and Alpha Limit Sets

### 14.1 Definitions (Definition 4.17)

- $\omega(x) = \{y : \exists t_n \to \infty \text{ with } \varphi(t_n, x) \to y\}$ — accumulation points as $t \to \infty$.
- $\alpha(x) = \{y : \exists t_n \to -\infty \text{ with } \varphi(t_n, x) \to y\}$ — accumulation points as $t \to -\infty$.

### 14.2 Alternative Characterisation (Proposition 4.19)

$$\omega(x) = \bigcap_{t \geq 0} \overline{O^+(\varphi(t, x))}$$

**Conceptually:** $\omega(x)$ is "what's left after throwing away the transient parts" — you're looking at the tail of the trajectory.

### 14.3 Properties (Proposition 4.21)

- $\omega(x)$ is always invariant.
- If $O^+(x)$ is bounded with $\overline{O^+(x)} \subset D$, then $\omega(x)$ is non-empty and compact.

### 14.4 Common Exam Scenarios

**Scenario 1:** trajectory converges to an equilibrium $x^*$. Then $\omega(x) = \{x^*\}$.

**Scenario 2:** trajectory approaches a periodic orbit $\Gamma$. Then $\omega(x) = \Gamma$.

**Scenario 3 (2022 Q5e):** for a linear system with eigenvalues in the open right half-plane, orbits escape to infinity, so $\omega(x) = \emptyset$ for $x \neq 0$.

**Scenario 4 (2024 Q6):** if a Lyapunov function shows $\omega(x) \subseteq \{y : \dot V(y) = 0\}$ which contains only $x^*$, then $\omega(x) = \{x^*\}$, giving convergence.

---

## 15. Poincaré–Bendixson Theorem

**Exam relevance:** essential for Q6 in 2D systems. Directly tested in 2022 Q6 (show periodic orbit exists).

### 15.1 Statement (Theorem 4.34)

**For a 2D autonomous system** $\dot{x} = f(x)$ with $f \in C^1$: if $O^+(x)$ lies in a compact set $K$ containing only finitely many equilibria, then $\omega(x)$ is one of:
1. A single equilibrium.
2. A periodic orbit.
3. A union of equilibria and homoclinic/heteroclinic orbits connecting them.

**Why 2D matters:** in 3D or higher, you can get chaos (e.g. Lorenz). The Jordan curve theorem makes 2D special.

### 15.2 Standard Corollary (Finding Periodic Orbits)

**Exam template — "show a periodic orbit exists":**

1. **Find a positively invariant compact region $M$** (typically a rectangle or region bounded by nullclines).
2. **Show there are no equilibria in $M$**, or show any equilibria in $M$ are repulsive (e.g. sources, not sinks).
3. **Pick any $x \in M$.** Then $O^+(x)$ stays in compact $M$, so by Poincaré–Bendixson, $\omega(x)$ is non-empty.
4. **Rule out options (1) and (3)** in the theorem: since there are no stable equilibria (or no equilibria at all) in $M$, $\omega(x)$ must be a periodic orbit.

**The crucial check in exam questions:** that the equilibrium inside the region is **repulsive** (unstable with all eigenvalues having positive real part). Then trajectories can't converge to it, so must spiral into a periodic orbit.

### 15.3 Ruling Out Periodic Orbits (2024 Q6b)

Three common techniques:
1. **Bendixson–Dulac criterion (not explicitly in notes but implied):** if $\nabla \cdot f \neq 0$ throughout a simply-connected region, no periodic orbits there.
2. **Strict Lyapunov function:** if $\dot V < 0$ everywhere except at equilibrium, no periodic orbit possible (see §13.8).
3. **Invariant line/curve argument:** if the only invariant sets are unions of equilibria, and none form closed curves.

---

## 16. Variation of Constants Formula

**Exam relevance:** comes up for inhomogeneous problems (Problem Sheet 6) and in the proof of Theorem 4.10.

### 16.1 Statement (Proposition 3.10)

For $\dot{x} = Ax + g(t)$:
$$\lambda(t, t_0, x_0) = e^{A(t-t_0)} x_0 + \int_{t_0}^t e^{A(t-s)} g(s)\,ds.$$

**Intuition:** first term is the homogeneous solution; second term accumulates the "forcing" $g$, propagating each pulse through time via $e^{A(t-s)}$.

### 16.2 1D Scalar Version

For $\dot{x} = a(t) x + g(t)$:
$$\lambda(t, t_0, x_0) = e^{\int_{t_0}^t a(s)\,ds} x_0 + \int_{t_0}^t e^{\int_s^t a(\tau)\,d\tau} g(s)\,ds.$$

---

## 17. Putting It Together: Exam Structure Summary

| Exam Q | Topic | Sub-topics repeatedly tested |
|---|---|---|
| **Q4** | Existence/uniqueness | Picard–Lindelöf (which version applies), Picard iterates explicit computation, separation of variables for maximal solution, blow-up verification, global existence via a priori bound |
| **Q5** | Linear systems | $e^{At}$ computation, stability of trivial equilibrium, eigenvalue-based bounds $\|e^{At}\| \leq Ke^{\gamma t}$, periodic orbits in linear systems (answer: only when pure imaginary eigenvalues), Lyapunov exponents |
| **Q6** | Nonlinear systems | Find equilibria, linearise and classify, check positive invariance of a region, build Lyapunov functions, domain of attraction arguments, Poincaré–Bendixson for periodic orbits |

### What the Examiners Flag as Common Errors

- Confusing global and local Lipschitz (applying wrong Picard–Lindelöf version)
- Getting a complex (not real) matrix exponential in Q5 — check your work at the end, $e^{At}$ must be real for real $A$
- Not verifying *maximality* of a solution (just computing it isn't enough — you must show it can't be extended, typically via blow-up or boundary-approach)
- Rushing positive invariance checks — you need to check the vector field on *every* piece of the boundary
- Writing down eigenvalues without computing trace and determinant first — 2x2 shortcut saves time
- Forgetting that for 1D autonomous systems, there are no periodic orbits — a useful True/False trick

---

## 18. Key Problem Sheet Exercises for Revision

These are the exercises most directly relevant to exam-style questions:

### Problem Sheets 1–3: Picard iterates, Lipschitz, existence/uniqueness
- **PS1 Exercise 1** (Variation of constants) — elementary technique that underlies Chapter 3.
- **PS1 Exercise 3** (Picard iterates) — direct prep for 2022 Q4, 2025 Q4(a).
- **PS1 Exercise 4** (Comparison with differential inequality) — gives the Gronwall-style argument used later.
- **PS2 Exercise 7** (Lipschitz continuity) — foundational.
- **PS2 Exercise 8** (Global existence of solutions) — prep for 2023 Q4(b), 2025 Q4(b).
- **PS3 Exercise 11** (Unique/non-unique solutions) — trains pattern recognition for when uniqueness fails.
- **PS3 Exercise 12** (Maximal solutions) — *the* practice for Q4 separation-of-variables questions.
- **PS3 Exercise 14** (Criterion for non-global existence) — relates blow-up to unboundedness.

### Problem Sheet 4: Flows and orbits
- **PS4 Exercise 16** (Monotone/constant solutions of autonomous ODEs) — important for phase-plane arguments, 1D no-periodic-orbits fact.
- **PS4 Exercise 17** (Half-orbits) — conceptually important for invariance.
- **PS4 Exercise 19** (Variational equation) — foundational for perturbation analysis.

### Problem Sheet 5: Matrix exponential and linear systems
- **PS5 Exercise 21** (Properties of $e^{At}$) — identical to Proposition 3.4, must know cold.
- **PS5 Exercise 22** (Transformation of phase portraits) — shows how Jordan form results transfer to general systems.
- **PS5 Exercise 23** (Computation of matrix exponentials) — *the* direct practice for Q5(a). Includes the technique of writing $A = D + P$ for defective cases.
- **PS5 Exercise 24** (Bounded solutions of linear systems) — tests understanding of stability.

### Problem Sheet 6: Stability basics
- **PS6 Exercise 26** (Inhomogeneous linear systems) — variation of constants practice.
- **PS6 Exercise 27** (Stability for 1D ODEs) — trains distinction between attractivity and stability, see Example 4.8.
- **PS6 Exercise 29** (Exponential stability of linear systems) — proves the first part of Theorem 4.5 yourself.

### Problem Sheet 7: Phase-plane analysis
- **PS7 Exercise 31** (First steps in phase plane analysis) — foundational for Q6 in exams.
- **PS7 Exercise 32** (Stability under linear transformations) — conceptual.
- **PS7 Exercise 33** (Openness of domain of attraction) — challenging, justifies claim used in 2025 Q6(d).
- **PS7 Exercise 34** (Boundary of an invariant set) — technical but useful.

### Problem Sheet 8: Lyapunov and Poincaré–Bendixson
- **PS8 Exercise 37** (Instability via positive orbital derivative) — time-reversed Lyapunov argument.
- **PS8 Exercises 38, 39** (Phase plane analysis I, II) — full-scale Q6-style practice.
- **PS8 Exercise 40, 41** (Lyapunov direct method I, II) — *the* direct prep for 2024 Q6(a), 2023 Q6.
- **PS8 Exercise 42** (Existence of periodic orbit) — *the* direct prep for 2022 Q6(d).
- **PS8 Exercise 43** (Omega limit set with equilibria and heteroclinics) — prep for Poincaré–Bendixson type (iii).

---

## Appendix A: Quick Reference — When to Use What

| Question asks... | Use this |
|---|---|
| "Does a unique solution exist locally?" | Check $f$ is $C^1$ or locally Lipschitz → Local Picard–Lindelöf |
| "Does a unique solution exist globally?" | Need global Lipschitz (or bounded + locally Lipschitz). |
| "Compute the maximal solution" | Separation of variables (usually). Verify blow-up. |
| "Show $I_{\max} = \mathbb{R}$" | Get a priori bound on $\|\lambda(t)\|$ that rules out blow-up. Use Gronwall if needed. |
| "Compute $e^{At}$" | Eigenvalues → identify case → $e^{Jt}$ → conjugate back with $T$. |
| "$\|e^{At}\| \leq Ke^{\gamma t}$ for which $\gamma$?" | Max real part of eigenvalues (strict if non-semi-simple, equal if semi-simple). |
| "Is $x^*$ stable?" (linear) | Eigenvalues: Re $\leq 0$ + semi-simple on imaginary axis. |
| "Is $x^*$ stable?" (nonlinear, hyperbolic) | Linearise. All eigenvalues Re $< 0$ → exp. stable. Any eigenvalue Re $> 0$ → unstable. |
| "Is $x^*$ stable?" (nonlinear, non-hyperbolic) | Build a Lyapunov function. |
| "Is region $M$ positively invariant?" | Check vector field on every boundary piece points inward. |
| "Does a periodic orbit exist?" | Find positively invariant region with repulsive/no equilibria, apply Poincaré–Bendixson. |
| "Is there no periodic orbit?" | Strict Lyapunov function with $\dot V < 0$ off equilibrium. |
| "Find domain of attraction" | Lyapunov function + compact sublevel set inside the basin. |

---

## Appendix B: Common Formulas and Identities

**Fundamental solution identities:**
- $\varphi(0, x) = x$
- $\varphi(t + s, x) = \varphi(t, \varphi(s, x))$ (group property)
- For linear: $\varphi(t, x) = e^{At} x$
- $e^{A(t+s)} = e^{At} e^{As}$, $e^{-At} = (e^{At})^{-1}$

**Matrix exponential for Jordan blocks (2D):**
- $\exp\left(\begin{pmatrix} a & 0 \\ 0 & b\end{pmatrix} t\right) = \begin{pmatrix} e^{at} & 0 \\ 0 & e^{bt}\end{pmatrix}$
- $\exp\left(\begin{pmatrix} a & 1 \\ 0 & a\end{pmatrix} t\right) = e^{at}\begin{pmatrix} 1 & t \\ 0 & 1\end{pmatrix}$
- $\exp\left(\begin{pmatrix} a & b \\ -b & a\end{pmatrix} t\right) = e^{at}\begin{pmatrix} \cos(bt) & \sin(bt) \\ -\sin(bt) & \cos(bt)\end{pmatrix}$

**Gronwall:** $u(t) \leq c + d\int_a^t u(s)\,ds \Rightarrow u(t) \leq c e^{d(t-a)}$.

**Variation of constants:** $\dot{x} = Ax + g(t) \Rightarrow \lambda(t) = e^{A(t-t_0)}x_0 + \int_{t_0}^t e^{A(t-s)}g(s)\,ds$.

**Mean Value Inequality:** $\|f(x) - f(y)\| \leq \|f'(\xi)\| \cdot \|x - y\|$ for some $\xi \in [x,y]$.

---

*Good luck!*