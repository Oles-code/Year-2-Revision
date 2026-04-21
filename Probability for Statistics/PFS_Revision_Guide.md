# MATH50010 — Probability for Statistics: Exam-Focused Revision Notes

> A conceptual summary of the lecture notes, weighted toward topics that recur in the past exams (2021–2025). Proofs are included for results that are (a) short enough to be asked on an exam and (b) likely to actually be asked. Long non-examinable proofs (e.g. strong law of large numbers, almost sure convergence theory) are either sketched or omitted. Key problem-sheet questions to practise are listed at the end.

---

## How the exam is structured (based on 5 recent papers)

Every exam in the last five years has the **same four-question skeleton**:

| Q | Topic | Typical structure |
|---|---|---|
| **Q1** | Sigma algebras & probability measures | Short definitional parts + "prove/disprove this set is a sigma algebra" |
| **Q2** | Random variables, transformations, convergence, MGFs | Varies year to year — one of: densities, CLT, convergence in probability/distribution |
| **Q3** | Bivariate / multivariate Normal | Compute covariance matrix, find conditional distribution, use MGFs of linear combinations |
| **Q4** | Markov chains | Draw diagram, find stationary distribution / hitting probabilities, classify states |

So roughly **25% of marks come from each of these four areas**. Prioritise accordingly.

---

# 1 — Probability Foundations

## 1.1 Sigma algebras

A **sigma algebra** $\mathcal F$ on a set $\Omega$ is a collection of subsets of $\Omega$ such that:

1. $\varnothing \in \mathcal F$,
2. $A \in \mathcal F \implies A^c \in \mathcal F$ (closure under complements),
3. $A_1, A_2, \dots \in \mathcal F \implies \bigcup_{i=1}^\infty A_i \in \mathcal F$ (closure under **countable** unions).

If (3) only holds for finite unions, it's an **algebra**, not a sigma algebra.

**The conceptual picture.** Think of $\Omega$ as the set of possible outcomes and $\mathcal F$ as the set of "questions you can answer" about the outcome. "Was an even number rolled?" is a question, and the corresponding set $\{2,4,6\}$ is an event. The axioms say: you can always ask "did something happen" (closure gives $\Omega$), you can negate any question (complements), and you can take countable "or"s of questions (countable unions). Taking **countable** rather than just finite unions is the one subtle bit — it's what lets us define things like $\{X \leq x\}$ as a limit of $\{X \leq x - 1/n\}$.

### Derived closure properties (very common exam fodder)

From the three axioms you can derive:

- **Closure under countable intersections.** By De Morgan: $\bigcap_i A_i = \big(\bigcup_i A_i^c\big)^c$. Both operations preserve $\mathcal F$.
- **Closure under set differences.** $A \setminus B = A \cap B^c$.
- **Closure under symmetric differences.** $A \triangle B = (A \setminus B) \cup (B \setminus A)$.
- **Closure under $A^c \cap B$** (a 2022 exam warm-up).

The pattern: write the operation you're given in terms of complements, unions, intersections — each of which keeps you inside $\mathcal F$. This is the template for most "is this set operation closed?" questions.

### The intersection of sigma algebras is a sigma algebra

**Proposition.** If $\{\mathcal F_i : i \in I\}$ is a non-empty collection of sigma algebras on $\Omega$, then $\bigcap_{i \in I} \mathcal F_i$ is a sigma algebra.

**Proof.** Three things to check.

1. $\varnothing \in \mathcal F_i$ for every $i$, so $\varnothing$ is in the intersection.
2. If $E$ is in the intersection, then $E \in \mathcal F_i$ for every $i$, hence $E^c \in \mathcal F_i$ for every $i$, hence $E^c$ is in the intersection.
3. If $E_1, E_2, \dots$ are all in the intersection, then they all lie in each $\mathcal F_i$, so $\bigcup_j E_j \in \mathcal F_i$ for every $i$, so the union is in the intersection. $\square$

This proof **really does come up in various forms** — the 2025 paper asks about an intersection of sigma algebras directly. The 30-second version of the proof is: "Each property is quantified over all $i$, so if it holds in every $\mathcal F_i$ it holds in the intersection." Write it out precisely for full marks.

> **Watch out.** The *union* of two sigma algebras is **not** generally a sigma algebra. Minimal counterexample: $\Omega = \{1,2,3\}$, $\mathcal F_1 = \{\varnothing, \{1\}, \{2,3\}, \Omega\}$, $\mathcal F_2 = \{\varnothing, \{2\}, \{1,3\}, \Omega\}$. Their union contains $\{1\}$ and $\{2\}$ but not $\{1,2\}$.

### The Borel sigma algebra

**$\mathcal B$ = the smallest sigma algebra on $\mathbb R$ containing every open interval.**

"Smallest" makes sense because (i) the power set of $\mathbb R$ is *one* sigma algebra containing all open intervals, so non-empty collection, (ii) the intersection of all such sigma algebras is itself a sigma algebra (previous prop), (iii) this intersection is by construction contained in every one of them.

Things you can build from open intervals using countable operations — hence Borel:
- closed intervals $[a,b] = \big((-\infty, a) \cup (b, \infty)\big)^c$;
- singletons $\{a\} = \bigcap_n (a - 1/n, a + 1/n)$;
- the integers $\mathbb Z = \bigcup_{n \in \mathbb Z} \{n\}$;
- any open or closed set.

Practically: *essentially everything you will ever write down on paper* is a Borel set. You only meet non-Borel sets through the axiom of choice.

## 1.2 Probability measures (Kolmogorov axioms)

A **probability measure** on $(\Omega, \mathcal F)$ is $\Pr: \mathcal F \to [0,1]$ satisfying:

1. $\Pr(A) \geq 0$ for all $A \in \mathcal F$,
2. $\Pr(\Omega) = 1$,
3. **Countable additivity:** for pairwise disjoint $A_1, A_2, \ldots \in \mathcal F$, $\Pr\!\big(\bigcup_i A_i\big) = \sum_i \Pr(A_i)$.

The triple $(\Omega, \mathcal F, \Pr)$ is a **probability space**.

Everything else follows: $\Pr(\varnothing) = 0$, $\Pr(A^c) = 1 - \Pr(A)$, monotonicity ($A \subseteq B \implies \Pr(A) \leq \Pr(B)$), inclusion–exclusion, and the law of total probability.

### The Continuity Property (examinable proof)

**Proposition.** If $(A_n)_{n \geq 1}$ is **increasing** ($A_1 \subseteq A_2 \subseteq \cdots$) and $A = \bigcup_n A_n$, then $\Pr(A) = \lim_n \Pr(A_n)$. Similarly if $(B_n)$ is decreasing with $B = \bigcap_n B_n$, then $\Pr(B) = \lim_n \Pr(B_n)$.

**Proof (increasing case).** The trick is to turn an increasing union into a *disjoint* union so we can use countable additivity:
$$A = A_1 \cup (A_2 \setminus A_1) \cup (A_3 \setminus A_2) \cup \cdots$$
These are pairwise disjoint. Applying countable additivity,
$$\Pr(A) = \Pr(A_1) + \sum_{i=1}^\infty \big[\Pr(A_{i+1}) - \Pr(A_i)\big] = \lim_{n} \Pr(A_n),$$
by telescoping. The decreasing case follows by taking complements. $\square$

This matters because it's the key technical tool that turns "limit of events" into "limit of probabilities" — it's what makes the whole theory go.

## 1.3 Independence — the distinction that catches people out

Events $A_1, \ldots, A_n$ are **mutually independent** if for *every* subcollection $A_{i_1}, \ldots, A_{i_k}$,
$$\Pr\!\Big(\bigcap_{j=1}^k A_{i_j}\Big) = \prod_{j=1}^k \Pr(A_{i_j}).$$

Two strictly weaker conditions exist, and the lecturer *wants* you to see they're weaker:

**(a) Pairwise independence is not enough.** Flip a fair coin twice. Let $A = \{HH, HT\}$ (first is H), $B = \{HH, TH\}$ (second is H), $C = \{HT, TH\}$ (exactly one H). All three pairs are independent (each pairwise intersection has probability $1/4 = \tfrac12 \cdot \tfrac12$), but $\Pr(A \cap B \cap C) = \Pr(\varnothing) = 0 \neq 1/8$.

**(b) Factorisation over the entire collection is not enough.** Flip three coins. Choose $A, B, C$ each with probability $1/2$ such that $\Pr(A \cap B \cap C) = 1/8$ but $\Pr(A \cap B) \neq 1/4$. Again: fails independence despite the "full" factorisation holding.

The moral: check *every* subcollection.

---

# 2 — Random Variables

## 2.1 Definition and measurability

A **random variable** on $(\Omega, \mathcal F, \Pr)$ is a function $X : \Omega \to \mathbb R$ such that $X^{-1}(B) \in \mathcal F$ for every Borel set $B \in \mathcal B$.

The motivation: we want $\Pr(X \in B)$ to be defined for Borel sets $B$, and that means the preimage $\{\omega : X(\omega) \in B\}$ needs to be an event in our sigma algebra.

### The "only check intervals" shortcut

**Proposition.** $X$ is a random variable iff $\{X \leq x\} \in \mathcal F$ for every $x \in \mathbb R$.

**Proof sketch.** $\Rightarrow$ is immediate since $(-\infty, x] \in \mathcal B$. For $\Leftarrow$: the collection $\mathcal A = \{B \in \mathcal B : X^{-1}(B) \in \mathcal F\}$ is a sigma algebra (easy to check). It contains every $(-\infty, x]$; from those, via countable operations, you get every open interval $(a, b) = (-\infty, b) \cap (a, \infty)$, where $(-\infty, b) = \bigcup_n (-\infty, b - 1/n]$ and $(a, \infty) = (-\infty, a]^c$. So $\mathcal A$ contains all open intervals, and since $\mathcal B$ is the *smallest* such sigma algebra, $\mathcal A = \mathcal B$. $\square$

This is the standard way to check measurability in practice: you almost never work with arbitrary Borel sets.

## 2.2 CDFs — the three defining properties

The **CDF** is $F_X(x) = \Pr(X \leq x)$. It satisfies:

1. **Non-decreasing.**
2. $\lim_{x \to -\infty} F_X(x) = 0$ and $\lim_{x \to \infty} F_X(x) = 1$.
3. **Right-continuous:** $\lim_{x \downarrow x_0} F_X(x) = F_X(x_0)$.

All three follow from the continuity property of $\Pr$ applied to increasing/decreasing sequences of events of the form $\{X \leq x_n\}$.

### Proof of property 3 (right-continuity) — short and examinable

Take any sequence $x_n \downarrow x_0$. Set $B_n = \{\omega : X(\omega) \leq x_n\}$. Then $(B_n)$ is a *decreasing* sequence of events with $\bigcap_n B_n = \{X \leq x_0\}$ (since $X(\omega) \leq x_n$ for all $n$ iff $X(\omega) \leq x_0$). By the continuity property, $\Pr(B_n) \to \Pr(\bigcap_n B_n)$, i.e. $F_X(x_n) \to F_X(x_0)$. $\square$

Analogous idea for the other two: left-hand limits use increasing sequences $\{X \leq x_n\}$ with $x_n \to -\infty$.

> **Subtle point** (came up in 2025 Q1): $F_X$ is always right-continuous but **not** necessarily left-continuous. The left-limit $F_X(x^-) = \lim_{y \uparrow x} F_X(y) = \Pr(X < x)$, which differs from $F_X(x) = \Pr(X \leq x)$ precisely by $\Pr(X = x)$.

## 2.3 Three types of random variable

- **Discrete:** $F_X$ is a weighted sum of point masses; equivalently there's a PMF $f_X(x) = \Pr(X = x)$ supported on a countable set.
- **Continuous:** $F_X$ is a continuous function. Implies $\Pr(X = x) = 0$ for every $x$.
- **Absolutely continuous:** there exists a PDF $f_X \geq 0$ with $F_X(x) = \int_{-\infty}^x f_X(t)\,dt$.

Absolutely continuous $\Rightarrow$ continuous, but not conversely (the Cantor distribution is the textbook counterexample — continuous CDF, no PDF, non-examinable but worth knowing exists so you distinguish the two in exam definitions).

**Lebesgue decomposition** (worth knowing the statement): any CDF splits uniquely as
$$F = \alpha F_c + \beta F_d + \gamma F_s$$
with $\alpha + \beta + \gamma = 1$ and $F_c, F_d, F_s$ the absolutely continuous, discrete, and singular parts.

**Why $\Pr(X = x) = 0$ for continuous $X$.** Take $x_n \uparrow x$. Then
$$\Pr(X = x) = \Pr(X \leq x) - \Pr(X < x) = F_X(x) - \lim_n F_X(x_n) = 0$$
by continuity of $F_X$.

## 2.4 Transformations — two methods

### Method 1: CDF method (always works)

For $Y = g(X)$, compute $F_Y(y) = \Pr(g(X) \leq y)$ by describing $\{g(X) \leq y\}$ as an event about $X$, then differentiate to get $f_Y$. This is bulletproof for non-monotonic $g$.

Classic example: $Y = \sin X$ with $X \sim \text{Unif}[0, 2\pi]$. Here $\{\sin X \leq y\}$ breaks into two sub-intervals of $X$, and you add their probabilities.

### Method 2: Jacobian formula (strictly monotonic $g$)

**Proposition.** If $X$ is absolutely continuous with density $f_X$ and $g: \mathbb R \to \mathbb R$ is strictly monotonic and differentiable, then $Y = g(X)$ has density
$$f_Y(y) = f_X(g^{-1}(y))\,\Big|\tfrac{d}{dy} g^{-1}(y)\Big|.$$

**Proof (increasing case).** $\{g(X) \leq y\} = \{X \leq g^{-1}(y)\}$, so $F_Y(y) = F_X(g^{-1}(y))$. Differentiate via chain rule. Decreasing case: flip the inequality, get $F_Y(y) = 1 - F_X(g^{-1}(y))$, and the absolute value handles the sign. $\square$

**The intuitive picture.** $f_X(x)\,dx$ is the probability mass in a tiny interval of width $dx$. Under $y = g(x)$, that interval gets mapped to one of width $|g'(x)|\,dx$. Probability mass must be preserved, so the density at $y$ is $f_X(x) / |g'(x)| = f_X(g^{-1}(y)) \cdot |dg^{-1}/dy|$.

### Probability Integral Transform

**Proposition.** If $U \sim \text{Unif}[0,1]$ and $F$ is a strictly increasing CDF, then $X := F^{-1}(U)$ has CDF $F$.

**Proof.** $\Pr(X \leq x) = \Pr(F^{-1}(U) \leq x) = \Pr(U \leq F(x)) = F(x)$, since for $0 \leq u \leq 1$, $F_U(u) = u$. $\square$

**Why you care.** This is the foundation of random-number generation: to simulate from any distribution with an invertible CDF, sample $U$ uniformly and apply $F^{-1}$. E.g. for $\text{Exp}(\beta)$ with $F(x) = 1 - e^{-x/\beta}$, you get $X = -\beta \log(1 - U)$. Comes up implicitly in 2024 Q2 (Exponential $\to$ Gumbel).

## 2.5 Location/scale families

If $Z$ has PDF $f_Z$, then:

- **Location:** $X = \mu + Z$ has PDF $f_Z(x - \mu)$.
- **Scale:** $X = \sigma Z$ has PDF $\tfrac{1}{\sigma} f_Z(x/\sigma)$, $\sigma > 0$.
- **Location-scale:** $X = \mu + \sigma Z$ has PDF $\tfrac{1}{\sigma} f_Z\big(\tfrac{x - \mu}{\sigma}\big)$.

For the Normal, $\mu$ is literally the mean and $\sigma$ the standard deviation; this is **not** true in general (e.g. for Gamma, $\sigma$ is a scale but not the standard deviation).

---

# 3 — Multivariate Random Variables

## 3.1 Joint distributions

Joint CDF: $F_{XY}(x,y) = \Pr(X \leq x, Y \leq y)$. If a joint density exists, the marginals are obtained by integrating one variable out:
$$f_X(x) = \int_{-\infty}^\infty f_{XY}(x, y)\,dy.$$

### Independence in terms of densities

$X, Y$ are independent iff the joint density **factorises for all $x, y$** as $f_{XY}(x, y) = g(x) h(y)$ where $g, h$ are *properly defined on all of $\mathbb R$*.

> **Trap.** A density like $f(x,y) = 18(1-u)\frac{v^2}{u^2}$ on $\{0 < v < u < 1\}$ *looks* like it factorises, but the **support** is not a product set — the constraint $v < u$ couples the two variables. Same goes for anything on a triangular or circular region. So $X, Y$ are not independent despite the algebraic factorisation.

## 3.2 Covariance, correlation, Cauchy–Schwarz

$$\text{Cov}(X, Y) = \mathbb E[(X - \mu_X)(Y - \mu_Y)] = \mathbb E(XY) - \mathbb E(X)\mathbb E(Y).$$
$$\rho(X, Y) = \frac{\text{Cov}(X, Y)}{\sqrt{\text{Var}(X)\,\text{Var}(Y)}} \in [-1, 1].$$

**Key facts:**
- Independent $\Rightarrow$ covariance zero. Converse false in general.
- The bound $|\rho| \leq 1$ follows from **Cauchy–Schwarz**: $|\mathbb E(XY)| \leq \sqrt{\mathbb E(X^2)\mathbb E(Y^2)}$, with equality iff $X$ and $Y$ are linearly related with probability 1.

**Counterexample to the converse** (a classic problem-sheet question): Let $X \sim \text{Unif}(-1,1)$ and $Y = X^2$. Then $Y$ is completely determined by $X$, so certainly not independent. But $\text{Cov}(X, Y) = \mathbb E(X^3) - \mathbb E(X)\mathbb E(X^2) = 0 - 0 = 0$ by symmetry.

## 3.3 Change of variables (Jacobian, bivariate form)

**Proposition.** If $T : D \subseteq \mathbb R^2 \to R \subseteq \mathbb R^2$ is an invertible smooth map, and $(U, V) = T(X, Y)$, then
$$f_{UV}(u, v) = f_{XY}(x(u, v), y(u, v)) \cdot |J(u, v)|,$$
where $J(u, v) = \det\!\begin{pmatrix} \partial x / \partial u & \partial x / \partial v \\ \partial y / \partial u & \partial y / \partial v \end{pmatrix}$ is the Jacobian of the *inverse* map.

The workflow:

1. Write down the mapping $(U, V) = T(X, Y)$.
2. Invert it to get $(X, Y)$ in terms of $(U, V)$.
3. Compute the Jacobian.
4. Substitute $f_{XY}$ and multiply by $|J|$.
5. **Determine the support of $(U, V)$** — this is where most mistakes happen. Map the boundary of the domain, not just the interior.

Step 5 is the step students lose marks on.

## 3.4 Conditional distributions

For continuous $X, Y$ with joint density $f_{XY}$:
$$f_{Y|X}(y \mid x) = \frac{f_{XY}(x, y)}{f_X(x)}, \quad \text{provided } f_X(x) > 0.$$

**Why you can't just condition on $\{X = x\}$ directly.** For continuous $X$, $\Pr(X = x) = 0$, so the naive definition $\Pr(\cdot \mid X = x) = \Pr(\cdot \cap \{X = x\})/\Pr(X = x)$ is $0/0$. The fix is to condition on $\{X \in (x, x + h)\}$ (positive probability!) and take $h \to 0$ using L'Hôpital's rule; the derivative of an integral over $(x, x + h)$ w.r.t. $h$ is the integrand evaluated at $x + h$.

## 3.5 Bivariate Normal — THE question 3 staple

The **standard** bivariate Normal with correlation $\rho \in (-1, 1)$ has density
$$f(x, y) = \frac{1}{2\pi\sqrt{1 - \rho^2}} \exp\!\Big(-\frac{1}{2(1 - \rho^2)}(x^2 - 2\rho xy + y^2)\Big).$$

### Key technique: completing the square

$$x^2 - 2\rho xy + y^2 = (x - \rho y)^2 + (1 - \rho^2) y^2.$$

This is the single most useful algebraic identity in Q3. Substituting:
$$f(x, y) = \frac{1}{\sqrt{2\pi}} e^{-y^2/2} \cdot \frac{1}{\sqrt{2\pi(1-\rho^2)}} \exp\!\Big(-\frac{(x - \rho y)^2}{2(1 - \rho^2)}\Big).$$

Reading this off:

- **Marginal of $Y$:** $Y \sim N(0, 1)$ (the first factor).
- **Conditional of $X$ given $Y = y$:** $X \mid Y = y \sim N(\rho y, 1 - \rho^2)$ (the second factor).

Both come for free once you've completed the square. The conditional mean is a linear function of $y$ with slope $\rho$, and the conditional variance $1 - \rho^2$ is *independent* of $y$ (a very special property of Normal distributions).

### Computing Cov$(X, Y) = \rho$ using the law of iterated expectations

$$\mathbb E(XY) = \mathbb E[Y \cdot \mathbb E(X \mid Y)] = \mathbb E[Y \cdot \rho Y] = \rho \mathbb E(Y^2) = \rho.$$

This is vastly cleaner than doing a double integral directly. **Always reach for iterated expectations** when you have a clean conditional distribution.

### General (non-standard) bivariate Normal

If you want arbitrary means and variances, write the density in vector form:
$$f_{\mathbf X}(\mathbf x \mid \boldsymbol\mu, \Sigma) = \frac{1}{2\pi\sqrt{\det \Sigma}} \exp\!\Big(-\tfrac12 (\mathbf x - \boldsymbol\mu)^T \Sigma^{-1} (\mathbf x - \boldsymbol\mu)\Big),$$
with $\Sigma = \begin{pmatrix} \sigma_X^2 & \rho \sigma_X \sigma_Y \\ \rho \sigma_X \sigma_Y & \sigma_Y^2 \end{pmatrix}$.

**Mahalanobis distance.** The quadratic form $(\mathbf x - \boldsymbol\mu)^T \Sigma^{-1}(\mathbf x - \boldsymbol\mu)$ is a generalised squared distance from the mean. Contours of equal density are **ellipses** whose principal axes are the eigenvectors of $\Sigma$.

```
        Y
        │     ╱── (eigenvector 1, eigenvalue λ₁: long axis)
        │    ╱
        │   ╱
        │  ╱   ╲
        │ ╱     ╲── (eigenvector 2, eigenvalue λ₂: short axis)
        │╱_______╲
        ●        ╲
       (μ_X,μ_Y)  
   ─────┼────────── X
        │
```

## 3.6 Multivariate Normal — general $d$ dimensions

$$f_{\mathbf X}(\mathbf x \mid \boldsymbol\mu, \Sigma) = \frac{1}{(2\pi)^{d/2} (\det \Sigma)^{1/2}} \exp\!\Big(-\tfrac12 (\mathbf x - \boldsymbol\mu)^T \Sigma^{-1} (\mathbf x - \boldsymbol\mu)\Big).$$

### Things you need to have at your fingertips

**1. Variance-covariance matrix** $\Sigma_{ij} = \text{Cov}(X_i, X_j)$. It's symmetric and **positive semi-definite**, because for any constant vector $\mathbf a$,
$$\mathbf a^T \Sigma \mathbf a = \text{Var}(\mathbf a^T \mathbf X) \geq 0.$$
Positive definite (invertible $\Sigma$) iff no non-trivial linear combination of the $X_i$ is almost-surely constant.

**2. Linear combinations of MVN are Normal.** For any $\mathbf a \in \mathbb R^d$:
$$\mathbf a^T \mathbf X \sim N(\mathbf a^T \boldsymbol\mu, \mathbf a^T \Sigma \mathbf a).$$

**3. Affine transformations.** If $\mathbf X \sim \text{MVN}_d(\boldsymbol\mu, \Sigma)$ and $A$ is invertible, then
$$\mathbf Y = A\mathbf X \sim \text{MVN}_d(A\boldsymbol\mu, A\Sigma A^T).$$

The proof is a Jacobian computation (see the lecture notes; the bookkeeping is somewhat long but mechanical). The **result** — $A\Sigma A^T$ — is what you need to remember.

**4. Decorrelation via eigendecomposition.** Because $\Sigma$ is symmetric and positive definite, there's an orthogonal matrix $Q$ (columns are eigenvectors of $\Sigma$) with
$$Q^T \Sigma Q = \text{diag}(\lambda_1, \ldots, \lambda_d).$$
Then $\mathbf Z = Q^T(\mathbf X - \boldsymbol\mu)$ has entries that are **independent** $N(0, \lambda_i)$ random variables. Geometrically, you're rotating to the principal-axis coordinates of the ellipsoid.

**5. Setting $M = \Sigma^{-1/2}$** (came up in 2025 Q3): if $\mathbf X \sim \text{MVN}(\boldsymbol\mu, \Sigma)$ and $M = \Sigma^{-1/2}$, then $M\mathbf X - M\boldsymbol\mu \sim \text{MVN}(\mathbf 0, I)$ — standardised to uncorrelated components with unit variance. The "$\Sigma^{-1/2}$" here means any matrix square root such that $M M^T = \Sigma^{-1}$; the symmetric one $\Sigma^{-1/2}$ is the standard choice.

## 3.7 Order statistics

Given an iid sample $X_1, \ldots, X_n \sim f_X$ with CDF $F_X$, arrange them in increasing order as $Y_1 < Y_2 < \cdots < Y_n$. The joint density is
$$f(y_1, \ldots, y_n) = n! \prod_{i=1}^n f_X(y_i), \quad y_1 < y_2 < \cdots < y_n,$$
and zero elsewhere (the $n!$ accounts for all the orderings of the original sample that map to the same ordered tuple).

### The two you must know by heart

**Minimum and maximum.**
- $\Pr(Y_1 > y) = \Pr(\text{all } X_i > y) = (1 - F_X(y))^n$, so $F_{Y_1}(y) = 1 - (1 - F_X(y))^n$.
- $\Pr(Y_n \leq y) = \Pr(\text{all } X_i \leq y) = F_X(y)^n$.

**Consequence for exponentials.** If $X_i \sim \text{Exp}(\lambda)$, then $Y_1 = \min X_i \sim \text{Exp}(n\lambda)$. This shows up constantly — waiting-time problems, competing risks, the first event of many Poisson processes.

### General $k$-th order statistic

$$f_{Y_k}(y) = k \binom{n}{k} f_X(y)\, F_X(y)^{k-1}\,(1 - F_X(y))^{n-k}.$$

**Intuition:** pick one of the $n$ samples to sit at $y$ (contributing $f_X(y)\,dy$), and out of the remaining $n - 1$, exactly $k - 1$ must be below $y$ and $n - k$ above — a Binomial count. Combining:
$$k \binom{n}{k} = n \binom{n-1}{k-1}, \quad \text{(one way to see the combinatorics)}.$$

You can derive this by differentiating the CDF $F_{Y_k}(y) = \Pr(\text{at least }k\text{ of the }X_i \leq y)$, which is a Binomial tail sum $\sum_{j=k}^n \binom{n}{j} F_X(y)^j (1 - F_X(y))^{n-j}$.

---

# 4 — Convergence of Random Variables

Four modes of convergence appear:

```
  a.s.  ──►  in probability  ──►  in distribution
                   ▲                     │
                   └── (if limit is a ───┘
                         constant)
```

In general: almost sure $\Rightarrow$ in probability $\Rightarrow$ in distribution, and **none of the reverse implications hold**. The *one* exception: convergence in distribution to a constant implies convergence in probability. This is essential and comes up.

## 4.1 Definitions

- **In probability:** $X_n \xrightarrow{P} X$ if $\Pr(|X_n - X| \geq \varepsilon) \to 0$ for all $\varepsilon > 0$.
- **In distribution:** $X_n \xrightarrow{D} X$ if $F_n(x) \to F_X(x)$ at every continuity point of $F_X$.
- **Almost surely:** $X_n \xrightarrow{a.s.} X$ if $\Pr(\{\omega : X_n(\omega) \to X(\omega)\}) = 1$ (non-examinable in detail but expect to know the definition).

### Why "continuity points of $F_X$" in the distribution definition?

Consider $X_n \sim \text{Unif}(-1/n, 1/n)$. Intuitively this should converge to the constant $0$. But at $x = 0$, $F_n(0) = 1/2$ for every $n$, while $F_X(0) = 1$ (since $X = 0$ almost surely, $F$ jumps from $0$ to $1$ there). So $F_n(0) \not\to F_X(0)$. The fix: only require convergence at continuity points. $x = 0$ is a discontinuity of $F_X$, so it's excluded, and everything works.

## 4.2 Markov's inequality and Chebyshev — both examinable

**Markov:** For $X \geq 0$ and $a > 0$:
$$\Pr(X \geq a) \leq \frac{\mathbb E(X)}{a}.$$

**Proof.** Let $\mathbb 1_{\{X \geq a\}}$ be the indicator. Since $X \geq a \cdot \mathbb 1_{\{X \geq a\}}$ (check both cases pointwise), take expectations: $\mathbb E(X) \geq a \cdot \Pr(X \geq a)$. $\square$

**Chebyshev:** For $X$ with mean $\mu$ and variance $\sigma^2$:
$$\Pr(|X - \mu| \geq \varepsilon) \leq \frac{\sigma^2}{\varepsilon^2}.$$

**Proof.** Apply Markov to the non-negative variable $Y = (X - \mu)^2$ with threshold $\varepsilon^2$:
$$\Pr(|X - \mu| \geq \varepsilon) = \Pr((X - \mu)^2 \geq \varepsilon^2) \leq \frac{\mathbb E[(X - \mu)^2]}{\varepsilon^2} = \frac{\sigma^2}{\varepsilon^2}. \quad \square$$

These are your go-to tools for "prove convergence in probability". The pattern is: bound $\Pr(|X_n - \mu| \geq \varepsilon)$ by something that goes to zero.

## 4.3 Weak Law of Large Numbers — examinable proof

**Theorem.** If $X_1, X_2, \ldots$ are iid with mean $\mu$ and variance $\sigma^2 < \infty$, then $\bar X_n \xrightarrow{P} \mu$.

**Proof.** By linearity $\mathbb E(\bar X_n) = \mu$, and by independence $\text{Var}(\bar X_n) = \sigma^2/n$. Chebyshev gives
$$\Pr(|\bar X_n - \mu| \geq \varepsilon) \leq \frac{\sigma^2}{n\varepsilon^2} \to 0. \quad \square$$

That's the entire proof. Memorise it — it's about 4 lines and it's almost free marks.

> Note: independence isn't really needed; *uncorrelated* with finite variance suffices, because that's all we use when computing $\text{Var}(\bar X_n)$.

## 4.4 Convergence in probability $\Rightarrow$ in distribution (examinable)

**Proposition.** If $X_n \xrightarrow{P} X$, then $X_n \xrightarrow{D} X$.

**Proof idea.** Let $x$ be a continuity point of $F_X$. The crucial set inclusion:
$$\{X_n \leq x\} \subseteq \{X \leq x + \varepsilon\} \cup \{|X_n - X| > \varepsilon\}$$
because if $X_n \leq x$ and $|X_n - X| \leq \varepsilon$, then $X \leq X_n + \varepsilon \leq x + \varepsilon$. By a union bound:
$$F_n(x) \leq F_X(x + \varepsilon) + \Pr(|X_n - X| > \varepsilon).$$
Similarly $F_X(x - \varepsilon) - \Pr(|X_n - X| > \varepsilon) \leq F_n(x)$. Taking $n \to \infty$, the probability terms vanish, and we're sandwiched between $F_X(x - \varepsilon)$ and $F_X(x + \varepsilon)$. Since $x$ is a continuity point, send $\varepsilon \to 0$. $\square$

The reverse implication fails in general (take $U \sim \text{Unif}(-1, 1)$ and $U_n = -U$: identically distributed, so $U_n \xrightarrow{D} U$ trivially, but $|U_n - U| = 2|U|$ doesn't go to zero).

**Special case that does reverse:** if $X_n \xrightarrow{D} c$ for a constant $c$, then $X_n \xrightarrow{P} c$. The proof uses that the CDF of a constant has a single discontinuity, so you can squeeze $\Pr(|X_n - c| \geq \varepsilon) \leq F_n(c - \varepsilon) + 1 - F_n(c + \varepsilon/2) \to 0 + 0$.

## 4.5 Limit events, {$A_n$ i.o.}, Borel–Cantelli

For a sequence of events $(A_n)$:
$$\{A_n \text{ i.o.}\} := \limsup_n A_n = \bigcap_{N=1}^\infty \bigcup_{n=N}^\infty A_n.$$
$$\{A_n \text{ a.a.}\} := \liminf_n A_n = \bigcup_{N=1}^\infty \bigcap_{n=N}^\infty A_n.$$

"Infinitely often" = "for every $N$, $A_n$ happens for some $n \geq N$". "Almost always" = "there exists $N$ such that $A_n$ happens for every $n \geq N$".

### Borel–Cantelli Lemmas (very common)

**BC1.** If $\sum_n \Pr(A_n) < \infty$, then $\Pr(\{A_n \text{ i.o.}\}) = 0$.

**BC2.** If $\sum_n \Pr(A_n) = \infty$ **and** the $A_n$ are independent, then $\Pr(\{A_n \text{ i.o.}\}) = 1$.

**Proof of BC1.** Set $B_N = \bigcup_{n \geq N} A_n$. Then $(B_N)$ is decreasing and $\{A_n \text{ i.o.}\} = \bigcap_N B_N$. By continuity,
$$\Pr(\{A_n \text{ i.o.}\}) = \lim_N \Pr(B_N) \leq \lim_N \sum_{n = N}^\infty \Pr(A_n) = 0,$$
since the tail of a convergent series tends to zero. $\square$

**Proof sketch for BC2.** Show $\Pr(\{A_n \text{ i.o.}\}^c) = \Pr(\{A_n^c \text{ a.a.}\}) = 0$. The key inequality:
$$\Pr\!\Big(\bigcap_{n=N}^M A_n^c\Big) = \prod_{n=N}^M (1 - \Pr(A_n)) \leq \exp\!\Big(-\sum_{n=N}^M \Pr(A_n)\Big) \to 0$$
as $M \to \infty$, using $1 - x \leq e^{-x}$ and divergence of the sum. $\square$

**Asymmetry alert.** BC1 requires **nothing** beyond $\sum \Pr(A_n) < \infty$. BC2 crucially needs **independence** — without it, you can make the sum diverge and still have $\{A_n \text{ i.o.}\}$ of probability $0$ (e.g. $A_n = A$ for a fixed event with $\Pr(A) = 0.5$; the "infinitely often" event is just $A$ with probability $0.5$).

---

# 5 — Central Limit Theorem

## 5.1 Moment Generating Functions

$M_X(t) = \mathbb E(e^{tX})$, defined for $t$ in some neighbourhood of $0$.

### Properties (all short, all examinable)

**1. Affine transformation.** If $Y = aX + b$, then $M_Y(t) = e^{bt} M_X(at)$.
$$M_Y(t) = \mathbb E(e^{t(aX + b)}) = e^{bt} \mathbb E(e^{(at)X}) = e^{bt} M_X(at). \quad \square$$

**2. Sum of independents.** If $X, Y$ independent and $Z = X + Y$, then $M_Z(t) = M_X(t) M_Y(t)$.
$$M_Z(t) = \mathbb E(e^{t(X+Y)}) = \mathbb E(e^{tX} e^{tY}) = \mathbb E(e^{tX}) \mathbb E(e^{tY}) = M_X(t) M_Y(t). \quad \square$$

**3. Derivatives give moments.** $\dfrac{d^k}{dt^k} M_X(t) \Big|_{t=0} = \mathbb E(X^k)$. Expand $e^{tX}$ as a power series and interchange expectation and sum (justified where the MGF exists).

**4. Uniqueness.** If two RVs have the same MGF in a neighbourhood of $0$, they have the same distribution. *(Stated without proof; use freely.)*

**5. Continuity.** If $M_{X_n}(t) \to M_X(t)$ pointwise on some interval containing $0$, then $X_n \xrightarrow{D} X$. *(Stated without proof; this is the engine behind the CLT proof.)*

### MGFs you should be able to reproduce

- $Z \sim N(0, 1)$: $M_Z(t) = e^{t^2/2}$ — by completing the square inside the integral.
- $X \sim N(\mu, \sigma^2)$: $M_X(t) = e^{\mu t + \sigma^2 t^2 / 2}$ — by the affine property.
- $X \sim \text{Gamma}(\alpha, \lambda)$ (rate parameterisation): $M_X(t) = \big(\frac{\lambda}{\lambda - t}\big)^\alpha$ for $t < \lambda$.
- $X \sim \text{Exp}(\lambda)$: $M_X(t) = \frac{\lambda}{\lambda - t}$ (the $\alpha = 1$ case above).
- $X \sim \text{Poisson}(\lambda)$: $M_X(t) = \exp(\lambda(e^t - 1))$.

**Classic trick:** sum of $n$ iid Exp($\lambda$)s has MGF $\big(\frac{\lambda}{\lambda-t}\big)^n$ = MGF of Gamma($n, \lambda$). By uniqueness, the sum *is* Gamma$(n, \lambda)$.

### Cumulant Generating Function

$K_X(t) = \log M_X(t)$. Nice because $K'_X(0) = \mathbb E(X)$ and $K''_X(0) = \text{Var}(X)$ (came up on PS6 and 2022 Q3). Also: $K_{X + Y}(t) = K_X(t) + K_Y(t)$ when $X, Y$ independent — additivity of cumulants.

**Derivation of $K''_X(0) = \text{Var}(X)$.** $K' = M'/M$, so $K'' = M''/M - (M'/M)^2$. At $t = 0$: $K''(0) = \mathbb E(X^2) - \mathbb E(X)^2 = \text{Var}(X)$. $\square$

## 5.2 The Central Limit Theorem — with full proof

**Theorem (CLT).** Let $X_1, X_2, \ldots$ be iid with mean $\mu$, variance $\sigma^2 < \infty$, and MGF existing on a neighbourhood of $0$. Then
$$\frac{\sqrt{n}(\bar X_n - \mu)}{\sigma} \xrightarrow{D} Z \sim N(0, 1).$$

**Proof.** Let $Y_i = (X_i - \mu)/\sigma$, so the $Y_i$ are standardised: $\mathbb E(Y_i) = 0$, $\text{Var}(Y_i) = 1$. Let $M(t)$ be the common MGF of the $Y_i$. Define
$$Z_n = \frac{\sqrt{n}(\bar X_n - \mu)}{\sigma} = \frac{1}{\sqrt{n}}\sum_{i=1}^n Y_i.$$

By independence and the affine-transformation property of MGFs:
$$M_{Z_n}(t) = \prod_{i=1}^n \mathbb E\!\left[\exp\!\left(\frac{t Y_i}{\sqrt{n}}\right)\right] = M\!\left(\tfrac{t}{\sqrt{n}}\right)^n.$$

Taylor-expand $M$ near $0$. Using $M(0) = 1$, $M'(0) = \mathbb E(Y_i) = 0$, $M''(0) = \mathbb E(Y_i^2) = 1$:
$$M(s) = 1 + \tfrac{s^2}{2} M''(r_s), \qquad |r_s| \leq |s|,$$
where we've used the Lagrange form of the remainder after absorbing the linear term (which vanishes).

Substituting $s = t/\sqrt{n}$:
$$M_{Z_n}(t) = \left(1 + \frac{t^2}{2n} M''(r_n)\right)^n, \qquad |r_n| \leq \tfrac{|t|}{\sqrt n}.$$

As $n \to \infty$, $r_n \to 0$, so by continuity of $M''$ (which exists because MGF is finite in a neighbourhood), $M''(r_n) \to M''(0) = 1$. Using the limit $(1 + a_n/n)^n \to e^a$ whenever $a_n \to a$:
$$M_{Z_n}(t) \to e^{t^2/2}.$$

This is the MGF of $N(0, 1)$, so by MGF continuity, $Z_n \xrightarrow{D} N(0, 1)$. $\square$

### Why the CLT is remarkable

The limit distribution is $N(0, 1)$ **regardless of the parent distribution**, as long as variance is finite. You could start from Bernoulli coin flips, Exponential waiting times, Uniform samples — the sample mean's standardised fluctuations always look Normal for large $n$. This is why Normal approximations are used everywhere in statistics.

### When the CLT fails

If $\sigma^2 = \infty$ (or even $\mathbb E|X| = \infty$), the CLT does not apply. The **Cauchy distribution** is the canonical failure case: $f(x) = \frac{1}{\pi(1 + x^2)}$ has no defined mean, and in fact $\bar X_n$ is *itself* Cauchy-distributed (never concentrates). Sample medians still have a CLT, though — that's model-robustness in action.

---

# 6 — Markov Chains (the every-year Q4 topic)

## 6.1 Definition

A sequence $(X_n)_{n \geq 0}$ on a finite/countable state space $E$ is a **Markov chain** if
$$\Pr(X_n = x_n \mid X_{n-1} = x_{n-1}, \ldots, X_0 = x_0) = \Pr(X_n = x_n \mid X_{n-1} = x_{n-1})$$
for all $n$ and all states. The future depends on the past *only* through the present.

**Time-homogeneous** means the transition probabilities don't depend on $n$: $\Pr(X_{n+1} = j \mid X_n = i) = \Pr(X_1 = j \mid X_0 = i) =: p_{ij}$.

The **transition matrix** $P = (p_{ij})$ is a **stochastic matrix**: rows are non-negative and sum to 1.

To specify the process you also need the **initial distribution** $\lambda_j = \Pr(X_0 = j)$.

## 6.2 Chapman–Kolmogorov Equations (examinable proof)

The **$n$-step transition matrix** $P(n)$ has entries $p_{ij}(n) = \Pr(X_n = j \mid X_0 = i)$.

**Theorem.** $P(m + n) = P(m) P(n)$. Consequently $P(n) = P^n$.

**Proof.** For $m \geq 0, n \geq 1$:
$$p_{ij}(m + n) = \Pr(X_{m+n} = j \mid X_0 = i) = \sum_{l \in E} \Pr(X_{m+n} = j, X_m = l \mid X_0 = i).$$
Break the joint probability using conditioning:
$$= \sum_l \Pr(X_{m+n} = j \mid X_m = l, X_0 = i) \Pr(X_m = l \mid X_0 = i).$$
Apply the Markov property (the first factor becomes $\Pr(X_{m+n} = j \mid X_m = l)$) and time homogeneity ($= \Pr(X_n = j \mid X_0 = l) = p_{lj}(n)$):
$$= \sum_l p_{il}(m) p_{lj}(n).$$
That's exactly the $(i,j)$ entry of the matrix product $P(m) P(n)$. Iterating, $P(n) = P^n$. $\square$

The useful practical corollary: **to find the distribution of $X_n$ starting from distribution $\lambda$, compute $\lambda P^n$**.

## 6.3 Classification of states

### Accessibility and communication

- $i \to j$ (**$j$ accessible from $i$**) if $p_{ij}(n) > 0$ for some $n \geq 0$.
- $i \leftrightarrow j$ (**communicate**) if both $i \to j$ and $j \to i$.

Communication is an **equivalence relation** (reflexive trivially, symmetric by definition, transitive by Chapman–Kolmogorov: if $i \to j$ in $m$ steps and $j \to k$ in $n$ steps, then $p_{ik}(m + n) \geq p_{ij}(m) p_{jk}(n) > 0$).

This partitions $E$ into **communicating classes**. A class $C$ is **closed** if you can't escape: $p_{ij} = 0$ for $i \in C$, $j \notin C$. A chain is **irreducible** if all of $E$ is one communicating class.

### Periodicity

The **period** of state $i$: $d(i) = \gcd\{n \geq 1 : p_{ii}(n) > 0\}$. Aperiodic if $d(i) = 1$, periodic otherwise.

**Key fact (proof is short and examinable):** all states in one communicating class share the same period. The proof uses Chapman–Kolmogorov to construct loops $i \to j \to j \to i$ of length $a + m + b$, forcing divisibility.

### Recurrence and transience

Let $T_j = \min\{n \geq 1 : X_n = j\}$ be the **first passage time** to $j$ (with $T_j = \infty$ if we never reach $j$). Define
$$f_{ij} = \Pr(T_j < \infty \mid X_0 = i),$$
the probability of ever reaching $j$ starting from $i$.

**State $i$ is:**
- **Recurrent** if $f_{ii} = 1$ (almost surely returns to $i$, hence visits $i$ infinitely often).
- **Transient** if $f_{ii} < 1$ (positive probability of never returning; number of visits is geometric, hence finite a.s.).

**Equivalent characterisation (very useful in practice):**
$$i \text{ recurrent} \iff \sum_{n=1}^\infty p_{ii}(n) = \infty.$$

**Proof sketch.** The total number of visits $N_i = \sum_n \mathbb 1_{\{X_n = i\}}$ has $\mathbb E(N_i \mid X_0 = i) = \sum_n p_{ii}(n)$. If $i$ is transient, $N_i$ has a geometric distribution with success probability $1 - f_{ii}$, giving finite expectation $f_{ii}/(1 - f_{ii})$. If recurrent, $N_i = \infty$ a.s. $\square$

**Class property:** if $i \leftrightarrow j$, they are either both recurrent or both transient. *(Useful in exams: establish one state's behaviour, the whole class follows.)*

**Closed+recurrent:** every recurrent communicating class is closed. Equivalently: if you leave a class, it was transient.

**State space decomposition.** $E = T \cup C_1 \cup C_2 \cup \cdots$ where $T$ is the set of transient states and each $C_i$ is a closed recurrent class. For finite state spaces, chains can only spend finitely long in $T$, so at least one $C_i$ exists.

### Positive vs null recurrence

For a recurrent state, define the **mean recurrence time** $\mu_i = \mathbb E(T_i \mid X_0 = i)$.

- **Positive recurrent:** $\mu_i < \infty$.
- **Null recurrent:** $\mu_i = \infty$ (recurrent, but returns take infinitely long on average).

Null recurrence is an infinite-state-space phenomenon — in a finite irreducible chain, all states are positive recurrent.

## 6.4 Hitting probabilities — the $h_i^A$ linear system

Let $H^A = \min\{n \geq 0 : X_n \in A\}$ and $h_i^A = \Pr(H^A < \infty \mid X_0 = i)$ = probability of ever hitting $A$ from $i$.

**Proposition.** The vector $(h_i^A)_{i \in E}$ is the **minimal non-negative solution** of
$$h_i^A = \begin{cases} 1 & i \in A \\ \sum_{j \in E} p_{ij} h_j^A & i \notin A. \end{cases}$$

**Why "minimal"?** The system can have multiple non-negative solutions (e.g. the constant vector $\mathbf 1$ always solves it). The actual hitting-probability vector is the smallest one. In practice this means: find the general solution to the linear system, then impose any boundary conditions to pick out the minimal.

**Workflow:**

1. Identify the target set $A$.
2. For each non-target state $i$, apply the law of total probability by conditioning on $X_1$:
$$h_i^A = \sum_j p_{ij} h_j^A.$$
3. Set $h_i^A = 1$ for $i \in A$.
4. Solve. If unique non-negative solution, you're done; otherwise take the minimal.

**Example pattern (gambler's ruin / absorbing boundaries).** Chain with two absorbing states (say $1$ and $N$), and you want $h_i^N$ = probability of absorption at $N$ starting from $i$. You get a second-order linear recurrence $h_i = p\, h_{i+1} + q\, h_{i-1}$ with $h_1 = 0$, $h_N = 1$. Solve using the characteristic equation.

## 6.5 Stationary Distributions

A **stationary distribution** $\pi = (\pi_j)_{j \in E}$ satisfies:
1. $\pi_j \geq 0$ and $\sum_j \pi_j = 1$ (probability distribution),
2. $\pi P = \pi$ (left eigenvector of $P$ with eigenvalue $1$).

If $X_n \sim \pi$ then $X_{n+1} \sim \pi$ — the distribution is preserved by the dynamics.

### Why you care

For an **irreducible, aperiodic, positive recurrent** Markov chain, the distribution of $X_n$ converges to $\pi$ as $n \to \infty$, independently of the starting distribution. Moreover:
$$\pi_i = \frac{1}{\mu_i},$$
the reciprocal of the mean recurrence time — so $\pi_i$ is literally the long-run fraction of time spent in state $i$.

### How to find $\pi$

Solve the system $\pi P = \pi$ together with $\sum_j \pi_j = 1$. For a $K$-state chain this gives $K + 1$ equations in $K$ unknowns — but one of the $\pi P = \pi$ equations is redundant (they sum to $0$ on each side), so you get a unique answer.

**Small trick (symmetry):** if the chain has an obvious symmetry (e.g. a random walk on a ring), guess the uniform distribution and check it satisfies the stationarity equations.

**Theorem (random walk on a graph).** For a simple symmetric random walk on a finite connected graph, if state $i$ has degree $d_i$, then
$$\pi_i = \frac{d_i}{\sum_j d_j}.$$
More popular vertices get proportionally more probability mass — intuitive.

### Two-state chain worked example

$P = \begin{pmatrix} 1 - \alpha & \alpha \\ \beta & 1 - \beta \end{pmatrix}$, $\alpha, \beta \in (0, 1)$.

Solve $(\pi_1, \pi_2) P = (\pi_1, \pi_2)$ with $\pi_1 + \pi_2 = 1$:
$$\pi_1 (1 - \alpha) + \pi_2 \beta = \pi_1 \implies \pi_2 \beta = \pi_1 \alpha.$$

Combined with normalisation: $\pi_1 = \beta/(\alpha + \beta)$, $\pi_2 = \alpha/(\alpha + \beta)$.

The "balance" interpretation: at stationarity, the flow from $1 \to 2$ (which is $\pi_1 \alpha$) equals the flow from $2 \to 1$ (which is $\pi_2 \beta$). This is **detailed balance** and it's a handy sanity check on any stationarity calculation.

## 6.6 Simple random walk — conceptual anchor

The SRW on $\mathbb Z$ with $\Pr(+1) = p$, $\Pr(-1) = 1-p$. Some facts that come up in various forms:

- $X_n = \sum_{k=1}^n Z_k$ where $Z_k \in \{-1, +1\}$ iid.
- $\Pr(X_n = j \mid X_0 = 0) = \binom{n}{(n+j)/2} p^{(n+j)/2} (1-p)^{(n-j)/2}$ when $n + j$ is even; zero otherwise.
- Period $2$: you can only be back where you started at even times.
- **Recurrent iff $p = 1/2$** (symmetric walk). Otherwise transient — biased walks drift to $\pm\infty$.
- In $\mathbb Z^d$: recurrent for $d = 1, 2$; transient for $d \geq 3$ ("Pólya's theorem" / a drunk walks home, a drunk bird may not — not examinable in detail).

---

# 7 — Key Problem-Sheet Questions to Revise

The following questions are high-value: they either mirror exam-style problems, contain techniques that are reused, or cover the handful of trickier topics.

## Sigma algebras & probability (Q1 material)

- **PS1 Q1(a)–(e):** All the basic sigma-algebra existence/closure proofs. (a) and (d) especially are templates for Q1 in the exam.
- **PS1 Q2:** Sample space for repeated coin flips — good for understanding countable vs finite state spaces.
- **PS2 Q2(a):** Show $\mathcal F_X = \{X^{-1}(B) : B \in \mathcal B\}$ is a sigma algebra — exactly the proof pattern you'll need to apply to verify measurability.

## Random variables and transformations (Q2 material)

- **PS2 Q1:** Determining whether a function is a random variable w.r.t. a given sigma algebra. This is *the* type of question for the measurability part of Q1.
- **PS3 Q1:** Full suite of density transformations ($Y = X^4$, $e^X$, $\log X$, $(X - 0.5)^2$) — hits the monotonic + CDF methods and a non-monotonic case.
- **PS3 Q6:** Joint density, marginals, conditional densities, iterated expectation — a one-stop shop for the Q3 techniques.
- **PS4 Q7, PS4 Q8:** Order statistics derivations (min/max CDFs, and rescaling for convergence in distribution — direct precursor to 2024 Q2 and PS6 Q1).

## Bivariate / multivariate Normal (Q3 material)

- **PS4 Q1:** Given a bivariate density with a quadratic form in the exponent, read off the covariance matrix by matching coefficients — the core 2023 Q3 skill.
- **PS4 Q9:** Trivariate Normal conditional distribution — extends the bivariate technique.
- **PS6 Q5:** Cumulant generating function → mean and variance. Same computation as 2022 Q3.

## Convergence, CLT, MGFs (Q2 material)

- **PS4 Q8 / PS5 Q1:** Maximum of uniforms, convergence to 1 in probability, and then $n(1 - M_n) \xrightarrow{D} \text{Exp}(1)$. Appears in multiple exams re-packaged.
- **PS5 Q5:** $\mathbb E\!\left[\frac{|X_n|}{1 + |X_n|}\right] \to 0$ iff $X_n \xrightarrow{P} 0$ — lovely iff-style argument that consolidates the Markov-inequality technique.
- **PS5 Q6 (Slutsky):** Not examinable as theorem but great practice for convergence reasoning.
- **PS6 Q1, Q3:** Exponential to Gumbel via max + rescaling — directly the 2024 Q2 setup.

## Markov chains (Q4 material)

- **PS7 Q2:** Prove $P(n) = P^n$ by induction, and show that a stochastic matrix has eigenvalue $1$. **Both parts** get used in exams.
- **PS7 Q3:** Identify stochastic matrices, draw transition diagrams. Low-tech but tests fundamentals.
- **PS7 Q7:** An "almost-Markov" chain that satisfies Chapman–Kolmogorov but isn't Markov — good for understanding that Chapman–Kolmogorov is necessary but not sufficient.
- **PS8 Q1:** Computing multi-step distributions by matrix multiplication.
- **PS8 Q2 (gambler's ruin):** Set up and solve the hitting-probability linear system; directly relevant to the absorbing-boundary hitting-time calculations in 2021 Q4 and 2025 Q4.
- **PS8 Q3:** Periodicity of random walks on $\mathbb Z^d$ / cycles — classification-of-states practice.
- **PS8 Q4:** Birth-death chain hitting probability $h_i$ satisfying a second-order recurrence — the general technique behind any gambler's-ruin-style Q4.

---

# Appendix — Standard distributions quick reference

| Distribution | PDF/PMF | Mean | Var | MGF |
|---|---|---|---|---|
| Bernoulli$(p)$ | $p^x(1-p)^{1-x}$, $x \in \{0,1\}$ | $p$ | $p(1-p)$ | $1 - p + pe^t$ |
| Binomial$(n, p)$ | $\binom{n}{x} p^x(1-p)^{n-x}$ | $np$ | $np(1-p)$ | $(1 - p + pe^t)^n$ |
| Poisson$(\lambda)$ | $\frac{\lambda^x e^{-\lambda}}{x!}$ | $\lambda$ | $\lambda$ | $\exp(\lambda(e^t - 1))$ |
| Uniform$(a, b)$ | $\frac{1}{b-a}$ | $\frac{a+b}{2}$ | $\frac{(b-a)^2}{12}$ | $\frac{e^{tb} - e^{ta}}{t(b-a)}$ |
| Exp$(\lambda)$ | $\lambda e^{-\lambda x}$ | $1/\lambda$ | $1/\lambda^2$ | $\frac{\lambda}{\lambda - t}$, $t < \lambda$ |
| Gamma$(\alpha, \lambda)$ | $\frac{\lambda^\alpha x^{\alpha-1} e^{-\lambda x}}{\Gamma(\alpha)}$ | $\alpha/\lambda$ | $\alpha/\lambda^2$ | $\big(\frac{\lambda}{\lambda - t}\big)^\alpha$ |
| $N(\mu, \sigma^2)$ | $\frac{1}{\sqrt{2\pi\sigma^2}} e^{-(x-\mu)^2/2\sigma^2}$ | $\mu$ | $\sigma^2$ | $e^{\mu t + \sigma^2 t^2/2}$ |
| Beta$(\alpha, \beta)$ | $\frac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha,\beta)}$ | $\frac{\alpha}{\alpha+\beta}$ | $\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}$ | (no closed form) |
| Geometric$(p)$ | $(1-p)^{x-1} p$, $x \geq 1$ | $1/p$ | $(1-p)/p^2$ | $\frac{pe^t}{1 - (1-p)e^t}$, $t < -\log(1-p)$ |

*(Beware: Gamma/Exp can be parametrised by rate $\lambda$ or scale $\sigma = 1/\lambda$. 2025 Q3 used the scale version — always check which your problem is using.)*

---

**Final exam checklist** (things to have at your fingertips as you walk in):

- [ ] Three sigma-algebra axioms and how to use them to show closure under any reasonable set operation.
- [ ] Continuity property of probability measures and its proof.
- [ ] The "check intervals" shortcut for measurability.
- [ ] Three properties of a CDF.
- [ ] Monotonic-transformation formula with Jacobian and how to invert it.
- [ ] Probability integral transform.
- [ ] Complete-the-square identity for the bivariate Normal.
- [ ] Conditional distribution formula for bivariate Normal: $N(\rho y, 1 - \rho^2)$.
- [ ] Linear combination of MVN is MVN: $A\mathbf X \sim N(A\boldsymbol\mu, A\Sigma A^T)$.
- [ ] Markov's and Chebyshev's inequalities + proofs.
- [ ] WLLN + proof.
- [ ] Both Borel–Cantelli lemmas + the asymmetry in their conditions.
- [ ] MGF properties: affine, sum of independents, uniqueness, continuity.
- [ ] MGFs of $N(0,1)$, Exp$(\lambda)$, Gamma$(\alpha, \lambda)$, Poisson$(\lambda)$.
- [ ] CLT proof via MGFs.
- [ ] Chapman–Kolmogorov equations + proof.
- [ ] Classification: accessibility, communicating class, closed class, periodicity, recurrence/transience.
- [ ] $\sum_n p_{ii}(n) = \infty \iff i$ recurrent.
- [ ] Hitting probability linear system (with *minimal* solution caveat).
- [ ] Stationary distribution equations $\pi P = \pi$, $\sum \pi_i = 1$.

Good luck!