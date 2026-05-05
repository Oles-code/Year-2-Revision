# MATH50011 Statistical Modelling 1 — Exam Analysis & Revision Summary

**Based on: 5 past exams (2021–2025), 7 problem sheets, lecture notes, and all our chat history**

---

## Part A: Topics You've Found Difficult (from our conversations)

### 1. Confidence Interval Construction from Random Intervals
You asked for help on **2021 Q2(a)(ii)** (CI using $S_n^{\max}$) and **2023 Q3(d)** (verifying a random interval is a valid CI for Binomial). The struggle was with the general methodology: given a random interval $[L(X), U(X)]$, how do you verify it has coverage $\geq 1-\alpha$ for **all** $\theta$?

**The recipe you should internalise:**
1. Split over the possible values of $X$ using the law of total probability
2. For each value $X = x$, determine whether $\theta \in [L(x), U(x)]$ — this is a deterministic yes/no
3. Sum $P_\theta(X = x)$ over the "yes" cases to get coverage $P_\theta(\theta \in [L,U])$ as a function of $\theta$
4. Check this function is $\geq 1-\alpha$ for every $\theta$

### 2. The Two Forms of the F-Statistic (RSS form vs matrix quadratic form)
You were confused on **2021 Q4(b)** about how the mark scheme gets from the RSS-based $F = \frac{RSS_0 - RSS}{RSS}\cdot\frac{n-r}{r-s}$ to the matrix form $Q = \frac{(A\hat\beta)^T(A(X^TX)^{-1}A^T)^{-1}A\hat\beta}{\hat\sigma^2}$. These are equivalent — the matrix form arises from writing $H_0: A\beta = 0$ and using the fact that $A\hat\beta \sim N(0, \sigma^2 A(X^TX)^{-1}A^T)$ under $H_0$, then forming the standard quadratic form $W^T\Sigma^{-1}W \sim \chi^2$.

### 3. Dummy Variables and the Dummy Variable Trap
On **2021 Q4(a)** you asked why spring/spring×rain terms are excluded. The key: with $k$ categories and an intercept, you only include $k-1$ dummies. The omitted category becomes the baseline captured by $\beta_0$. Including all $k$ causes linear dependence (the dummies sum to the intercept column), making $X^TX$ singular. You can choose **any** category to drop — the fitted values are identical, only the interpretation of individual $\beta_i$ changes.

### 4. Multivariate Normal Independence
You worked through **Sheet 6 Q6** and needed guidance on the conditional distribution of a bivariate normal (**Sheet 6 Q7**). The core issue: for a jointly normal vector, uncorrelated $\Leftrightarrow$ independent (this is **not** true for general random vectors — Sheet 6 Q5 gives a counterexample). The exam question pattern is: given $Z = AX + b$ with $X \sim N(\mu, \Sigma)$, compute $\Sigma_Z = A\Sigma A^T$ and read off independence from the zero entries.

### 5. Fisher Information Conventions (per-observation vs full-sample)
On **2024 Q2(d)(ii)** you were confused about why the solutions divide the Fisher information by $n$. There are two conventions: $I_1(\theta)$ (single observation) and $I_n(\theta) = nI_1(\theta)$ (full sample). The asymptotic normality theorem is usually stated as $\sqrt{n}(\hat\theta - \theta) \xrightarrow{d} N(0, I_1(\theta)^{-1})$, so the variance of $\hat\theta$ itself is $[nI_1(\theta)]^{-1}$. When building a confidence region via the quadratic form $(\hat\theta - \theta)^T [\text{Cov}(\hat\theta)]^{-1} (\hat\theta - \theta)$, the $n$ reappears naturally.

### 6. Consistency via Var → 0 vs WLLN
On **2023 Q2(b)** you asked whether showing unbiased + $\text{Var}(T_n) \to 0$ is valid for proving consistency (instead of the WLLN approach in the solutions). **Yes, it is valid** — the lemma that "asymptotically unbiased + variance → 0 implies consistency" is a perfectly acceptable route, using Markov/Chebyshev inequality.

### 7. The Wilks Theorem / LRT Asymptotics
You asked for a full explanation of **Theorem 6** (Wilks): under $H_0$, $2\log t(Y) \xrightarrow{d} \chi^2_r$ where $r$ = number of independent restrictions. The key conceptual point is that this is an **asymptotic** result — it tells you the distribution of the LRT statistic without having to derive it exactly, but it requires the models to be nested ($\Theta_0 \subset \Theta$).

---

## Part B: Key Theorems (Precise Statements from the Notes)

### Theorem 1 — Cramér-Rao Lower Bound
For an unbiased estimator $T$ of $\theta$ (scalar, regular model):
$$\text{Var}_\theta(T) \geq \frac{1}{I(\theta)}, \qquad I(\theta) = -E_\theta\!\left[\frac{\partial^2}{\partial\theta^2}\log f_\theta(X)\right]$$
For iid samples: $I_n(\theta) = nI_1(\theta)$, so the bound is $1/(nI_1(\theta))$.

### Theorem 3 — The Delta Method
If $\sqrt{n}(T_n - \theta) \xrightarrow{d} N(0, \sigma^2(\theta))$ and $g$ is differentiable at $\theta$ with $g'(\theta) \neq 0$, then:
$$\sqrt{n}(g(T_n) - g(\theta)) \xrightarrow{d} N(0,\; g'(\theta)^2 \sigma^2(\theta))$$

### Theorem 6 — Wilks' Theorem (LRT Asymptotics)
Under $H_0$ (nested, with regularity conditions): $2\log t(Y) \xrightarrow{d} \chi^2_r$ as $n \to \infty$, where $r$ = number of independent restrictions defining $H_0$.

### Theorem 7 — The Gauss-Markov Theorem
Assume (FR) and (SOA). Let $c \in \mathbb{R}^p$ and let $\hat\beta$ be the LSE. Then $c^T\hat\beta$ has the **smallest variance** among all linear unbiased estimators of $c^T\beta$.

### Theorem 8 — Unbiasedness of $\hat\sigma^2$
$\hat\sigma^2 := \text{RSS}/(n-r)$ is an unbiased estimator of $\sigma^2$, where $r = \text{rank}(X)$.

### Lemma 11 — Projection Matrix Characterisation
$P$ is a projection matrix $\Leftrightarrow$ $P^T = P$ (symmetric) and $P^2 = P$ (idempotent).

### Lemma 14 — Independence in Multivariate Normal
If $Z \sim N(\mu, \Sigma)$ with $\Sigma$ block-diagonal (zero cross-covariances between blocks), then the corresponding sub-vectors of $Z$ are **independent**.

### Lemma 21 — Distribution of RSS
Under (FR) and (NTA): $\text{RSS}/\sigma^2 \sim \chi^2_{n-r}$ where $r = \text{rank}(X)$.

### Lemma 22 — t-pivot for $c^T\beta$
Under (FR) and (NTA): $\frac{c^T\hat\beta - c^T\beta}{\sqrt{c^T(X^TX)^{-1}c \cdot \text{RSS}/(n-p)}} \sim t_{n-p}$

### Lemma 23 — The F-test
Under $H_0: EY \in \text{span}(X_0)$:
$$F = \frac{RSS_0 - RSS}{RSS}\cdot\frac{n-r}{r-s} \sim F_{r-s,\, n-r}$$
where $r = \text{rank}(X)$, $s = \text{rank}(X_0)$.

### MLE Functional Invariance
If $\hat\theta$ is the MLE of $\theta$, then $g(\hat\theta)$ is the MLE of $g(\theta)$ for any function $g$.

### Slutsky's Lemma
If $X_n \xrightarrow{d} X$ and $Y_n \xrightarrow{P} c$ (constant), then $X_nY_n \xrightarrow{d} cX$, $X_n + Y_n \xrightarrow{d} X + c$, $X_n/Y_n \xrightarrow{d} X/c$ (if $c \neq 0$).

---

## Part C: Key Proofs / Recurring Proof Ideas

### Proofs that have been asked verbatim on exams:

| Proof | Asked in |
|-------|----------|
| Delta Method (full statement + proof) | 2024 Q1(d) |
| Gauss-Markov Theorem | 2021 Q3, 2024 Q3(c) |
| $\text{RSS}/\sigma^2 \sim \chi^2_{n-r}$ | 2023 Q3(c), 2025 (derive) |
| Markov's Inequality | 2024 Q2(c) |
| Projection matrix $\Leftrightarrow$ symmetric + idempotent | 2025 |
| MSE = Var + Bias² | Standard bookwork |
| Asymptotic normality $\Rightarrow$ consistency | 2024 Q1(b) |
| $\hat\sigma^2 = \text{RSS}/(n-r)$ is unbiased (Theorem 8) | Used everywhere |

### Recurring proof techniques:

**1. "Write it as a projection, then use Fisher-Cochran"** — This is how RSS/$\sigma^2 \sim \chi^2_{n-r}$ and the F-test are proved. The skeleton: write $Q = I - P$, show $Q$ is a projection matrix, then $Z^TQZ \sim \chi^2_{\text{rank}(Q)}$ by Lemma 18, and independence of numerator/denominator by Lemma 19 ($A_1A_2 = 0$).

**2. "Unbiasedness of any competitor minus $c^T\hat\beta$ forces $D^TX = 0$"** — The Gauss-Markov proof trick. If $\hat\gamma = L^TY$ is any other linear unbiased estimator, define $D^T = L^T - c^T(X^TX)^{-1}X^T$. Then $\text{Var}(\hat\gamma) = \text{Var}(c^T\hat\beta) + \text{Var}(D^TY) + 2\text{Cov}(\cdots)$. The cross-covariance vanishes because $D^TX = 0$ (forced by unbiasedness of both estimators). So $\text{Var}(\hat\gamma) \geq \text{Var}(c^T\hat\beta)$.

**3. "Taylor expand + Slutsky"** — The Delta Method. Taylor expand $g(T_n)$ around $\theta$, multiply by $\sqrt{n}$, use consistency to handle the remainder, apply Slutsky.

**4. "Markov/Chebyshev → consistency"** — For $P(|T_n - \theta| \geq \varepsilon) \leq \text{MSE}/\varepsilon^2 \to 0$. This route works whenever you can show bias → 0 and variance → 0.

---

## Part D: Exam Frequency Table (2021–2025)

| Topic | '21 | '22 | '23 | '24 | '25 | Frequency |
|-------|:---:|:---:|:---:|:---:|:---:|:---------:|
| MLE computation + invariance | ✓ | ✓ | ✓ | ✓ | ✓ | **5/5** |
| CI construction + test-from-CI | ✓ | ✓ | ✓ | ✓ | ✓ | **5/5** |
| Linear model matrix setup + LSE | ✓ | ✓ | ✓ | ✓ | ✓ | **5/5** |
| Asymptotic normality / Delta method | ✓ | ✓ | — | ✓* | ✓ | **4/5** |
| Hypothesis testing / power function | ✓ | ✓ | ✓ | ✓ | — | **4/5** |
| Numerical LSE from summary statistics | — | ✓ | ✓ | ✓ | ✓ | **4/5** |
| Consistency (definition/proof/use) | — | ✓ | ✓ | ✓* | ✓ | **4/5** |
| Gauss-Markov (stated/proved) | ✓ | ✓ | — | ✓* | — | **3/5** |
| Multivariate Normal + independence | ✓ | — | — | ✓ | ✓ | **3/5** |
| RSS distribution proof/use | — | — | ✓* | ✓ | ✓* | **3/5** |
| F-test | ✓ | — | ✓ | — | ✓ | **3/5** |
| Identifiability / dummy trap | — | ✓ | — | ✓ | ✓ | **3/5** |
| Confidence regions (ellipse/F) | — | ✓ | — | ✓ | — | 2/5 |
| Reparametrisation | — | ✓ | ✓ | — | — | 2/5 |
| R² / adjusted R² | — | ✓ | — | ✓ | — | 2/5 |
| LRT (Wilks) | — | — | ✓ | — | — | 1/5 |
| Projection matrices (proof) | — | — | — | — | ✓* | 1/5 |
| WLS | — | ✓ | — | — | — | 1/5 |

*starred = asked to prove, not just use*

---

## Part E: Common & Difficult Exam Question Types

### TYPE 1: "Compute the MLE, then use functional invariance / Delta method for a function of $\theta$" (Every year)

**Example pattern (2024 Q2(d)):** Given $(Y_i, X_i)$ with $Y_i \sim \text{Poisson}(\lambda)$, $X_i | Y_i \sim \text{Binomial}(Y_i, p)$. Find the MLE of $(\lambda, p)$, then construct an asymptotic confidence region using the Fisher information matrix.

**What makes it hard:** You need to compute the Fisher information matrix (possibly multivariate), recognise the quadratic form $n(\hat\theta - \theta)^T I_1(\theta)(\hat\theta - \theta) \xrightarrow{d} \chi^2_k$, and connect this to a confidence region. The per-observation vs full-sample convention trips people up.

### TYPE 2: "Given summary statistics, compute $\hat\beta$, RSS, $R^2$, then do a t-test" (2022, 2023, 2024, 2025)

**Example pattern (2024 Q4):** Model $Y_i = \beta_1 + w_i\beta_2 + \varepsilon_i$ with $n = 83$. Given $\frac{1}{n}\sum y_i = 2$, $\frac{1}{n}\sum w_i = 2$, $\frac{1}{n}\sum w_i^2 = 5$, $\frac{1}{n}\sum w_iy_i = 5$, $\frac{1}{n}\sum y_i^2 = 30$. Show $\hat\beta = (0, 1)^T$, compute $\hat\sigma^2$, $R^2$, and do a one-sided t-test for $\beta_2$.

**The mechanical recipe:**
1. Normal equations: $X^TX\hat\beta = X^TY$. For simple linear regression:
$$X^TX = n\begin{pmatrix} 1 & \bar{w} \\ \bar{w} & \overline{w^2}\end{pmatrix}, \quad X^TY = n\begin{pmatrix} \bar{y} \\ \overline{wy}\end{pmatrix}$$
2. Solve $2\times 2$ system.
3. $\text{RSS} = \sum y_i^2 - \hat\beta^T X^T Y$ (useful shortcut: $\text{RSS} = Y^TY - \hat{Y}^T\hat{Y}$).
4. $\hat\sigma^2 = \text{RSS}/(n-p)$.
5. $R^2 = 1 - \text{RSS}/\text{TSS}$ where $\text{TSS} = \sum(y_i - \bar{y})^2 = n(\overline{y^2} - \bar{y}^2)$.

### TYPE 3: "Which components of $Z = AX + b$ are independent?" (2024 Q3(a), 2025 Q4(a), Sheet 6 Q6)

**Example pattern (2024 Q3(a)):** $X \sim N\left(\binom{2}{3}, I_2\right)$, $Z = \begin{pmatrix}1&1\\0&1\\1&0\end{pmatrix}X + \begin{pmatrix}-1\\-3\\2\end{pmatrix}$.

**Recipe:**
1. $\mu_Z = A\mu + b$
2. $\Sigma_Z = A\Sigma_X A^T$ (here $= AA^T$ since $\Sigma_X = I$)
3. $Z_i \perp Z_j$ iff $(\Sigma_Z)_{ij} = 0$

This is **nearly identical** across 2024, 2025, and Sheet 6 Q6. Practise Sheet 6 Q6 and you've essentially done the exam version.

### TYPE 4: "Construct a test from a confidence interval" (2021, 2023, 2024)

**Example pattern (2024 Q3(d)):** $X \sim \text{Bernoulli}(\theta)$. CI is $[0, 1-\alpha]$ when $X = 0$ and $[\alpha, 1]$ when $X = 1$. Construct a level-$\alpha$ test for $H_0: \theta = \theta_0$ and draw the power function.

**The principle:** Reject $H_0: \theta = \theta_0$ when $\theta_0 \notin [L(X), U(X)]$. The power function is $\beta(\theta) = P_\theta(\text{reject}) = P_\theta(\theta_0 \notin [L, U])$, which you compute by splitting over the possible $X$ values.

**What makes it hard:** The power function has different formulas in different regions of $\theta$ (often 3 cases), and drawing the graph correctly requires care.

### TYPE 5: "Prove the Gauss-Markov theorem" (2021, 2024)

This is a bookwork proof asked repeatedly. The trick is the decomposition $\hat\gamma = c^T\hat\beta + D^TY$ where unbiasedness of both forces $D^TX = 0$, which kills the cross-covariance term.

### TYPE 6: "Prove or derive the RSS/$\sigma^2 \sim \chi^2_{n-r}$ result" (2023, 2024, 2025)

Write $\text{RSS} = Y^TQY$ where $Q = I - P$. Then $\text{RSS}/\sigma^2 = Z^TQZ$ where $Z = Y/\sigma \sim N(X\beta/\sigma, I)$. Since $Q$ is a projection matrix of rank $n-r$, apply Lemma 18. The non-centrality is zero because $QX = 0$.

### TYPE 7: "Categorical predictors / dummy variables / F-test for sub-model" (2021 Q4, 2023 Q4, 2025 Q3)

**Example pattern (2021 Q4):** Seasonal model with SUMMER, AUTUMN, WINTER dummies + RAIN + interactions. Test $H_0$: all seasonal/interaction effects are zero using the F-test. Compute degrees of freedom, relate RSS to $\hat\sigma^2$.

**Key facts:**
- With $k$ categories + intercept, use $k-1$ dummies
- $\text{RSS} = (n-r)\hat\sigma^2$ where $r = \text{rank}(X)$
- F-test degrees of freedom: numerator = $r - s$ (parameters removed), denominator = $n - r$

### TYPE 8: "Delta Method proof" (2024 Q1(d))

Full statement-and-proof. The skeleton: Taylor expand → rearrange → consistency of $T_n$ implies $\tilde\theta_n \xrightarrow{P} \theta$ → continuity gives $g'(\tilde\theta_n) \xrightarrow{P} g'(\theta)$ → Slutsky.

---

## Part F: High-Priority Problem Sheet Questions

These are the sheet questions that most directly map to exam questions:

| Sheet Question | Maps to Exam Question(s) | Why |
|---------------|--------------------------|-----|
| **Sheet 6 Q6** | 2024 Q3(a), 2025 Q4(a) | Nearly verbatim "which components independent?" |
| **Sheet 6 Q1** | 2023 Q4, 2025 Q3 | 3-group dummy model + reparametrisation + fitted values |
| **Sheet 5 Q1** | 2022 Q3, 2024 Q4 | Dummy variable model for two-group comparison |
| **Sheet 5 Q5** | Bedrock for all LSE questions | Derivation of simple linear regression LSE from scratch |
| **Sheet 4 Q8** | 2022 Q2, 2023 Q1 | Two-sample t-test derived as LRT |
| **Sheet 4 Q9** | 2023 Q1(c) | LRT for Geometric — nearly identical |
| **Sheet 4 Q3–4** | 2024 Q4 | One-sample t-test, CI/test duality |
| **Sheet 4 Q5** | 2022 Q2(d) | Power function & sample-size calculation |
| **Sheet 7 Q4** | Conceptual | Equivalence of t-test and F-test for single coefficient |
| **Sheet 7 Q5** | 2021 Q4, 2024 Q4 | Interaction terms, testing main effects vs interactions |
| **Sheet 6 Q2** | Foundation | Computing $\text{Cov}(e)$ from $Q = I - P$ |
| **Sheet 6 Q7** | Challenging | Conditional distribution of bivariate Normal |
| **Sheet 7 Q6** | 2022 Q4 | Under/over-fitting bias and variance analysis |

---

## Part G: Quick-Reference Cheat Sheet

- $I_n(\theta) = nI_1(\theta)$ — Fisher information scales linearly with $n$
- MLE asymptotic variance: $1/(nI_1(\theta))$
- $\hat\beta = (X^TX)^{-1}X^TY$, requires $X$ full column rank
- $E(\hat\beta) = \beta$, $\text{Cov}(\hat\beta) = \sigma^2(X^TX)^{-1}$
- $\hat\sigma^2 = \text{RSS}/(n-r)$ is unbiased; $r = \text{rank}(X)$
- Under NTA: $\hat\beta \sim N(\beta, \sigma^2(X^TX)^{-1})$, $\text{RSS}/\sigma^2 \sim \chi^2_{n-r}$, and they are **independent**
- $P = X(X^TX)^{-1}X^T$, $\text{trace}(P) = \text{rank}(P) = r$
- Projection $\Leftrightarrow$ symmetric + idempotent
- Jointly Normal + uncorrelated $\Rightarrow$ independent (NOT true in general)
- Gauss-Markov needs only SOA + FR; normality is NOT required
- Normality IS required for exact $t$ and $F$ distributions
- Adding columns to $X$ never increases RSS, never decreases $R^2$
- $z_{0.025} \approx 1.96$, $z_{0.05} \approx 1.645$, $z_{0.10} \approx 1.28$
- Asymptotic normality $\Rightarrow$ consistency (via Slutsky: $T_n - \theta = \frac{1}{\sqrt{n}} \cdot \sqrt{n}(T_n - \theta) \xrightarrow{P} 0$)
- Unbiased does NOT imply $h(T)$ unbiased for $h(\theta)$ (unless $h$ is linear)