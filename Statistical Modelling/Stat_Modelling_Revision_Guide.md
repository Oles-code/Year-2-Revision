# Statistical Modelling 1 — Revision Summary

**MATH50011 · Imperial College London**
*Focused on topics that recur in exams, especially in the places people lose marks.*

---

## How to read this document

I've gone through all five past papers (2021–2025) and the seven problem sheets and identified what actually comes up. The structure below follows the course, but the **emphasis** is deliberately tilted:

- Topics that appear on **every exam** (MLE, asymptotic normality, LSE, CI construction, hypothesis testing) get more space and more worked intuition.
- **Bookwork proofs that have been asked verbatim** (Gauss-Markov, $\text{RSS}/\sigma^2 \sim \chi^2_{n-r}$, Delta method, Markov's inequality, projection matrix iff symmetric + idempotent) are given in full.
- The recurring **numerical patterns** (compute $\hat\beta$ from a list of $\frac{1}{n}\sum$ values; figure out independence of a linear transformation of a Normal vector) get their own worked templates.
- A final "problem sheet checklist" points at the specific questions I'd redo if I had limited time.

I use the notation from the lecture notes throughout ($\theta$ for a parameter, $T$ for an estimator, $\beta$ for the regression vector, etc.).

---

## Part I — Point estimation

### 1. The three properties of an estimator

You have a model $\{P_\theta : \theta \in \Theta\}$, observe data from $P_{\theta_0}$ (where $\theta_0$ is "the truth"), and want a statistic $T = t(Y_1, \dots, Y_n)$ that is close to $\theta_0$. Three numbers measure "close":

- **Bias:** $\text{bias}_\theta(T) = \mathbb{E}_\theta(T) - \theta$. Zero bias ⇔ unbiased.
- **Variance (standard error):** $\text{Var}_\theta(T)$; $\text{SE}_\theta(T) = \sqrt{\text{Var}_\theta(T)}$.
- **MSE:** $\text{MSE}_\theta(T) = \mathbb{E}_\theta[(T - \theta)^2]$.

The single most useful identity in the course:

$$\boxed{\text{MSE}_\theta(T) = \text{Var}_\theta(T) + \text{bias}_\theta(T)^2}$$

**Proof** (worth knowing — it's been asked in 2025). Let $\mu = \mathbb{E}_\theta T$. Then
$$
\mathbb{E}_\theta[(T-\theta)^2] = \mathbb{E}_\theta[(T - \mu + \mu - \theta)^2] = \mathbb{E}_\theta[(T-\mu)^2] + 2(\mu - \theta)\underbrace{\mathbb{E}_\theta(T-\mu)}_{=0} + (\mu-\theta)^2
$$
and the three terms are $\text{Var}_\theta(T)$, $0$, $\text{bias}_\theta(T)^2$. $\square$

**Conceptual take-away.** Bias tells you where your estimator is centred; variance tells you how wobbly it is. MSE packages both into one number, and the trade-off between them (you can sometimes *add a little bias* to get *less variance* overall) is the motivation for half the model-selection material at the end of the course.

> **Worked example of the trade-off.** With $X \sim \text{Bin}(n, p)$, compare $S = X/n$ (unbiased) with $T = (X+1)/(n+2)$ (biased). When $p$ is close to $1/2$, $T$ has smaller MSE; when $p$ is near 0 or 1, $S$ wins. So "best" depends on the truth — this is the whole reason we need **Cramér-Rao** and other optimality notions.

### 2. A small but useful fact

$T$ unbiased for $\theta$ **does not** imply $h(T)$ unbiased for $h(\theta)$. Counterexample: $Y_i \sim N(\mu, \sigma^2)$, $s^2$ is unbiased for $\sigma^2$ but $s$ is **not** unbiased for $\sigma$ (Problem Sheet 1 Q5b). This mistake is a classic exam trap.

Exception: if $h$ is *linear*, $h(T) = aT + b$, then $\mathbb{E}(aT + b) = a\theta + b = h(\theta)$ so unbiasedness is preserved. (2024 Q2(a) asked exactly this.)

### 3. Cramér-Rao lower bound

For an unbiased estimator $T$ of $\theta$ (scalar, regular model):
$$
\text{Var}_\theta(T) \ge \frac{1}{I(\theta)}, \qquad I(\theta) = \mathbb{E}_\theta\!\left[\left(\frac{\partial}{\partial \theta}\log f_\theta(X)\right)^{\!2}\right] = -\mathbb{E}_\theta\!\left[\frac{\partial^2}{\partial\theta^2}\log f_\theta(X)\right].
$$

**Two things to remember:**

1. **For iid samples:** $I_n(\theta) = n\,I_1(\theta)$ — Fisher information scales linearly with $n$, so the lower bound on variance scales like $1/n$. This is why "$\sqrt{n}$-consistent" keeps appearing.
2. **Practical recipe.** To find the CRLB for an iid sample from $f_\theta$:
   - compute $\ell = \log f_\theta(x)$;
   - compute $\partial^2 \ell / \partial\theta^2$;
   - take $\mathbb{E}_\theta$ and negate to get $I_1(\theta)$;
   - lower bound is $\frac{1}{n\,I_1(\theta)}$.

The full proof (Cauchy-Schwarz trick) is in the notes; it's long and not usually asked to reproduce.

> **The CRLB is a floor, not a guarantee.** There doesn't have to exist an unbiased estimator that *achieves* it. When one does exist (e.g. $\bar{X}$ for the mean of a Bernoulli), we call it **efficient**.

---

## Part II — Asymptotic tools (these pervade the whole course)

The exam uses five results constantly. Know their **statements** and when to apply each.

### Convergence types

- **In probability** $T_n \xrightarrow{P} c$: for all $\varepsilon > 0$, $P(|T_n - c| > \varepsilon) \to 0$. "Eventually close with high probability."
- **In distribution** $T_n \xrightarrow{d} T$: $F_{T_n}(x) \to F_T(x)$ at continuity points. "Shape of the distribution settles down."

Weak law of large numbers: $\bar{X}_n \xrightarrow{P} \mu$. Central limit theorem: $\sqrt{n}(\bar{X}_n - \mu) \xrightarrow{d} N(0, \sigma^2)$.

### Consistency

$T_n$ is **consistent** for $g(\theta)$ if $T_n \xrightarrow{P_\theta} g(\theta)$ for every $\theta$.

The cheap route to consistency:

**Lemma.** If $T_n$ is **asymptotically unbiased** ($\mathbb{E}_\theta T_n \to g(\theta)$) and $\text{Var}_\theta T_n \to 0$, then $T_n$ is consistent.

**Proof** (via Markov's inequality — also asked as standalone bookwork in 2024 Q2(c)):

*Markov's inequality.* For $X \ge 0$ and $a > 0$: $P(X \ge a) \le \mathbb{E}(X)/a$. To see this, $a\,\mathbb{1}\{X \ge a\} \le X$, so taking expectations, $a\,P(X \ge a) \le \mathbb{E}(X)$. $\square$

*Applying it to $|T_n - g(\theta)|$:*
$$
P_\theta(|T_n - g(\theta)| \ge \varepsilon) \le \frac{\mathbb{E}_\theta[(T_n - g(\theta))^2]}{\varepsilon^2} = \frac{\text{MSE}}{\varepsilon^2} = \frac{\text{Var} + \text{bias}^2}{\varepsilon^2} \to 0. \square
$$

> **Warning: the converse is false.** A consistent estimator need not be asymptotically unbiased, and an asymptotically unbiased estimator need not be consistent (Problem Sheet 1 Q11, Q12). Don't mix these up in an exam.

### Asymptotic normality

$T_n$ is **asymptotically normal** for $\theta$ if
$$\sqrt{n}(T_n - \theta) \xrightarrow{d} N(0, \sigma^2(\theta)).$$

**Asymptotic normality $\Rightarrow$ consistency** (2024 Q1(b) asked for this proof): if $\sqrt{n}(T_n-\theta)$ has a limit distribution, then $T_n - \theta = \frac{1}{\sqrt{n}} \cdot W_n$ where $W_n \xrightarrow{d} W$. Since $\frac{1}{\sqrt{n}} \to 0$, by Slutsky (see below) $T_n - \theta \xrightarrow{d} 0$, and convergence in distribution to a constant implies convergence in probability. $\square$

### Slutsky's lemma

If $X_n \xrightarrow{d} X$ and $Y_n \xrightarrow{P} c$ (constant), then:
$$X_n + Y_n \xrightarrow{d} X + c, \quad X_n Y_n \xrightarrow{d} cX, \quad X_n/Y_n \xrightarrow{d} X/c \text{ if } c \ne 0.$$

This is the "replace a constant by a consistent estimator of it without breaking anything" lemma. Used heavily for:
- Converting $\sigma(\theta)$ to $\hat\sigma_n$ in asymptotic CIs (so you get a usable pivot).
- Proving that the $t$-statistic $\sqrt{n}(\bar X - \mu)/S_n \xrightarrow{d} N(0,1)$.

### The Delta method

If $\sqrt{n}(T_n - \theta) \xrightarrow{d} N(0, \sigma^2(\theta))$ and $g$ is differentiable at $\theta$ with $g'(\theta) \ne 0$, then
$$\sqrt{n}(g(T_n) - g(\theta)) \xrightarrow{d} N(0, g'(\theta)^2 \sigma^2(\theta)).$$

**This was asked as a full statement-and-proof in 2024.** The proof is short:

By Taylor, $g(T_n) = g(\theta) + g'(\tilde\theta)(T_n - \theta)$ for some $\tilde\theta$ between $T_n$ and $\theta$. Since $T_n \xrightarrow{P} \theta$ (by the same asymptotic-normality-⇒-consistency argument above), $\tilde\theta \xrightarrow{P} \theta$, and by continuity of $g'$, $g'(\tilde\theta) \xrightarrow{P} g'(\theta)$. Then
$$\sqrt{n}(g(T_n) - g(\theta)) = g'(\tilde\theta) \cdot \sqrt{n}(T_n - \theta)$$
and Slutsky gives $g'(\theta) \cdot N(0, \sigma^2(\theta)) = N(0, g'(\theta)^2 \sigma^2(\theta))$. $\square$

> **Why you care.** You usually know the asymptotic distribution of a simple thing ($\bar X$), but you want the asymptotic distribution of a transformed quantity (odds = $p/(1-p)$, log mean, etc.). The Delta method is the lever. 2021 Q1(b) used it for $\Phi(-\hat\mu/\sigma)$; Problem Sheet 3 Q6 uses it for odds; Problem Sheet 2 Q3 for the geometric mean.

### Continuous mapping theorem

If $X_n \xrightarrow{d} X$ (or $\xrightarrow{P}$, or almost surely) and $g$ is continuous, then $g(X_n) \xrightarrow{d} g(X)$ (resp. the other modes). Use this when you need to apply a continuous function inside a convergence statement.

---

## Part III — Maximum likelihood estimation (MLE)

### Recipe

Given data $Y$ with joint density $f_\theta(y)$:
1. Write the **likelihood** $L(\theta) = f_\theta(Y)$ (as a function of $\theta$, with data fixed).
2. Take the **log-likelihood** $\ell(\theta) = \log L(\theta)$ — for iid data this turns the product into a sum.
3. Solve $\partial\ell/\partial\theta = 0$, check it's a maximum (second derivative negative).

For iid $Y_1, \dots, Y_n \sim f^{(1)}_\theta$: $\ell(\theta) = \sum_{i=1}^n \log f^{(1)}_\theta(Y_i)$.

### Key properties

**(a) Functional invariance.** If $\hat\theta$ is the MLE of $\theta$ and $g$ is a bijection, then $\widehat{g(\theta)} = g(\hat\theta)$ is the MLE of $g(\theta)$. (Extended to non-bijective $g$ via induced likelihood.)

This is powerful: rather than re-maximising, you just transform. Example from 2021 Q1(b): to find the MLE of $\theta = P(Y_i \le 0) = \Phi(-\mu/\sigma)$ in a $N(\mu, \sigma^2)$ model with $\sigma^2$ known, just take $\hat\theta_Y = \Phi(-\hat\mu/\sigma) = \Phi(-\bar Y / \sigma)$. No re-derivation needed.

**(b) Not necessarily unbiased.** For $Y_i \sim N(\mu, \sigma^2)$ iid, the MLE of $\sigma^2$ is $\frac{1}{n}\sum(Y_i - \bar Y)^2$, biased by a factor $(n-1)/n$. Consistency is what MLEs *do* have, not unbiasedness.

**(c) Asymptotic distribution (big result).** Under regularity,
$$\sqrt{n}(\hat\theta_n - \theta_0) \xrightarrow{d} N\!\left(0, I_f(\theta_0)^{-1}\right).$$

Meaning: the MLE is asymptotically unbiased, is $\sqrt{n}$-consistent, and **achieves the Cramér-Rao lower bound** asymptotically. So MLEs are "as good as it gets" for large $n$.

**Multivariate version:** $\sqrt{n}(\hat\theta_n - \theta_0) \xrightarrow{d} N(0, I_f(\theta_0)^{-1})$ where the information matrix is $I_f(\theta) = -\mathbb{E}_\theta[\nabla^2 \log f(X;\theta)]$.

### How to estimate $I_f(\theta_0)$ from the data

In practice we don't know $\theta_0$, so we can't use $I_f(\theta_0)$ directly. Three consistent plug-in estimators:
1. $I_f(\hat\theta)$ — plug MLE into the formula.
2. $\frac{1}{n}\sum_i \left[\frac{\partial}{\partial\theta}\log f(x_i;\theta)\right]^2_{\theta=\hat\theta}$ — "observed information, outer product version".
3. $-\frac{1}{n}\sum_i \frac{\partial^2}{\partial\theta^2} \log f(x_i;\theta)\big|_{\theta=\hat\theta}$ — "observed information, Hessian version".

> **Exam checklist for any "find the MLE" question:**
> 1. Write down the log-likelihood for *your* data.
> 2. Differentiate, set to zero, solve.
> 3. State the MLE.
> 4. If asked for its distribution: either exact (e.g. $\bar Y \sim N(\mu, \sigma^2/n)$ for Normal data) or asymptotic ($\sqrt{n}(\hat\theta - \theta) \xrightarrow{d} N(0, I_f(\theta)^{-1})$).
> 5. If asked for the MLE of a transformed parameter $g(\theta)$: invariance + Delta method.

### Identifiability

A model is **identifiable** if $\theta_1 \ne \theta_2 \Rightarrow P_{\theta_1} \ne P_{\theta_2}$. Without identifiability, the MLE isn't even well-defined (multiple $\theta$ give the same likelihood). Example from 2024 Q2(b): the model $N(\theta^\alpha, 1)$ is identifiable iff $\alpha$ is **odd**, because $\theta \mapsto \theta^\alpha$ must be injective on $\mathbb{R}$.

In linear models, the analogue is the **full rank (FR) assumption** — see Part VI.

---

## Part IV — Confidence intervals

### The pivotal quantity recipe

A **pivot** is a function $t(Y, \theta)$ of data and parameter whose distribution does **not** depend on $\theta$ (or any nuisance parameter).

Once you have a pivot, CIs are easy:
1. Find $a_1, a_2$ with $P(a_1 \le t(Y, \theta) \le a_2) \ge 1-\alpha$ (using the known distribution).
2. Solve the inequalities for $\theta$ to get $h_1(Y) \le \theta \le h_2(Y)$.
3. $[h_1(y), h_2(y)]$ is a $(1-\alpha)$ CI.

### The two big standard cases

**(i) $Y_i \sim N(\mu, \sigma^2)$, $\sigma^2$ known:** Pivot is $\frac{\bar Y - \mu}{\sigma/\sqrt n} \sim N(0,1)$, giving CI $\bar y \pm z_{\alpha/2} \sigma/\sqrt{n}$.

**(ii) $Y_i \sim N(\mu, \sigma^2)$, both unknown:** Pivot is $\frac{\bar Y - \mu}{S/\sqrt n} \sim t_{n-1}$, giving CI $\bar y \pm t_{n-1, \alpha/2} s/\sqrt{n}$.

For $\sigma^2$, use $\frac{\sum(Y_i - \bar Y)^2}{\sigma^2} \sim \chi^2_{n-1}$ as the pivot.

### Asymptotic CIs from asymptotic normality

If $\sqrt{n}(T_n - \theta) \xrightarrow{d} N(0, \sigma^2(\theta))$ and $\hat\sigma_n \xrightarrow{P} \sigma(\theta)$ (e.g. $\hat\sigma_n = \sqrt{\sigma^2(\hat\theta)}$), then by Slutsky
$$\sqrt{n}\cdot\frac{T_n - \theta}{\hat\sigma_n} \xrightarrow{d} N(0,1),$$
giving the **approximate** $(1-\alpha)$ CI
$$T_n \pm z_{\alpha/2} \cdot \frac{\hat\sigma_n}{\sqrt{n}} = T_n \pm z_{\alpha/2}\cdot \widehat{\text{SE}}(T_n).$$

This is the workhorse for everything that isn't exactly Normal. It's what gets used for Binomial proportions, odds ratios, sample medians, you name it.

### Simultaneous CIs — Bonferroni

To build a joint CI for $(\theta_1, \dots, \theta_k)$: make a $(1 - \alpha/k)$ CI for each $\theta_i$ separately; the Cartesian product of these is a $(1-\alpha)$ simultaneous CI. Proof is a one-line union bound.

**Warning:** Bonferroni is **conservative**. If the intervals are actually independent, the true coverage is $(1-\alpha/k)^k > 1-\alpha$. For linear models we often have a better (ellipsoidal) option — see Part VII.

### Building a test from a CI (very common exam question)

**Test recipe.** Given a $(1-\alpha)$ CI/region $A(Y)$ for $\theta$, to test $H_0: \theta \in \Theta_0$:
> Reject $H_0$ iff $\Theta_0 \cap A(y) = \emptyset$.

For a simple null $H_0: \theta = \theta_0$, this reduces to "reject iff $\theta_0 \notin A(y)$".

**Why it has level $\alpha$:** for any $\theta \in \Theta_0$,
$$P_\theta(\text{reject}) = P_\theta(\Theta_0 \cap A(Y) = \emptyset) \le P_\theta(\theta \notin A(Y)) \le \alpha. \square$$

This pops up in 2023 Q1(d), 2024 Q3(d), 2025 Q1(c). You should be able to write the proof from memory.

---

## Part V — Hypothesis testing

### Key vocabulary

- **$H_0$ vs $H_1$**, **Type I** (reject true $H_0$), **Type II** (accept false $H_0$).
- **Level $\alpha$:** $\sup_{\theta \in \Theta_0} P_\theta(\text{reject}) \le \alpha$.
- **Power function:** $\beta(\theta) = P_\theta(\text{reject } H_0)$. You want it **small on $\Theta_0$**, **large on $\Theta_1$**.
- **p-value:** $p = \sup_{\theta \in \Theta_0} P_\theta(T \ge t_{\text{obs}})$ (for a right-tailed test). Reject iff $p \le \alpha$.

### Why tests are symmetric around the null for two-sided alternatives

A question from 2022: "why do we pick $\pm t_{n-1, \alpha/2}$ symmetrically?" The answer: symmetry of the $t$-distribution (or Normal) around zero means symmetric critical values give the **shortest** rejection region at a fixed level, and equivalently the **narrowest** confidence interval. Other splits (e.g. $\alpha/4$ left, $3\alpha/4$ right) would give a valid but wider interval.

### Power function computation template

With $H_0: \mu = \mu_0$ vs $H_1: \mu \ne \mu_0$ and test statistic $T = \sqrt{n}(\bar X - \mu_0)$, under $N(\mu, 1)$ data:

Under general $\mu$, $T \sim N(\sqrt{n}(\mu - \mu_0), 1)$. So
$$\beta(\mu) = P_\mu(|T| \ge z_{\alpha/2}) = 1 - \Phi(z_{\alpha/2} - \sqrt{n}(\mu - \mu_0)) + \Phi(-z_{\alpha/2} - \sqrt{n}(\mu - \mu_0)).$$

Sketch (the classic symmetric bathtub):

```
β(µ)
 1 |                 __....——....__
   |              _.·              ·._
   |           .·                      ·.
 α +--------·----------.----------·--------
   |     .·           µ₀           ·.
 0 |__·__________________________________·__
           ←  reject →   ←  don't  →   ← reject →
```

Power is $\alpha$ at $\mu_0$, rises symmetrically as $\mu$ moves away from $\mu_0$, approaches 1 far out.

### The Likelihood Ratio Test (LRT)

The LRT statistic is
$$t(y) = \frac{\sup_{\theta \in \Theta} L(\theta; y)}{\sup_{\theta \in \Theta_0} L(\theta; y)} \ge 1.$$

**Wilks' theorem (big result):** Under $H_0$ and mild regularity (including that $\Theta_0$ is a lower-dimensional subspace of $\Theta$),
$$2\log t(Y) \xrightarrow{d} \chi^2_r,$$
where $r$ = (dimension of $\Theta$) − (dimension of $\Theta_0$) = number of independent restrictions imposed by $H_0$.

**Sketch of the proof** (for simple null $\Theta_0 = \{\theta_0\}$):
By Taylor around $\hat\theta$ (where $\partial \log L / \partial\theta = 0$):
$$\log L(\theta_0) \approx \log L(\hat\theta) + \tfrac{1}{2}(\theta_0 - \hat\theta)^\top \!\left[-\frac{\partial^2 \log L}{\partial\theta\partial\theta^\top}\bigg|_{\hat\theta}\right] (\theta_0 - \hat\theta).$$
So $2\log t(Y) \approx (\hat\theta - \theta_0)^\top I_f(\theta_0) (\hat\theta - \theta_0) \cdot n$. Since $\sqrt{n}(\hat\theta - \theta_0) \xrightarrow{d} Z \sim N(0, I_f^{-1})$, this is $Z^\top I_f Z$, which by distribution theory (see Part VIII) is $\chi^2_r$ where $r = \dim \Theta$. $\square$

> **When to reach for the LRT.** Two settings where it's the go-to: (i) you want to compare nested models (restricted vs unrestricted), and (ii) simpler pivots aren't obvious. Problem Sheet 4 Q8 derives the two-sample $t$-test as an LRT — a lovely application that ties the material together.

---

## Part VI — Linear models (Second Order Assumptions)

### Setup

Data: $Y \in \mathbb{R}^n$. Model: $Y = X\beta + \varepsilon$, where
- $X \in \mathbb{R}^{n \times p}$ is the **design matrix** (known, fixed);
- $\beta \in \mathbb{R}^p$ is the unknown parameter vector;
- $\varepsilon$ is an $n$-vector of unobservable errors with $\mathbb{E}\varepsilon = 0$.

**Three key assumptions** (know how to state each — asked in 2024 Q3(b), 2025 Q3(a)):
- **(FR)** Full rank: $\text{rank}(X) = p$.
- **(SOA)** Second-order assumption: $\text{Cov}(\varepsilon) = \sigma^2 I_n$. (Errors uncorrelated, equal variance.)
- **(NTA)** Normal theory assumption: $\varepsilon \sim N(0, \sigma^2 I_n)$. Implies SOA.

### Least-squares estimator (LSE)

Minimise $S(\beta) = \|Y - X\beta\|^2$. Setting $\partial S/\partial\beta = -2X^\top Y + 2X^\top X\beta = 0$:

$$\boxed{X^\top X \hat\beta = X^\top Y \quad\text{(normal equations)}}$$

Under (FR), $X^\top X$ is invertible (Lemma 7 of notes: $\text{rank}(X^\top X) = \text{rank}(X)$), so
$$\hat\beta = (X^\top X)^{-1} X^\top Y.$$

### Properties of the LSE

Under (FR) and (SOA):

**(i) Linear in $Y$:** $\hat\beta = AY$ with $A = (X^\top X)^{-1}X^\top$.

**(ii) Unbiased:**
$$\mathbb{E}\hat\beta = (X^\top X)^{-1}X^\top \mathbb{E}Y = (X^\top X)^{-1}X^\top X\beta = \beta.$$

**(iii) Covariance:**
$$\text{Cov}(\hat\beta) = A\cdot\text{Cov}(Y)\cdot A^\top = \sigma^2 (X^\top X)^{-1}X^\top X(X^\top X)^{-1} = \sigma^2 (X^\top X)^{-1}.$$

### The Gauss-Markov theorem

**Statement.** Under (FR) and (SOA), for any $c \in \mathbb{R}^p$, $c^\top \hat\beta$ has the smallest variance among all **linear unbiased** estimators of $c^\top \beta$ (i.e. $\hat\beta$ is **BLUE**).

**Proof** (often asked — 2024 Q3(c), 2021 Q3(e) used its corollary).

$c^\top \hat\beta = c^\top (X^\top X)^{-1}X^\top Y$ is linear in $Y$. Let $\tilde\gamma = L^\top Y$ be any other linear unbiased estimator of $c^\top\beta$. Define
$$D^\top := L^\top - c^\top (X^\top X)^{-1}X^\top, \qquad \text{so } L^\top Y = c^\top\hat\beta + D^\top Y.$$

**Step 1:** $\tilde\gamma$ unbiased $\Rightarrow D^\top X = 0^\top$.
$$c^\top\beta = \mathbb{E}\tilde\gamma = \mathbb{E}(c^\top\hat\beta + D^\top Y) = c^\top\beta + D^\top X\beta \quad \forall\beta.$$
So $D^\top X\beta = 0$ for all $\beta$, forcing $D^\top X = 0^\top$.

**Step 2:** Variance decomposition.
$$\text{Var}(\tilde\gamma) = \text{Var}(c^\top\hat\beta) + \text{Var}(D^\top Y) + 2\text{Cov}(c^\top\hat\beta, D^\top Y).$$
Compute the cross term:
$$\text{Cov}(c^\top\hat\beta, D^\top Y) = c^\top (X^\top X)^{-1}X^\top \cdot \sigma^2 I \cdot D = \sigma^2 c^\top(X^\top X)^{-1}\underbrace{(D^\top X)^\top}_{=0} = 0.$$
Hence $\text{Var}(\tilde\gamma) = \text{Var}(c^\top\hat\beta) + \underbrace{\text{Var}(D^\top Y)}_{\ge 0} \ge \text{Var}(c^\top\hat\beta)$. $\square$

> **Intuition.** Any linear unbiased estimator is $c^\top\hat\beta + (\text{orthogonal noise})$. The noise only inflates variance, never reduces it. So the LSE is optimal in this restricted but natural class.

### The dummy-variable trap / identifiability

When encoding a categorical variable with $k$ levels using $k$ indicator columns **plus** an intercept, the columns are linearly dependent (each row sums to 1 = intercept column), so $X$ does not have full rank and $\beta$ is not identifiable. (This is literally the twins-and-treatments example in the notes.)

**Fix:** drop one indicator column (the "reference level") or drop the intercept.

> This happens in 2022 Q3, 2024 Q2(b), 2025 Q3(c). If you see a model with "male/female dummy" + intercept, alarms should go off. Always be explicit about which category is the baseline.

### Projection matrices

$P \in \mathbb{R}^{n\times n}$ is the **projection matrix onto a subspace $L \subseteq \mathbb{R}^n$** if $Px = x$ for $x \in L$ and $Px = 0$ for $x \in L^\perp$.

**Key characterisation** (asked to prove in 2025 Q3(d)):
$$P \text{ is a projection matrix} \iff P^\top = P \text{ and } P^2 = P.$$

**Proof.**

(⇒) Any $x \in \mathbb{R}^n$ decomposes uniquely as $x = x_L + x_{L^\perp}$.
- *Idempotent:* $P^2 x = P(Px) = P x_L = x_L = Px$, so $P^2 = P$.
- *Symmetric:* For all $x, y$, $x^\top P^\top y = (Px)^\top y = x_L^\top y = x_L^\top y_L$ (since $x_L \perp y_{L^\perp}$), and also $x^\top P y = x^\top y_L = x_L^\top y_L$. So $P^\top = P$.

(⇐) Let $L = \text{span}(P) = \{Pz : z \in \mathbb{R}^n\}$.
- If $x \in L$, say $x = Pz$, then $Px = P^2 z = Pz = x$.
- If $x \in L^\perp$, then for all $y$: $y^\top P x = (P^\top y)^\top x = (Py)^\top x = 0$ (since $Py \in L$ and $x \in L^\perp$). As this holds for all $y$, $Px = 0$. $\square$

**Formula:** If columns of $X$ span $L$ and $X$ has full column rank, then
$$P = X(X^\top X)^{-1}X^\top.$$

**Useful facts about projection matrices** (worth memorising):
- Eigenvalues are 0 or 1.
- $\text{rank}(P) = \text{trace}(P) = \dim L$.
- $I - P$ is the projection onto $L^\perp$.

### Fitted values and residuals

- **Fitted values:** $\hat Y = X\hat\beta = PY$ where $P$ projects onto $\text{span}(X)$.
- **Residuals:** $e = Y - \hat Y = (I - P)Y = QY$, where $Q = I - P$.
- **Residual sum of squares:** $\text{RSS} = e^\top e = Y^\top Q Y$.

**Geometric picture.** $Y$ sits in $\mathbb{R}^n$; the model says $\mathbb{E}Y = X\beta$ lives in the $p$-dimensional subspace $\text{span}(X)$. $\hat Y$ is the nearest point of $\text{span}(X)$ to $Y$ (the orthogonal projection), and $e$ is the part of $Y$ that's orthogonal to the model subspace.

```
          Y
           \
            \  e  (residual, ⊥ span(X))
             \
              ·——————— span(X)
             Ŷ  (fitted values)
```

**Unbiased estimate of $\sigma^2$:**
$$\hat\sigma^2 = \frac{\text{RSS}}{n - r}, \qquad r = \text{rank}(X).$$

The degrees-of-freedom correction $n - r$ is not arbitrary — it's exactly $\text{trace}(Q)$, which is why $\mathbb{E}[\text{RSS}] = \sigma^2 \cdot \text{trace}(Q) = \sigma^2(n-r)$. (For simple iid data with a mean, $r = 1$, giving the familiar $s^2 = \sum(Y_i - \bar Y)^2 / (n-1)$.)

### The coefficient of determination

For models containing an intercept:
$$R^2 = 1 - \frac{\text{RSS}}{\sum_{i=1}^n(Y_i - \bar Y)^2}.$$

Interpretation: fraction of variance in $Y$ "explained" by the fitted model. $R^2 \in [0, 1]$, with $R^2 = 1$ for a perfect fit.

**Big caveat.** Adding *any* column to $X$ will weakly decrease RSS, so $R^2$ can never go down. Therefore $R^2$ on its own is **not** a model-selection criterion. Penalise for model size:

$$R^2_{\text{adj}} = 1 - (1 - R^2)\cdot\frac{n-1}{n-\ell},$$

where $\ell$ is the number of parameters. This can go down when you add a useless predictor. 2022 Q4(b) asked exactly this.

---

## Part VII — Linear models under Normal theory (NTA)

Now assume $\varepsilon \sim N(0, \sigma^2 I_n)$, equivalently $Y \sim N(X\beta, \sigma^2 I_n)$.

### MLE under NTA = LSE

The log-likelihood is
$$\log L(\beta, \sigma^2) = -\frac{n}{2}\log(2\pi\sigma^2) - \frac{1}{2\sigma^2}\underbrace{(Y-X\beta)^\top(Y-X\beta)}_{=S(\beta)}.$$

For fixed $\sigma^2$, maximising over $\beta$ is equivalent to **minimising $S(\beta)$** — same as least squares. So $\hat\beta_{\text{MLE}} = \hat\beta_{\text{LS}}$.

For $\sigma^2$: $\hat\sigma^2_{\text{MLE}} = \text{RSS}/n$ (biased by factor $n/(n-r)$; the unbiased version uses $n-r$).

### The central distributional results

**Lemma 21.** Under (FR) and (NTA): $\text{RSS}/\sigma^2 \sim \chi^2_{n-r}$.

**Proof** (short and beautiful; asked in 2023 Q3(c)):

Let $P = X(X^\top X)^{-1}X^\top$, $Q = I - P$. Then $Q$ is a projection matrix and
$$\frac{\text{RSS}}{\sigma^2} = \left(\frac{Y}{\sigma}\right)^\top Q\left(\frac{Y}{\sigma}\right) = Z^\top Q Z,$$
where $Z = Y/\sigma \sim N(X\beta/\sigma, I)$.

By Lemma 18 of the notes (if $Z \sim N(\mu, I)$ and $A$ is a projection of rank $r$, then $Z^\top A Z \sim \chi^2_r(\delta^2)$ with $\delta^2 = \mu^\top A \mu$), we have $Z^\top Q Z \sim \chi^2_{\text{rank}(Q)}(\delta^2)$ with $\delta^2 = (X\beta/\sigma)^\top Q (X\beta/\sigma)$.

But $QX = (I-P)X = X - X = 0$ (since $P$ projects onto $\text{span}(X)$, $PX = X$), so $\delta^2 = 0$. And $\text{rank}(Q) = n - \text{rank}(P) = n - r$. Therefore $\text{RSS}/\sigma^2 \sim \chi^2_{n-r}$. $\square$

**Lemma 22 (the master $t$-pivot).** Under (FR) and (NTA), for any $c \in \mathbb{R}^p$:
$$\frac{c^\top\hat\beta - c^\top\beta}{\sqrt{c^\top(X^\top X)^{-1}c \cdot \text{RSS}/(n-p)}} \sim t_{n-p}.$$

This is the engine for CIs and tests about *any* linear combination of $\beta$: pick $c = e_i$ for $\beta_i$; pick $c = (1, x_0^\top)^\top$ for $\mathbb{E}(Y | x = x_0)$ in simple linear regression; etc.

**Why it works:** numerator is $N(0, \sigma^2 c^\top(X^\top X)^{-1}c)$, so divided by $\sqrt{\sigma^2 c^\top(X^\top X)^{-1}c}$ it's $N(0,1)$. Denominator is $\sqrt{\text{RSS}/((n-p)\sigma^2)} = \sqrt{\chi^2_{n-p}/(n-p)}$. By independence of $\hat\beta$ and RSS (via Lemma 17: $(X^\top X)^{-1}X^\top Q = 0$), the ratio is by definition $t_{n-p}$.

### The F-test (nested-model testing)

**Setting.** We have a big model $\mathbb{E}Y = X\beta$ with $\text{rank}(X) = r$, and a sub-model $\mathbb{E}Y = X_0 \beta_0$ with $\text{span}(X_0) \subset \text{span}(X)$ and $\text{rank}(X_0) = s$. Test:
$$H_0: \mathbb{E}Y \in \text{span}(X_0) \quad \text{vs} \quad H_1: \mathbb{E}Y \in \text{span}(X) \setminus \text{span}(X_0).$$

**Test statistic:**
$$F = \frac{(\text{RSS}_0 - \text{RSS})/(r - s)}{\text{RSS}/(n-r)} = \frac{(\text{RSS}_0 - \text{RSS})}{\text{RSS}} \cdot \frac{n-r}{r-s}.$$

**Distribution under $H_0$:** $F \sim F_{r-s,\, n-r}$. Reject $H_0$ for large $F$.

**Why this is the right statistic** (Fisher-Cochran). With $Z = Y/\sigma$ and $A_1 = I-P$, $A_2 = P - P_0$, $A_3 = P_0$:
- $A_1 + A_2 + A_3 = I$.
- $A_1, A_3$ are projections; $A_2$ is a projection too because $PP_0 = P_0$ (since columns of $X_0$ are in span of $X$), giving $(P-P_0)^2 = P - P_0$.
- So by Fisher-Cochran: $Z^\top A_i Z$ are independent χ² with d.f. $\text{rank}(A_i) = n-r,\ r-s,\ s$ respectively. Non-centrality is zero under $H_0$ because $P_0(\mathbb{E}Z) = \mathbb{E}Z$ (since $\mathbb{E}Y \in \text{span}(X_0)$).

The numerator of $F$ is $Z^\top A_2 Z/(r-s)$, denominator is $Z^\top A_1 Z/(n-r)$, ratio of independent χ²'s scaled by d.f. = $F$-distribution by definition.

**Common use cases** (very frequent in exam):
- **Test if a single $\beta_i = 0$:** $r - s = 1$. Equivalent to the $t$-test, and $F = t^2$.
- **Test if a group of coefficients are all zero** (e.g. all seasonal dummies): $r - s = $ number of tested coefficients.
- **Test for interaction:** $H_0$ drops the interaction columns only.

> **Intuition.** "How much does RSS *increase* when I'm forced to stay in the smaller model?" A big increase means the extra columns in the bigger model were pulling real weight; a small increase means they added noise. The scaling by $\text{RSS}/(n-r)$ standardises against the noise level.

### Confidence regions for the whole $\beta$ vector

Under (FR), (NTA):
$$\frac{(\hat\beta - \beta)^\top X^\top X (\hat\beta - \beta)}{p \cdot \text{RSS}/(n-p)} \sim F_{p,\,n-p}.$$

So a $(1-\alpha)$ confidence region for $\beta$ is
$$R = \left\{\gamma \in \mathbb{R}^p : (\hat\beta - \gamma)^\top X^\top X (\hat\beta - \gamma) \le p\cdot\frac{\text{RSS}}{n-p}\cdot F_{p,n-p,\alpha}\right\},$$

which is an **ellipsoid** centred at $\hat\beta$. For $p=2$, this is an ellipse in the $(\beta_1, \beta_2)$-plane. Its axes are aligned with the eigenvectors of $X^\top X$, with half-lengths $\sqrt{c/d_i}$ where $d_i$ are the eigenvalues of $X^\top X$ and $c = p\cdot\text{RSS}/(n-p)\cdot F_{p,n-p,\alpha}$.

```
   β₂
   │        ╱ ⟵ confidence ellipse
   │      _╱______
   │    ╱   ·     ╲
   │   │ β̂ =(β̂₁,β̂₂)│
   │    ╲_________╱
   │
   └────────────── β₁
```

**Why you care.** This is **strictly better** than applying Bonferroni to individual CIs: the ellipse fits inside the Bonferroni rectangle, so has smaller area.

---

## Part VIII — Multivariate normal distribution

This is the toolkit that makes the NTA proofs work. It also appears directly on exams (2024 Q3(a), 2025 Q4(a), Problem Sheet 6 Q6 — these are near-identical).

### Three equivalent definitions

$Z \sim N(\mu, \Sigma)$ with $\mu \in \mathbb{R}^n$, $\Sigma$ positive semi-definite, iff any of:
1. Every linear combination $a^\top Z$ is univariate Normal.
2. $Z = \mu + AX$ for some matrix $A$ and $X \sim N(0, I_r)$ iid standard Normal, with $AA^\top = \Sigma$.
3. The characteristic function $\varphi_Z(t) = \exp\!\big(i\mu^\top t - \tfrac12 t^\top \Sigma t\big)$.

### Essential properties

If $Z \sim N(\mu, \Sigma)$ and $A, b$ are deterministic of matching size:
$$AZ + b \sim N(A\mu + b,\, A\Sigma A^\top).$$

**Use this every time you see $Z = AX + b$ with $X$ Normal.**

### Independence via block-diagonality

**Lemma 14.** If $Z \sim N(\mu, \Sigma)$ with blocks $Z = (Z_1^\top, \dots, Z_k^\top)^\top$ and $\Sigma$ is **block-diagonal** (covariances between different blocks are zero), then $Z_1, \dots, Z_k$ are **independent**.

For a general random vector, uncorrelated ≠ independent. But *jointly Normal* + uncorrelated ⇒ independent. This is the single most-tested fact about multivariate Normal in this course.

### Template for the "which components of $Z$ are independent?" exam question

(Appears in 2024 Q3(a), 2025 Q4(a), Problem Sheet 6 Q6.)

**Recipe.**
1. Compute the mean: $\mu_Z = A\mu_X + b$.
2. Compute the covariance: $\Sigma_Z = A \Sigma_X A^\top$.
3. Look at the off-diagonal entries of $\Sigma_Z$. Two components $Z_i, Z_j$ are independent iff $(\Sigma_Z)_{ij} = 0$.

**Worked mini-example.** $X \sim N\!\left(\begin{pmatrix}a\\b\end{pmatrix}, I_2\right)$, $Z = A X + c$ where $A$ is a $3\times 2$ matrix.
$\Sigma_Z = A \cdot I_2 \cdot A^\top = AA^\top$. If $A = \begin{pmatrix}1 & 1\\ 0 & 1\\ 1 & 0\end{pmatrix}$:
$$AA^\top = \begin{pmatrix}2 & 1 & 1\\ 1 & 1 & 0\\ 1 & 0 & 1\end{pmatrix}.$$
So $(\Sigma_Z)_{23} = 0$ ⇒ $Z_2 \perp Z_3$. All other off-diagonals nonzero ⇒ no other independencies.

### Derived distributions

Let $Z_1, \dots, Z_n$ iid $N(0,1)$.

- **$\chi^2_n$:** $\sum Z_i^2$. Non-central version $\chi^2_n(\delta)$: $Z \sim N(\mu, I_n)$, then $Z^\top Z \sim \chi^2_n(\delta)$ with $\delta^2 = \mu^\top\mu$.
- **$t_n$:** $X/\sqrt{U/n}$ where $X \sim N(0,1)$ independent of $U \sim \chi^2_n$. Non-central $t_n(\delta)$: take $X \sim N(\delta, 1)$.
- **$F_{m,n}$:** $(W_1/m)/(W_2/n)$ where $W_1 \sim \chi^2_m, W_2 \sim \chi^2_n$ independent.

### Key independence lemmas (used in the F-test and $t$-statistic)

**Lemma 17.** If $X \sim N(\mu, I)$, $A$ is symmetric PSD, and $BA = 0$, then $X^\top A X$ and $BX$ are independent.

**Lemma 19.** If $X \sim N(\mu, I)$ and $A_1, A_2$ are projection matrices with $A_1 A_2 = 0$, then $X^\top A_1 X$ and $X^\top A_2 X$ are independent.

**Fisher-Cochran Theorem.** If $A_1, \dots, A_k$ are $n \times n$ projection matrices with $\sum A_i = I_n$ and $Z \sim N(\mu, I_n)$, then $Z^\top A_1 Z, \dots, Z^\top A_k Z$ are independent and
$$Z^\top A_i Z \sim \chi^2_{r_i}(\delta_i), \quad r_i = \text{rank}(A_i),\ \delta_i^2 = \mu^\top A_i \mu.$$

This is the single most-used heavy hammer in the NTA chapter.

---

## Part IX — Diagnostics and selection (lighter-touch exam topics)

These don't show up as full questions often but can cameo as short parts.

### Residual anatomy

- **Standardised residual:** $r_i = e_i / \sqrt{\hat\sigma^2 (1 - P_{ii})}$, where $P_{ii}$ is the $i$-th diagonal of the hat matrix. These should look roughly $N(0,1)$-ish if the model is correct.
- **Leverage:** $P_{ii}$ itself. Rule of thumb: "high leverage" if $P_{ii} > 2r/n$. High leverage means the point has disproportionate control over its own fitted value — often an $x_i$ far from the others.
- **Cook's distance:** $D_i = r_i^2 \cdot P_{ii}/[(1-P_{ii})\cdot r]$. Combines large residual AND high leverage into one "how influential is this observation" number. Concern if close to 1.

### Diagnostic plots

- **Residuals vs fitted** / **residuals vs $x_j$:** should look like noise. A funnel shape ⇒ non-constant variance (heteroscedasticity). A curve ⇒ missing non-linear term.
- **QQ plot:** plot sample quantiles vs theoretical $\Phi^{-1}$ quantiles. Should be a straight line if Normality holds.

### Under-/over-fitting

- **Under-fit** (missing a predictor $Z$): LSE of $\beta$ is *biased* unless $X^\top Z = 0$. MSE can sometimes be **smaller** than the full model if $\gamma$ is small — adding bias to cut variance.
- **Over-fit** (including an unnecessary predictor): LSE remains unbiased, but variance of individual coefficients goes *up*.

This is the **bias-variance trade-off** driving model selection:

$$\text{AIC} = -2\log L + 2p, \qquad \text{BIC} = -2\log L + p\log n.$$

Smaller AIC/BIC = better. For linear models, $\text{AIC} = n\log(\text{RSS}/n) + 2p + \text{const}$.

### Weighted least squares (WLS)

If $\text{Cov}(Y) = \sigma^2 V$ (known symmetric PD $V$, not $\sigma^2 I$), then SOA fails. Transform: find $T$ with $T^\top V T = I$; apply OLS to $Z = T^\top Y$ with design $\tilde X = T^\top X$:
$$\hat\beta_{\text{WLS}} = (X^\top V^{-1} X)^{-1} X^\top V^{-1} Y.$$

This is **BLUE** in the transformed model (Gauss-Markov applied to the transformed model). 2022 Q4(a) asked for this.

---

## Part X — Numerical-calculation templates

These appear every year (2023 Q4, 2024 Q4, 2025 Q3/Q4) and are mechanical once you know the moves.

### Template A: compute $\hat\beta$ from provided sample means

**Scenario.** Model $Y_i = \beta_1 + \beta_2 w_i + \varepsilon_i$; you're given $\overline{y}, \overline{w}, \overline{w^2}, \overline{wy}, \overline{y^2}$ (where $\overline{f(y,w)} = \frac{1}{n}\sum f(y_i, w_i)$).

**Moves:**
1. $X^\top X = \begin{pmatrix}n & n\overline{w} \\ n\overline{w} & n\overline{w^2}\end{pmatrix}$, $X^\top Y = \begin{pmatrix}n\overline{y} \\ n\overline{wy}\end{pmatrix}$.
2. $\det(X^\top X) = n^2(\overline{w^2} - \overline{w}^2) = n^2 \cdot \frac{1}{n}\sum(w_i - \overline{w})^2 = n\sum(w_i - \overline{w})^2$.
3. $\hat\beta = (X^\top X)^{-1} X^\top Y$.

The closed form you'll keep rediscovering:
$$\hat\beta_2 = \frac{\overline{wy} - \overline{w}\,\overline{y}}{\overline{w^2} - \overline{w}^2}, \qquad \hat\beta_1 = \overline{y} - \hat\beta_2 \overline{w}.$$

### Template B: compute RSS without computing residuals one-by-one

Use the identity
$$\text{RSS} = Y^\top Y - \hat Y^\top \hat Y = n\overline{y^2} - \hat\beta^\top X^\top Y.$$

Then $\hat\sigma^2 = \text{RSS}/(n-p)$, and SE of $\hat\beta_i$ is $\sqrt{\hat\sigma^2 \cdot [(X^\top X)^{-1}]_{ii}}$.

### Template C: compute $R^2$

$$R^2 = 1 - \frac{\text{RSS}}{n(\overline{y^2} - \overline{y}^2)}.$$

### Template D: one-sided $t$-test for $\beta_i$

$$t_{\text{obs}} = \frac{\hat\beta_i - 0}{\text{SE}(\hat\beta_i)} = \frac{\hat\beta_i}{\sqrt{\hat\sigma^2 \cdot [(X^\top X)^{-1}]_{ii}}}.$$

Reject $H_0: \beta_i \le 0$ at level $\alpha$ iff $t_{\text{obs}} > t_{n-p, \alpha}$. For $\alpha = 0.025$ and large $n-p$, threshold ≈ 1.96. For $\alpha = 0.05$, threshold ≈ 1.64.

### Template E: linear transformation of a Normal

Given $X \sim N(\mu_X, \Sigma_X)$, find distribution and independencies of $Z = AX + b$:
1. $\mu_Z = A\mu_X + b$.
2. $\Sigma_Z = A\Sigma_X A^\top$.
3. Read off zero off-diagonal entries of $\Sigma_Z$ to identify independent pairs.

### Template F: reparametrisation of a linear model

If you rewrite $Y = X\beta + \varepsilon$ in terms of new parameters $\gamma$ with $\beta = M\gamma$ (or $\gamma = M^{-1}\beta$), the new design matrix is $\tilde X = XM$, and $\hat\gamma = M^{-1}\hat\beta$.

**Example** (2023 Q4(d) style): if $v_i = 2 + w_i + z_i$ is a linear combination of existing covariates, the model with $v$ replacing $z$ is obtained by a reparametrisation — algebraically equivalent, not a *new* model. You can either set up the new normal equations or translate between $\hat\beta$ and $\hat\gamma$ directly.

---

## Part XI — Problem-sheet checklist

The exam reuses patterns from the problem sheets heavily. If you're short on time, prioritise these (ranked by how likely the pattern recurs):

### Highest priority — do these first

**Sheet 1 (Weeks 1–2)**
- **Q3** — Prove $\text{MSE} = \text{Var} + \text{bias}^2$. Short, asked verbatim in 2025.
- **Q5** — Prove $\bar Y^2$ is biased for $\mu^2$, and $S$ biased for $\sigma$. Classic trap.
- **Q8** — Compute CRLB for Exponential, Normal, Bernoulli, Poisson. Pure practice.

**Sheet 2 (Week 3)**
- **Q2** — Slutsky giving asymptotic Normality of the $t$-statistic. Core technique.
- **Q3** — Delta method on geometric mean. Good practice for the machinery.
- **Q6** — MLE for Bernoulli, Poisson, Exponential. Absolute bookwork.

**Sheet 3 (Week 4)**
- **Q5** — $2\lambda\sum Y_i \sim \chi^2_{2n}$ for exponential data; building an exact CI. Pivots in action.
- **Q6** — Approximate CI for odds via Delta method. Templates the 2021 Q1 moves.
- **Q7** — Bonferroni CI for $(\mu, \sigma^2)$.

**Sheet 4 (Week 5)**
- **Q3, Q4** — Classical one-sample $t$-test, CIs, and test/CI duality. 2024 Q4 asks this pattern.
- **Q5** — Power function & sample-size calculation. Asked in 2022 Q2(d).
- **Q8** — Two-sample $t$-test **derived** as an LRT. Beautiful, and shows up in 2022 Q2 / 2023 Q1.
- **Q9** — LRT for a Geometric — near-identical to 2023 Q1(c).

**Sheet 5 (Weeks 6–7)**
- **Q1** — Dummy-variable model for two-group comparison. Primes 2022 Q3, 2024 Q4.
- **Q5** — Derivation of simple-linear-regression LSE from scratch. Bedrock.
- **Q7** — Interpretation of $\beta$ under log transforms (geometric-mean story).
- **Q11** — Construct projection matrices onto given subspaces. Concrete geometric practice.

**Sheet 6 (Week 8) — CRITICAL**
- **Q1** — 3-group linear model with dummies + reparametrisations + fitted values. Essentially 2023 Q4 / 2025 Q3.
- **Q2** — Compute $\text{Cov}(e)$ from $Q = I - P$. Foundational.
- **Q6** — **"Which components of $Z$ are independent?"** — this is almost verbatim 2024 Q3(a) and 2025 Q4(a). **Do this one.**
- **Q7** — Conditional distribution of bivariate Normal.

**Sheet 7 (Week 9)**
- **Q1** — Properties of non-central $\chi^2$ (mean, variance, sum). Useful bookwork.
- **Q4** — Equivalence of $t$-test and $F$-test for a single coefficient, $p$-values equal. Conceptually key.
- **Q5** — Interaction terms, testing for main-effects vs interactions. Asked in 2021 Q4 / 2024 Q4.
- **Q6** — Under-/over-fitting bias and variance analysis. Matches the WLS/BLUE discussions from 2022 Q4.

### Medium priority

- Sheet 1 Q10 (censored Bernoulli CRLB), Q6 (linear combinations of estimators)
- Sheet 2 Q8 (MLE of $P(X = 0)$ via functional invariance)
- Sheet 3 Q1 (general estimating-equations asymptotics — extends MLE story)
- Sheet 4 Q6 (chi-squared test for two Binomials)
- Sheet 5 Q9 (error-in-variables — bias of LSE with mismeasured covariates)
- Sheet 6 Q3 (residuals orthogonal to columns + constants)
- Sheet 7 Q2 (normal approximation to $\chi^2_n$ for large $n$)

### Lower priority (skip if short on time)

- Challenging problems on Uniform$(0,\theta)$ (Sheet 1 Q13, Sheet 2 Q4) — interesting, but these distributions rarely appear in exams
- Sheet 2 Q10 (asymptotic relative efficiency in censoring context)
- Sheet 3 Q1, Q2 (general estimating-equations asymptotics, one-step estimator proof)
- R-lab questions — mark scheme doesn't generally assess R code on exam

---

## Part XII — Cheat-sheet of "ask-if-I-forget" facts

- **Fisher information scales with $n$:** $I_n(\theta) = n I_1(\theta)$.
- **MLE asymptotic variance:** $1/(n I_1(\theta))$ — achieves the CRLB.
- **$X^\top X$ invertible** iff $X$ has full column rank.
- **$\hat\beta$ unbiased, $\text{Cov}(\hat\beta) = \sigma^2 (X^\top X)^{-1}$**.
- **$\text{RSS}/(n-r)$ unbiased** for $\sigma^2$; $r = \text{rank}(X)$.
- **Under NTA:** $\hat\beta \sim N(\beta, \sigma^2(X^\top X)^{-1})$, $\text{RSS}/\sigma^2 \sim \chi^2_{n-r}$, and they're **independent**.
- **$P = X(X^\top X)^{-1}X^\top$** is the projection matrix onto $\text{span}(X)$; $\text{trace}(P) = \text{rank}(X) = r$.
- **Projection characterisation:** $P^2 = P$ and $P^\top = P$.
- **$t_n \to N(0,1)$** as $n \to \infty$; $F_{m,n}/m \to \chi^2_m/m$ as $n \to \infty$ for fixed $m$.
- **Adding columns to $X$** never increases RSS, never decreases $R^2$ → use adjusted $R^2$ or AIC.
- **Gauss-Markov** needs only SOA + FR; normality is not required for BLUE.
- **Normality is required** for exact $t$ and $F$ distributions — otherwise they're asymptotic/approximate.
- **Jointly Normal + uncorrelated ⇒ independent.** (Not true for general random variables.)
- **$z_{0.025} \approx 1.96$, $z_{0.05} \approx 1.64$, $z_{0.10} \approx 1.28$.**

---

*Good luck. In the exam: read the question twice, state what assumption you're using (FR/SOA/NTA), and when in doubt, reach for a projection matrix.*