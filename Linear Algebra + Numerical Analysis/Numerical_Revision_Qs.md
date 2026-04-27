# MATH50003 Numerical Analysis — Practice Problems

These problems mirror the style of the official revision sheet and recent exam papers (2024, 2025), staying strictly within the in-scope topics:

- Calculus on a computer (Taylor-error bounds for finite differences, dual numbers)
- Floating-point error analysis (matrix-vector / dot products)
- Cholesky factorisation
- QR factorisation via Householder reflections
- Fourier series and the DFT
- Orthogonal polynomials & weighted inner products

**Excluded** (per your request): II.1 Reals, III.1 Structured Matrices, III.2 LU/PLU, III.5.1 Reduced QR / Gram–Schmidt, IV.2 SVD, V.3 Discrete Convolutions, VI.1.3 Stability of Chebyshev interpolation, and Gauss quadrature.

> **How to use this sheet.** Try each problem fully before reading the solution. The mechanical calculations (Cholesky, QR, DFT) reward practice; the proofs (1, 2, 7) reward thinking about *why* each step works.

---

## Problem 1 — Centred difference error bound

Show that
$$\frac{f(x+h) - f(x-h)}{2h} = f'(x) + \delta$$
where
$$|\delta| \le \frac{Mh^2}{6}$$
for $M = \sup_{x-h \le t \le x+h} |f'''(t)|$.

### Solution

Taylor expand each term to **third** order — we need cubic terms because the linear terms will cancel and $f''$ will too:
$$f(x+h) = f(x) + f'(x)h + \tfrac{1}{2} f''(x) h^2 + \tfrac{1}{6} f'''(t_1) h^3$$
$$f(x-h) = f(x) - f'(x)h + \tfrac{1}{2} f''(x) h^2 - \tfrac{1}{6} f'''(t_2) h^3$$
for some $t_1 \in [x, x+h]$, $t_2 \in [x-h, x]$.

Subtracting kills both the $f(x)$ and $f''(x)$ terms:
$$f(x+h) - f(x-h) = 2 f'(x) h + \tfrac{1}{6}(f'''(t_1) + f'''(t_2)) h^3.$$

Dividing by $2h$:
$$\frac{f(x+h) - f(x-h)}{2h} = f'(x) + \tfrac{1}{12}(f'''(t_1) + f'''(t_2)) h^2.$$

Hence $\delta = \tfrac{1}{12}(f'''(t_1) + f'''(t_2)) h^2$ and
$$|\delta| \le \tfrac{1}{12}(M + M) h^2 = \tfrac{Mh^2}{6}. \quad \blacksquare$$

> **Conceptual note.** This *centred* difference is **second-order accurate** ($O(h^2)$), whereas the asymmetric scheme in the revision sheet's Problem 1 is only first-order ($O(h)$). Symmetric subtraction kills the $f''$ term exactly, so the leading error comes from $f'''$ instead. This is exactly the same idea behind why centred-difference stencils dominate in practical numerical PDE codes.

---

## Problem 2 — Floating-point dot product error

Let $x, y \in \mathbb{F}_{\sigma,Q,S}^n$ be vectors of normal floating-point numbers. Define
$$\mathrm{dotmul}(x, y) := (x_1 \otimes y_1) \oplus (x_2 \otimes y_2) \oplus \cdots \oplus (x_n \otimes y_n).$$

Show that $\mathrm{dotmul}(x, y) = x^\top y + \varepsilon$ where, with $E_{n,\epsilon} := \tfrac{n\epsilon}{1 - n\epsilon}$ and $n\epsilon_m < 2$,
$$|\varepsilon| \le 2 E_{n, \epsilon_m/2}\, \|x\|_\infty \|y\|_1.$$

You may use that $x_1 \oplus \cdots \oplus x_n = x_1 + \cdots + x_n + \sigma_n$ with $|\sigma_n| \le \|x\|_1 E_{n-1,\epsilon_m/2}$.

### Solution

By the IEEE rounding axiom, each individual product satisfies
$$x_j \otimes y_j = x_j y_j (1 + \delta_j), \qquad |\delta_j| \le \tfrac{\epsilon_m}{2}.$$

Now sum, using the given summation bound. Writing $z_j := x_j y_j (1 + \delta_j)$:
$$\mathrm{dotmul}(x, y) = \sum_{j=1}^n z_j + \sigma_n = \sum_{j=1}^n x_j y_j (1+\delta_j) + \sigma_n,$$
where
$$|\sigma_n| \le \|z\|_1 E_{n-1, \epsilon_m/2} = \Big(\sum_j |x_j y_j||1 + \delta_j|\Big) E_{n-1, \epsilon_m/2} \le 2 \|x\|_\infty \|y\|_1 E_{n-1, \epsilon_m/2},$$
using $|1+\delta_j| \le 2$ and the trivial $\sum_j |x_j y_j| \le \|x\|_\infty \|y\|_1$.

Therefore
$$\mathrm{dotmul}(x, y) = x^\top y + \underbrace{\sum_j x_j y_j \delta_j + \sigma_n}_{=\,\varepsilon}.$$

Bounding $\varepsilon$:
$$|\varepsilon| \le \sum_j |x_j y_j| |\delta_j| + |\sigma_n| \le \tfrac{\epsilon_m}{2} \|x\|_\infty \|y\|_1 + 2 E_{n-1, \epsilon_m/2} \|x\|_\infty \|y\|_1.$$

Now the same algebraic trick from the revision sheet's Problem 2(b):
$$\tfrac{\epsilon_m}{2} + 2 E_{n-1, \epsilon_m/2} = \frac{\epsilon_m/2 - (n-1)\epsilon_m^2/4 + 2(n-1)\epsilon_m/2}{1 - (n-1)\epsilon_m/2} \le \frac{2 n \epsilon_m/2}{1 - n \epsilon_m/2} = 2 E_{n, \epsilon_m/2}.$$

Hence $|\varepsilon| \le 2 E_{n, \epsilon_m/2}\, \|x\|_\infty \|y\|_1.$ $\quad\blacksquare$

> **Conceptual note.** The structure is identical to the matrix-vector case: every floating-point operation introduces a relative error $(1+\delta)$ with $|\delta| \le \epsilon_m/2$, and these accumulate **linearly in $n$** (that's what the $E_{n,\epsilon}$ bound captures). The two contributions to $\varepsilon$ are: (i) errors from each multiplication, (ii) error from the summation cascade.

---

## Problem 3 — Dual extension of $\exp$

What is the dual extension of $f(x) = e^x$? That is, what should $\exp(a + b\epsilon)$ equal?

### Solution

Use the Taylor series and the defining property $\epsilon^2 = 0$:
$$\exp(a + b\epsilon) = e^a \cdot e^{b\epsilon} = e^a \left( 1 + b\epsilon + \tfrac{(b\epsilon)^2}{2!} + \cdots \right) = e^a (1 + b\epsilon).$$

So
$$\boxed{\exp(a + b\epsilon) = e^a + e^a b\, \epsilon}. \quad \blacksquare$$

> **Conceptual note.** The pattern is always $f(a + b\epsilon) = f(a) + f'(a) b \epsilon$ — that's the whole point of dual numbers. Differentiation becomes algebra: the second component of $f(a + 1\cdot\epsilon)$ is exactly $f'(a)$. This is the engine behind forward-mode automatic differentiation. For a quick conceptual primer see Wikipedia's [Automatic differentiation page](https://en.wikipedia.org/wiki/Automatic_differentiation) — the "Forward accumulation" section formalises exactly this construction.

---

## Problem 4 — Cholesky factorisation

Use the Cholesky algorithm to determine whether the matrix
$$A = \begin{pmatrix} 4 & 2 & 2 \\ 2 & 5 & 3 \\ 2 & 3 & 6 \end{pmatrix}$$
is symmetric positive definite, and if so produce the lower-triangular $L$ with $A = LL^\top$.

### Solution

**Step 1.** $\alpha_1 = 4 > 0$, $v = (2, 2)^\top$.

$$A_2 = \begin{pmatrix} 5 & 3 \\ 3 & 6 \end{pmatrix} - \tfrac{1}{4} \begin{pmatrix} 2 \\ 2 \end{pmatrix} \begin{pmatrix} 2 & 2 \end{pmatrix} = \begin{pmatrix} 5 - 1 & 3 - 1 \\ 3 - 1 & 6 - 1 \end{pmatrix} = \begin{pmatrix} 4 & 2 \\ 2 & 5 \end{pmatrix}.$$

**Step 2.** $\alpha_2 = 4 > 0$, $v = (2)$.

$$A_3 = (5) - \tfrac{1}{4}(4) = (4).$$

**Step 3.** $\alpha_3 = 4 > 0$.

All pivots are positive, so a Cholesky decomposition exists and $A$ is SPD. The factor is
$$L = \begin{pmatrix} 2 & & \\ 1 & 2 & \\ 1 & 1 & 2 \end{pmatrix}.$$

**Sanity check:** $LL^\top = \begin{pmatrix} 4 & 2 & 2 \\ 2 & 5 & 3 \\ 2 & 3 & 6 \end{pmatrix} = A.$ ✓ $\quad\blacksquare$

> **Conceptual note.** The Cholesky algorithm is just Gaussian elimination *exploiting symmetry* — at each step you peel off a rank-1 outer product $\frac{1}{\alpha_k} v v^\top$ from the trailing block. Because $A$ is SPD, the Schur complement at each step is again SPD (key theorem!), so we never run out of positive pivots. If you ever hit $\alpha_k \le 0$, the matrix isn't SPD — the algorithm doubles as a positive-definiteness test, which is *much* cheaper than computing eigenvalues.

---

## Problem 5 — QR factorisation via Householder

Compute the QR factorisation of
$$A = \begin{pmatrix} 3 & 0 \\ 4 & 5 \end{pmatrix}.$$

### Solution

Take $x$ to be the first column and apply a Householder reflection that maps it onto $e_1$:
$$x = \begin{pmatrix} 3 \\ 4 \end{pmatrix}, \qquad \|x\| = 5.$$

Build the Householder vector:
$$y = x - \|x\| e_1 = \begin{pmatrix} -2 \\ 4 \end{pmatrix}, \qquad \|y\| = \sqrt{4 + 16} = 2\sqrt{5}.$$
$$w = \frac{y}{\|y\|} = \frac{1}{\sqrt{5}} \begin{pmatrix} -1 \\ 2 \end{pmatrix}.$$

Then the reflection
$$Q = I - 2 w w^\top = I - \tfrac{2}{5} \begin{pmatrix} 1 & -2 \\ -2 & 4 \end{pmatrix} = \tfrac{1}{5} \begin{pmatrix} 3 & 4 \\ 4 & -3 \end{pmatrix}.$$

(Note $Q$ is symmetric and orthogonal — both characteristic of Householder reflections.)

Multiply through:
$$R = Q A = \tfrac{1}{5} \begin{pmatrix} 3 & 4 \\ 4 & -3 \end{pmatrix} \begin{pmatrix} 3 & 0 \\ 4 & 5 \end{pmatrix} = \begin{pmatrix} 5 & 4 \\ 0 & -3 \end{pmatrix}.$$

So $A = Q^\top R$ with the matrices above. $\quad\blacksquare$

> **Conceptual note.** The reflection $Q$ is engineered so that $Qx = \|x\| e_1$ — it zeroes everything below the first entry of the first column, which is why $R$ is upper triangular. For larger matrices you iterate: at step $k$ you Householder-reflect within the trailing $(n-k+1) \times (n-k+1)$ block, then the product $Q = Q_1 Q_2 \cdots Q_{n-1}$ is what you store (often implicitly via the $w_k$ vectors). Householder gives **better numerical stability** than Gram–Schmidt because reflections preserve the 2-norm exactly, whereas projecting and subtracting (Gram–Schmidt) suffers from cancellation.

---

## Problem 6 — Fourier coefficients and DFT

For the function $f(\theta) = \cos 2\theta$, give explicit formulae for its Fourier coefficients
$$\hat f_k := \frac{1}{2\pi} \int_0^{2\pi} f(\theta) e^{-ik\theta} \, d\theta$$
and their discrete approximations
$$\hat f_k^n := \frac{1}{n} \sum_{j=0}^{n-1} f(\theta_j) e^{-i k \theta_j}, \qquad \theta_j = 2\pi j/n,$$
for all integers $k$ and all $n = 1, 2, \ldots$.

### Solution

Express $f$ in complex exponentials:
$$\cos 2\theta = \tfrac{1}{2} e^{2i\theta} + \tfrac{1}{2} e^{-2i\theta}.$$

Reading off coefficients: $\hat f_2 = \hat f_{-2} = \tfrac{1}{2}$, and $\hat f_k = 0$ otherwise.

For the DFT we use the **aliasing formula** $\hat f_k^n = \sum_{j \in \mathbb{Z}} \hat f_{k + nj}$. Only $\hat f_{\pm 2}$ are nonzero, so we just need to ask: for each residue class $k \pmod n$, do $-2$ or $+2$ fall in it?

$n = 1$: every integer is $\equiv 0 \pmod 1$, so $\hat f_k^1 = \hat f_{-2} + \hat f_2 = 1$ for all $k$.

$n = 2$:
$$\hat f_{2k}^2 = \hat f_{-2} + \hat f_2 = 1, \qquad \hat f_{2k+1}^2 = 0.$$

$n = 3$:
$$\hat f_{3k}^3 = 0, \qquad \hat f_{3k+1}^3 = \hat f_{-2} = \tfrac{1}{2}, \qquad \hat f_{3k+2}^3 = \hat f_{2} = \tfrac{1}{2}$$
(since $-2 \equiv 1 \pmod 3$ and $2 \equiv 2 \pmod 3$).

$n = 4$: the critical case — $-2 \equiv 2 \pmod 4$, so the modes alias on top of each other:
$$\hat f_{4k+2}^4 = \hat f_{-2} + \hat f_2 = 1, \qquad \hat f_{4k}^4 = \hat f_{4k+1}^4 = \hat f_{4k+3}^4 = 0.$$

For $n \ge 5$: the modes $k = 2$ and $k = -2$ live in *different* residue classes mod $n$ (since $|-2 - 2| = 4 < n$), so they no longer alias:
$$\hat f_{2 + nk}^n = \hat f_2 = \tfrac{1}{2}, \quad \hat f_{-2 + nk}^n = \hat f_{-2} = \tfrac{1}{2}, \quad \text{all other } \hat f_k^n = 0. \quad\blacksquare$$

> **Conceptual note.** The DFT *collapses* Fourier modes that differ by multiples of $n$ — that's **aliasing**. Once $n$ is large enough to separate the support of $\hat f$, the DFT recovers $\hat f_k$ exactly on each residue class. The general rule: if $\hat f$ is supported on $\{-K, \ldots, K\}$, you need $n \ge 2K+1$ to fully resolve $f$ (this is the discrete Nyquist criterion). Here $K = 2$, so $n \ge 5$ is the threshold — exactly where the answers stabilise.

---

## Problem 7 — Orthogonal polynomials and derivatives

The Legendre polynomials $\{P_n\}_{n \ge 0}$ are characterised (up to scaling) as the polynomials with $\deg P_n = n$ and
$$\int_{-1}^{1} P_n(x) P_m(x) \, dx = 0 \qquad (n \ne m).$$

### Part (a)

Prove that for any polynomial $p$ of degree $m < n$,
$$\int_{-1}^{1} P_n(x) p(x) \, dx = 0.$$

#### Solution

Since $\deg P_k = k$ for each $k$, the polynomials $P_0, P_1, \ldots, P_m$ are linearly independent and span the space of polynomials of degree $\le m$. So we can expand
$$p(x) = \sum_{k=0}^{m} c_k P_k(x)$$
for some coefficients $c_k$. Then
$$\int_{-1}^{1} P_n(x) p(x) \, dx = \sum_{k=0}^{m} c_k \underbrace{\int_{-1}^{1} P_n(x) P_k(x) \, dx}_{= \,0 \text{ since } k \ne n} = 0. \quad\blacksquare$$

### Part (b)

Use part (a) to show that the polynomials
$$Q_n(x) := P_{n+1}'(x)$$
are orthogonal with respect to the inner product
$$\langle f, g \rangle = \int_{-1}^{1} f(x) g(x) (1 - x^2) \, dx.$$

#### Solution

Take $m < n$ and integrate by parts, exploiting $(1 - x^2)\big|_{x = \pm 1} = 0$:
$$\int_{-1}^{1} Q_n(x) Q_m(x) (1-x^2) \, dx = \int_{-1}^{1} P_{n+1}'(x) \, P_{m+1}'(x) (1-x^2) \, dx$$
$$= \underbrace{\Big[ P_{n+1}(x)\, P_{m+1}'(x)(1-x^2) \Big]_{-1}^{1}}_{= \,0} - \int_{-1}^{1} P_{n+1}(x) \, \frac{d}{dx}\Big[ P_{m+1}'(x)(1-x^2) \Big] dx.$$

Compute the inner derivative:
$$\frac{d}{dx}\Big[ P_{m+1}'(x)(1-x^2) \Big] = \underbrace{P_{m+1}''(x)(1-x^2)}_{\deg \le m+1} \;+\; \underbrace{(-2x) P_{m+1}'(x)}_{\deg \le m+1}.$$

So this is a polynomial of degree at most $m + 1$.

Since $m < n$, we have $m + 1 \le n < n + 1 = \deg P_{n+1}$. By part (a), $P_{n+1}$ is orthogonal to every polynomial of degree less than $n+1$, hence
$$\int_{-1}^{1} P_{n+1}(x) \cdot (\text{polynomial of degree} \le m+1) \, dx = 0,$$
giving $\langle Q_n, Q_m \rangle = 0$. The case $m > n$ follows by symmetry (swap the integration by parts). $\quad\blacksquare$

> **Conceptual note.** This is exactly the technique used in the revision sheet's Problem 7 *and* it appeared as Q6(b)(ii) on the May 2025 paper — so it's worth knowing cold. The structural pattern is always:
>
> 1. **Integration by parts**, with the weight chosen so that the boundary term vanishes (here, $(1-x^2)$ kills the boundary at $\pm 1$).
> 2. **Compute the resulting derivative** and check its degree.
> 3. **Invoke part (a)** to conclude orthogonality of the original family kills the integral.
>
> The boundary-killing weight is always the punchline. For Chebyshev second-kind $U_n$, the analogous derived family $C_n = U_{n+1}'$ is orthogonal w.r.t. $(1-x^2)^{3/2}$, and the proof has the *same shape* — see the revision sheet's Problem 7(b).

---

## Quick-reference revision tips

1. **Floating-point bounds (Q1, Q2 territory).** The algebra of $E_{n,\epsilon} = \frac{n\epsilon}{1-n\epsilon}$ is fiddly but follows a fixed recipe: each operation introduces relative error $\le \epsilon_m/2$, errors accumulate at most $n$ times, and the master inequality
$$\tfrac{\epsilon_m}{2} + 2 E_{n-1, \epsilon_m/2} \le 2 E_{n, \epsilon_m/2}$$
glues together the per-multiplication and accumulated-summation errors.

2. **Cholesky and QR.** Always sanity-check by multiplying out: if $LL^\top \ne A$ or $Q^\top Q \ne I$, you've made an arithmetic slip. Cholesky pivots being positive ⇔ matrix is SPD.

3. **DFT problems.** Express $f$ in complex exponentials first, read off $\hat f_k$, then apply the aliasing formula $\hat f_k^n = \sum_j \hat f_{k+nj}$. Watch carefully for residue collisions like $-2 \equiv 2 \pmod 4$.

4. **Orthogonal-polynomial proofs.** The pattern is *orthogonality $\Rightarrow$ basis $\Rightarrow$ polynomial expansion $\Rightarrow$ orthogonality of derived family via integration by parts*. The boundary term vanishing is the linchpin — it's the reason the weight always picks up an extra factor like $(1-x^2)$ when you differentiate.

5. **Dual numbers.** Memorise the rule $f(a + b\epsilon) = f(a) + f'(a) b\, \epsilon$ — it's literally just a directional derivative dressed up algebraically. Great for catching out cubics, square roots, and exponentials in a couple of lines.

### Suggested further resources

- **Trefethen & Bau, *Numerical Linear Algebra*** — Chapters 10 (Householder), 16 (stability), 23 (Cholesky) are the gold standard for the linear-algebra half. Available cheaply secondhand and in Imperial library.
- **Heath, *Scientific Computing*** — clearer treatment of error analysis if you find Higham too dense.
- **3Blue1Brown's videos on the DFT and Fourier series** — for visual intuition before grinding through proofs.
- **Boyd, *Chebyshev and Fourier Spectral Methods*** (free online PDF) — comprehensive on Chebyshev / orthogonal polynomial machinery, with a friendly tone.