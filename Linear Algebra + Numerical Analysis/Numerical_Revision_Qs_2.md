# MATH50003 Numerical Analysis — Practice Problems (Sheet 2)

A second set of practice problems on the same in-scope topics as Sheet 1, with all-new questions. Topics:

- Calculus on a computer (finite differences, dual numbers)
- Floating-point error analysis
- Cholesky factorisation
- QR factorisation via Householder reflections
- Fourier series and the DFT
- Orthogonal polynomials (Stieltjes algorithm)

> **Tip.** Sheet 1's Q7 covered the *orthogonality of derivatives* technique that appeared as 2025 Q6(b)(ii). This sheet's Q7 covers *Stieltjes*, which appeared as 2025 Q6(b)(i). Together they should cover the orthogonal-polynomial half of an exam Q6.

---

## Problem 1 — Second derivative via central differences

Show that
$$\frac{f(x+h) - 2 f(x) + f(x-h)}{h^2} = f''(x) + \delta$$
where
$$|\delta| \le \frac{M h^2}{12}$$
for $M = \sup_{x-h \le t \le x+h} |f^{(4)}(t)|$.

### Solution

Taylor expand to **fourth** order — we need $f^{(4)}$ because the cubic terms will cancel:
$$f(x+h) = f(x) + f'(x) h + \tfrac{1}{2} f''(x) h^2 + \tfrac{1}{6} f'''(x) h^3 + \tfrac{1}{24} f^{(4)}(t_1) h^4,$$
$$f(x-h) = f(x) - f'(x) h + \tfrac{1}{2} f''(x) h^2 - \tfrac{1}{6} f'''(x) h^3 + \tfrac{1}{24} f^{(4)}(t_2) h^4,$$
for some $t_1 \in [x, x+h]$, $t_2 \in [x-h, x]$.

**Add** the two expressions — the odd-order terms cancel:
$$f(x+h) + f(x-h) = 2 f(x) + f''(x) h^2 + \tfrac{h^4}{24} \left( f^{(4)}(t_1) + f^{(4)}(t_2) \right).$$

Rearrange:
$$f(x+h) - 2 f(x) + f(x-h) = f''(x) h^2 + \tfrac{h^4}{24} \left( f^{(4)}(t_1) + f^{(4)}(t_2) \right).$$

Divide through by $h^2$:
$$\frac{f(x+h) - 2 f(x) + f(x-h)}{h^2} = f''(x) + \tfrac{h^2}{24}\left( f^{(4)}(t_1) + f^{(4)}(t_2) \right).$$

So $\delta = \tfrac{h^2}{24}( f^{(4)}(t_1) + f^{(4)}(t_2) )$ and
$$|\delta| \le \tfrac{h^2}{24} (M + M) = \tfrac{Mh^2}{12}. \quad \blacksquare$$

> **Conceptual note.** Same trick as the centred first-derivative scheme: symmetric stencils kill *odd-order* error terms exactly, so the leading error here jumps from $f'''$ to $f^{(4)}$. This is *the* standard 3-point stencil for $f''$ — you'll see it in every numerical PDE solver.

---

## Problem 2 — Floating-point cumulative product

Let $x_1, \ldots, x_n$ be normal floating-point numbers and define
$$p_n := x_1 \otimes x_2 \otimes \cdots \otimes x_n.$$

Show that
$$p_n = x_1 x_2 \cdots x_n \,(1 + \delta), \quad \text{where } |\delta| \le E_{n-1, \epsilon_m/2},$$
with $E_{k, \epsilon} := \tfrac{k\epsilon}{1 - k\epsilon}$ (assume $n\epsilon_m < 2$).

### Solution

Proceed by **induction on $n$**.

**Base case ($n = 1$):** $p_1 = x_1$ exactly, so $\delta = 0$ and $E_{0, \epsilon_m/2} = 0$. ✓

**Inductive step.** Suppose $p_k = x_1 \cdots x_k (1 + \delta_k)$ with $|\delta_k| \le E_{k-1, \epsilon_m/2}$. The next operation is a single floating-point multiply, which by the IEEE rounding axiom satisfies
$$p_{k+1} = p_k \otimes x_{k+1} = p_k \, x_{k+1} (1 + \eta_k), \qquad |\eta_k| \le \tfrac{\epsilon_m}{2}.$$

Substituting:
$$p_{k+1} = x_1 \cdots x_{k+1} \,(1 + \delta_k)(1 + \eta_k) = x_1 \cdots x_{k+1} \,(1 + \delta_{k+1}),$$
where $\delta_{k+1} := \delta_k + \eta_k + \delta_k \eta_k$.

Bound:
$$|\delta_{k+1}| \le |\delta_k| + |\eta_k| + |\delta_k| |\eta_k| \le E_{k-1,\epsilon_m/2} + \tfrac{\epsilon_m}{2} + E_{k-1,\epsilon_m/2} \cdot \tfrac{\epsilon_m}{2} = \tfrac{\epsilon_m}{2} + E_{k-1, \epsilon_m/2}\,(1 + \tfrac{\epsilon_m}{2}).$$

Now apply the standard algebraic identity (same trick as in the revision sheet's Problem 2(b)):
$$\tfrac{\epsilon_m}{2} + E_{k-1, \epsilon_m/2}(1 + \tfrac{\epsilon_m}{2}) = \frac{\tfrac{\epsilon_m}{2}\,(1 - (k-1)\tfrac{\epsilon_m}{2}) + (k-1)\tfrac{\epsilon_m}{2}\,(1 + \tfrac{\epsilon_m}{2})}{1 - (k-1)\tfrac{\epsilon_m}{2}} = \frac{k \tfrac{\epsilon_m}{2}}{1 - (k-1)\tfrac{\epsilon_m}{2}} \le E_{k, \epsilon_m/2}.$$

Hence $|\delta_{k+1}| \le E_{k, \epsilon_m/2}$, completing the induction. Setting $k+1 = n$ gives the result. $\quad\blacksquare$

> **Conceptual note.** Each multiplication contributes a single $(1 + \eta)$ factor with $|\eta| \le \epsilon_m/2$, and we do $n - 1$ of them — that's why the bound is $E_{n-1, \epsilon_m/2}$, not $E_{n, \epsilon_m/2}$. **Counting the multiplications** is the easy thing to lose marks on. (Compare with the dot product: there's $n$ multiplications and $n-1$ additions, so both contribute.)

---

## Problem 3 — Dual extension of $\log$

What is the dual extension of $f(x) = \log x$? That is, what should $\log(a + b\epsilon)$ equal for $a > 0$?

### Solution

Factor out $a$ and use the Taylor series:
$$\log(a + b\epsilon) = \log\!\Big( a \big( 1 + \tfrac{b}{a} \epsilon \big) \Big) = \log a + \log\!\big(1 + \tfrac{b}{a}\epsilon\big).$$

Now the standard expansion $\log(1+u) = u - u^2/2 + u^3/3 - \cdots$, applied with $u = b\epsilon/a$, collapses thanks to $\epsilon^2 = 0$:
$$\log\!\big(1 + \tfrac{b}{a}\epsilon\big) = \tfrac{b}{a}\epsilon - \tfrac{1}{2}\big(\tfrac{b}{a}\epsilon\big)^2 + \cdots = \tfrac{b}{a}\epsilon.$$

Therefore
$$\boxed{\log(a + b\epsilon) = \log a + \tfrac{b}{a} \epsilon}. \quad \blacksquare$$

(Sanity check: $\tfrac{d}{dx} \log x = \tfrac{1}{x}$, so the rule $f(a + b\epsilon) = f(a) + f'(a) b \epsilon$ gives $\log a + \tfrac{1}{a} b \epsilon$ ✓.)

> **Conceptual note.** Once you know the rule $f(a + b\epsilon) = f(a) + f'(a) b \epsilon$, dual-number questions are mechanical. The exam might ask you to *derive* it for a specific function (which is what's happening here), but you can always sanity-check by recalling the derivative.

---

## Problem 4 — Cholesky factorisation

Use the Cholesky algorithm to determine whether
$$A = \begin{pmatrix} 4 & 4 & 2 \\ 4 & 5 & 3 \\ 2 & 3 & 6 \end{pmatrix}$$
is symmetric positive definite, and if so produce $L$ with $A = L L^\top$.

### Solution

**Step 1.** $\alpha_1 = 4 > 0$, $v = (4, 2)^\top$.
$$A_2 = \begin{pmatrix} 5 & 3 \\ 3 & 6 \end{pmatrix} - \tfrac{1}{4}\begin{pmatrix} 4 \\ 2 \end{pmatrix}\begin{pmatrix} 4 & 2 \end{pmatrix} = \begin{pmatrix} 5 & 3 \\ 3 & 6 \end{pmatrix} - \begin{pmatrix} 4 & 2 \\ 2 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 5 \end{pmatrix}.$$

**Step 2.** $\alpha_2 = 1 > 0$, $v = (1)$.
$$A_3 = (5) - \tfrac{1}{1}(1) = (4).$$

**Step 3.** $\alpha_3 = 4 > 0$.

All pivots positive, so $A$ is SPD with
$$L = \begin{pmatrix} 2 & & \\ 2 & 1 & \\ 1 & 1 & 2 \end{pmatrix}.$$

**Sanity check:**
$$L L^\top = \begin{pmatrix} 2 & 0 & 0 \\ 2 & 1 & 0 \\ 1 & 1 & 2 \end{pmatrix} \begin{pmatrix} 2 & 2 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix} = \begin{pmatrix} 4 & 4 & 2 \\ 4 & 5 & 3 \\ 2 & 3 & 6 \end{pmatrix} = A. \;\checkmark \quad\blacksquare$$

> **Conceptual note.** Each Schur-complement step squashes a $k \times k$ trailing block into a $(k-1) \times (k-1)$ trailing block by subtracting a rank-1 outer product. The diagonal entry of $L$ at step $k$ is $\sqrt{\alpha_k}$, and the column below is $v / \sqrt{\alpha_k}$. If you ever see $\alpha_k \le 0$, stop — $A$ isn't SPD.

---

## Problem 5 — QR factorisation via Householder

Compute the QR factorisation of
$$A = \begin{pmatrix} 8 & 1 \\ 6 & 7 \end{pmatrix}.$$

### Solution

Take $x$ as the first column:
$$x = \begin{pmatrix} 8 \\ 6 \end{pmatrix}, \qquad \|x\| = \sqrt{64 + 36} = 10.$$

Build the Householder vector:
$$y = x - \|x\| e_1 = \begin{pmatrix} -2 \\ 6 \end{pmatrix}, \qquad \|y\| = \sqrt{4 + 36} = 2\sqrt{10}.$$
$$w = \frac{y}{\|y\|} = \frac{1}{\sqrt{10}}\begin{pmatrix} -1 \\ 3 \end{pmatrix}.$$

Then
$$Q = I - 2 w w^\top = I - \tfrac{2}{10}\begin{pmatrix} 1 & -3 \\ -3 & 9 \end{pmatrix} = \begin{pmatrix} 1 - \tfrac{1}{5} & \tfrac{3}{5} \\ \tfrac{3}{5} & 1 - \tfrac{9}{5} \end{pmatrix} = \tfrac{1}{5}\begin{pmatrix} 4 & 3 \\ 3 & -4 \end{pmatrix}.$$

Multiply through:
$$R = QA = \tfrac{1}{5}\begin{pmatrix} 4 & 3 \\ 3 & -4 \end{pmatrix}\begin{pmatrix} 8 & 1 \\ 6 & 7 \end{pmatrix} = \tfrac{1}{5}\begin{pmatrix} 50 & 25 \\ 0 & -25 \end{pmatrix} = \begin{pmatrix} 10 & 5 \\ 0 & -5 \end{pmatrix}.$$

So $A = Q^\top R$ (with $Q$ symmetric, $Q^\top = Q$). $\quad\blacksquare$

> **Conceptual note.** Householder reflections are particularly elegant here: $Q$ is *both* orthogonal and symmetric, so $Q^{-1} = Q^\top = Q$. The reflection sends $x = (8, 6)^\top$ to $(10, 0)^\top$ — exactly $\|x\| e_1$, which is what zeroes out the $(2,1)$ entry of $R$. Note $R$ has a negative diagonal entry; the Householder convention $y = x - \|x\| e_1$ doesn't enforce positive diagonals (unlike the Gram–Schmidt QR convention).

---

## Problem 6 — Fourier coefficients and DFT for $\cos^3 \theta$

For the function $f(\theta) = \cos^3\theta$, give explicit formulae for its Fourier coefficients
$$\hat f_k := \frac{1}{2\pi}\int_0^{2\pi} f(\theta) e^{-ik\theta}\,d\theta$$
and discrete approximations
$$\hat f_k^n := \frac{1}{n}\sum_{j=0}^{n-1} f(\theta_j) e^{-ik\theta_j}, \qquad \theta_j = 2\pi j/n,$$
for all $k \in \mathbb{Z}$, $n = 1, 2, \ldots$.

### Solution

Use the identity $\cos^3\theta = \tfrac{3 \cos\theta + \cos 3\theta}{4}$ (or expand $((e^{i\theta} + e^{-i\theta})/2)^3$):
$$\cos^3\theta = \tfrac{3}{8}e^{i\theta} + \tfrac{3}{8}e^{-i\theta} + \tfrac{1}{8}e^{3i\theta} + \tfrac{1}{8}e^{-3i\theta}.$$

Reading off coefficients:
$$\hat f_{\pm 1} = \tfrac{3}{8}, \quad \hat f_{\pm 3} = \tfrac{1}{8}, \quad \hat f_k = 0 \text{ otherwise}.$$

For the DFT, use $\hat f_k^n = \sum_{j \in \mathbb{Z}} \hat f_{k+nj}$ — work out which residue classes mod $n$ contain $\pm 1$ and $\pm 3$.

**$n = 1$:** all integers are equivalent; $\hat f_k^1 = 2 \cdot \tfrac{3}{8} + 2 \cdot \tfrac{1}{8} = 1$ for all $k$.

**$n = 2$:** $\pm 1, \pm 3$ are all odd, so all collapse onto residue $1$:
$$\hat f_{2k}^2 = 0, \qquad \hat f_{2k+1}^2 = 1.$$

**$n = 3$:** $-3 \equiv 0$, $-1 \equiv 2$, $1 \equiv 1$, $3 \equiv 0 \pmod{3}$:
$$\hat f_{3k}^3 = \hat f_3 + \hat f_{-3} = \tfrac{1}{4}, \qquad \hat f_{3k+1}^3 = \hat f_1 = \tfrac{3}{8}, \qquad \hat f_{3k+2}^3 = \hat f_{-1} = \tfrac{3}{8}.$$

**$n = 4$:** *aliasing collision* — $-3 \equiv 1$ and $3 \equiv -1 \pmod 4$:
$$\hat f_{4k+1}^4 = \hat f_1 + \hat f_{-3} = \tfrac{3}{8} + \tfrac{1}{8} = \tfrac{1}{2}, \quad \hat f_{4k+3}^4 = \hat f_{-1} + \hat f_3 = \tfrac{1}{2}, \quad \hat f_{4k}^4 = \hat f_{4k+2}^4 = 0.$$

**$n = 5$:** $-3 \equiv 2$, $-1 \equiv 4$, $1 \equiv 1$, $3 \equiv 3$ — all distinct:
$$\hat f_{5k}^5 = 0, \;\; \hat f_{5k+1}^5 = \tfrac{3}{8}, \;\; \hat f_{5k+2}^5 = \tfrac{1}{8}, \;\; \hat f_{5k+3}^5 = \tfrac{1}{8}, \;\; \hat f_{5k+4}^5 = \tfrac{3}{8}.$$

**$n = 6$:** another collision — $3 \equiv -3 \pmod 6$:
$$\hat f_{6k+3}^6 = \hat f_3 + \hat f_{-3} = \tfrac{1}{4}, \quad \hat f_{6k+1}^6 = \tfrac{3}{8}, \quad \hat f_{6k+5}^6 = \tfrac{3}{8}, \quad \text{rest } 0.$$

**$n \ge 7$:** the four nonzero modes $\pm 1, \pm 3$ now sit in four distinct residue classes (the smallest pairwise gap, $|3 - (-3)| = 6$, is less than $n$), so:
$$\hat f_{\pm 1 + nk}^n = \tfrac{3}{8}, \qquad \hat f_{\pm 3 + nk}^n = \tfrac{1}{8}, \qquad \text{all other } \hat f_k^n = 0. \quad\blacksquare$$

> **Conceptual note.** This problem has *two* aliasing collisions instead of one, at $n = 4$ (where $1 \equiv -3$) and $n = 6$ (where $3 \equiv -3$). The general rule: a collision happens between modes $k_1$ and $k_2$ at sampling rate $n$ exactly when $n \mid (k_1 - k_2)$. For the support $\{-3, -1, 1, 3\}$ this means $n \in \{2, 4, 6\}$ produce collisions, and $n \ge 7$ resolves everything cleanly.

---

## Problem 7 — Stieltjes algorithm

Use the Stieltjes algorithm to construct the first three monic orthogonal polynomials $p_0, p_1, p_2$ with respect to the inner product
$$\langle f, g \rangle := \int_0^1 f(x) g(x) \,x \,dx.$$

That is, $p_n$ is monic of degree $n$ and $\langle p_n, p_m \rangle = 0$ for $n \ne m$.

The Stieltjes recurrence for monic orthogonal polynomials reads
$$p_0(x) = 1, \quad p_1(x) = x - a_0, \quad p_{n+1}(x) = (x - a_n) p_n(x) - b_n p_{n-1}(x),$$
where
$$a_n = \frac{\langle x p_n, p_n \rangle}{\langle p_n, p_n \rangle}, \quad b_n = \frac{\langle p_n, p_n \rangle}{\langle p_{n-1}, p_{n-1} \rangle}.$$

### Solution

**Step 0.** $p_0(x) = 1$. Compute its norm:
$$\langle p_0, p_0 \rangle = \int_0^1 1 \cdot x \, dx = \tfrac{1}{2}.$$

**Step 1.** Compute $a_0$:
$$\langle x p_0, p_0 \rangle = \int_0^1 x \cdot x \, dx = \tfrac{1}{3}, \quad a_0 = \frac{1/3}{1/2} = \tfrac{2}{3}.$$

So $p_1(x) = x - \tfrac{2}{3}$.

Compute $\langle p_1, p_1 \rangle$:
$$\langle p_1, p_1 \rangle = \int_0^1 \big(x - \tfrac{2}{3}\big)^2 x \, dx = \int_0^1 \big(x^3 - \tfrac{4}{3} x^2 + \tfrac{4}{9} x\big) \, dx = \tfrac{1}{4} - \tfrac{4}{9} + \tfrac{2}{9} = \tfrac{1}{4} - \tfrac{2}{9} = \tfrac{1}{36}.$$

**Step 2.** Compute $a_1$:
$$\langle x p_1, p_1 \rangle = \int_0^1 x \big(x - \tfrac{2}{3}\big)^2 x \,dx = \int_0^1 \big(x^4 - \tfrac{4}{3} x^3 + \tfrac{4}{9} x^2\big)\,dx = \tfrac{1}{5} - \tfrac{1}{3} + \tfrac{4}{27} = \tfrac{2}{135}.$$
$$a_1 = \frac{2/135}{1/36} = \frac{2 \cdot 36}{135} = \tfrac{8}{15}.$$

Compute $b_1$:
$$b_1 = \frac{\langle p_1, p_1 \rangle}{\langle p_0, p_0 \rangle} = \frac{1/36}{1/2} = \tfrac{1}{18}.$$

Therefore
$$p_2(x) = \big(x - \tfrac{8}{15}\big)\big(x - \tfrac{2}{3}\big) - \tfrac{1}{18} = x^2 - \tfrac{6}{5} x + \tfrac{16}{45} - \tfrac{1}{18} = x^2 - \tfrac{6}{5} x + \tfrac{3}{10}.$$

(Algebra: $-\tfrac{8}{15} - \tfrac{2}{3} = -\tfrac{8}{15} - \tfrac{10}{15} = -\tfrac{18}{15} = -\tfrac{6}{5}$ for the $x$ coefficient; $\tfrac{16}{45} - \tfrac{1}{18} = \tfrac{32 - 5}{90} = \tfrac{27}{90} = \tfrac{3}{10}$ for the constant.)

**Summary:**
$$\boxed{p_0 = 1, \qquad p_1 = x - \tfrac{2}{3}, \qquad p_2 = x^2 - \tfrac{6}{5} x + \tfrac{3}{10}.}$$

**Sanity check** that $\langle p_0, p_2 \rangle = 0$:
$$\int_0^1 \big( x^2 - \tfrac{6}{5} x + \tfrac{3}{10} \big) x \, dx = \tfrac{1}{4} - \tfrac{6}{5} \cdot \tfrac{1}{3} + \tfrac{3}{10} \cdot \tfrac{1}{2} = \tfrac{1}{4} - \tfrac{2}{5} + \tfrac{3}{20} = \tfrac{5 - 8 + 3}{20} = 0. \;\checkmark \quad\blacksquare$$

> **Conceptual note.** The Stieltjes algorithm is just **Gram–Schmidt applied to the basis $\{1, x, x^2, \ldots\}$**, but exploiting the structural fact (provable from a three-term recurrence argument) that you only ever need to subtract off the *last two* polynomials, not all of them. That's why the recurrence has just $a_n$ and $b_n$, not coefficients indexed all the way back to $p_0$. This makes Stieltjes $O(n)$ per polynomial rather than $O(n^2)$.
>
> The 2025 exam used Stieltjes for the weight $w(x) = 1 - x^2$ on $[-1, 1]$ (computing $C_0, C_1, C_2$). The procedure is identical — only the integrals change.

---

## Quick-reference revision tips (companion to Sheet 1's tips)

1. **Counting floating-point operations matters.** $n$ multiplications $\Rightarrow$ bound $E_{n-1, \cdot}$ if multiplications-only (cumulative product), $E_{n, \cdot}$ if you also do additions on top (dot product, matrix-vector). One extra operation, one extra index — easy to miscount under exam pressure.

2. **For $f''(x)$ schemes, always Taylor expand to *fourth* order.** The cubic terms cancel by symmetry; you need the quartic terms to bound the error.

3. **Stieltjes intermediate quantities to track:**
   - The norms $\langle p_k, p_k \rangle$ (you need them for both $a_k$ and the next $b_k$).
   - The "shifted" inner products $\langle x p_k, p_k \rangle$ (only need them once each for $a_k$).
   - **Don't skip** sanity-checking that $\langle p_n, p_0 \rangle = 0$ for whichever $p_n$ you computed last — a cheap way to catch arithmetic slips.

4. **For Householder QR, two formulas to memorise:**
   - $y = x - \|x\| e_1$, $w = y / \|y\|$, $Q = I - 2 w w^\top$.
   - $Qx = \|x\| e_1$, $Q^\top = Q$ (symmetric), $Q^2 = I$ (idempotent up to sign).

5. **For DFT, the key sentence to internalise:** "$\hat f_k^n$ is the sum of $\hat f_{k'}$ over all $k' \equiv k \pmod n$." Aliasing collisions occur whenever two modes from the support of $\hat f$ land in the same residue class.

### Suggested further resources

- **Süli & Mayers, *An Introduction to Numerical Analysis*** — Imperial library has plenty of copies. Chapters 1, 6, and 9 are tightly aligned with this module.
- **Trefethen, *Approximation Theory and Approximation Practice*** (free PDF on his website) — gold-standard treatment of Chebyshev/orthogonal polynomial machinery.
- **Stewart, *Afternotes on Numerical Analysis*** — a friendly companion focused on the *intuition* behind floating-point error analysis and orthogonalisation algorithms; very readable.
- **Parlett, *The Symmetric Eigenvalue Problem*** — overkill for this module but the chapter on Stieltjes/Lanczos is the clearest explanation in print of why the three-term recurrence works.