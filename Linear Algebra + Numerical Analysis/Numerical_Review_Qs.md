# MATH50003 Numerical Analysis — Practice Questions

*Styled after the 2025–26 revision sheet and past exam papers. Covers the topics you've been working on: floating point rounding, divided differences, dual numbers, backward error analysis, Cholesky/LU/QR, Fourier coefficients, Chebyshev polynomials ($T_n$ and $U_n$), orthogonal polynomial proofs, and Gauss quadrature.*

---

## Question 1 — Divided Differences (revision sheet style)

Show that

$$\frac{f(x + 3h) - f(x - 2h)}{5h} = f'(x) + \delta$$

where

$$|\delta| \leq \frac{13Mh}{10}$$

for $M = \sup_{x-2h \leq t \leq x+3h} |f''(t)|$.

**Hint:** Expand $f(x+3h)$ and $f(x-2h)$ separately using Taylor's theorem with remainder, then subtract.

---

## Question 2 — Backward Error Analysis (revision sheet style)

**(a)** Consider an upper bidiagonal matrix with floating point entries:

$$U = \begin{bmatrix} u_{11} & u_{12} \\ & u_{22} & u_{23} \\ & & u_{33} \end{bmatrix} \in F_{\sigma,Q,S}^{3 \times 3}$$

and a vector $\mathbf{x} \in F_{\sigma,Q,S}^3$. Denoting matrix-vector multiplication implemented using floating point arithmetic as

$$\mathbf{b} := \texttt{upperbidiagmul}(U, \mathbf{x})$$

express the entries $b_k := \mathbf{e}_k^\top \mathbf{b}$ in terms of $u_{kj}$ and $x_j := \mathbf{e}_j^\top \mathbf{x}$, using rounded floating point operations $\oplus$ and $\otimes$.

**(b)** Assuming all operations involve normal floating point numbers, show that your approximation has the form

$$U\mathbf{x} = \texttt{upperbidiagmul}(U, \mathbf{x}) + \boldsymbol{\epsilon}$$

where

$$\|\boldsymbol{\epsilon}\|_1 \leq \frac{3}{2}\epsilon_m \|U\|_1 \|\mathbf{x}\|_1$$

using $\|A\|_1 := \max_j \sum_{k=1}^n |a_{kj}|$ and $\|\mathbf{x}\|_1 := \sum_k |x_k|$.

You may assume $\epsilon_m < 1$.

---

## Question 3 — Dual Numbers

**(a)** What is the dual extension of $1/x$? That is, what should $(a + b\epsilon)^{-1}$ equal for $a \neq 0$?

*Derive this by requiring $(a+b\epsilon)(a+b\epsilon)^{-1} = 1$.*

**(b)** Consider a "triple" dual number $x + a\epsilon_1 + b\epsilon_2 + c\epsilon_3$ where

$$\epsilon_1^2 = \epsilon_2, \quad \epsilon_2^2 = 0, \quad \epsilon_1 \epsilon_2 = \epsilon_3, \quad \epsilon_3^2 = \epsilon_2 \epsilon_3 = \epsilon_1 \epsilon_3 = 0.$$

(i) Derive the formula for the product of two triple dual numbers

$$(x + a\epsilon_1 + b\epsilon_2 + c\epsilon_3)(y + d\epsilon_1 + e\epsilon_2 + f\epsilon_3).$$

(ii) Show for all polynomials $p(x)$ that

$$p(x + a\epsilon_1 + b\epsilon_2 + c\epsilon_3) = p(x) + ap'(x)\epsilon_1 + \left(bp'(x) + \frac{a^2}{2}p''(x)\right)\epsilon_2 + \left(cp'(x) + abp''(x) + \frac{a^3}{6}p'''(x)\right)\epsilon_3.$$

(iii) Explain how $a$, $b$, and $c$ can be chosen to compute $p'''(x)$.

---

## Question 4 — Cholesky Factorisation

Use the Cholesky factorisation to determine whether the following matrix is symmetric positive definite:

$$A = \begin{bmatrix} 4 & 2 & -2 \\ 2 & 5 & 1 \\ -2 & 1 & 6 \end{bmatrix}$$

If it is, give $L$ such that $A = LL^\top$.

---

## Question 5 — QR Factorisation by Householder

Compute the QR factorisation of

$$A = \begin{bmatrix} 1 & 2 \\ 2 & 1 \end{bmatrix}$$

using a Householder reflection. That is, find an orthogonal $Q$ and upper triangular $R$ with $A = QR$.

---

## Question 6 — Fourier Coefficients and the DFT

For the function $f(\theta) = \cos^2 \theta$, state explicit formulae for its Fourier coefficients

$$\hat{f}_k := \frac{1}{2\pi} \int_0^{2\pi} f(\theta) e^{-ik\theta}\, d\theta$$

and their discrete approximation:

$$\hat{f}_k^n := \frac{1}{n} \sum_{j=0}^{n-1} f(\theta_j) e^{-ik\theta_j}$$

for all integers $k$, $n = 1, 2, \ldots$, where $\theta_j = 2\pi j / n$.

---

## Question 7 — Chebyshev Polynomials of the First Kind

**(a)** The Chebyshev polynomials of the first kind are defined by $T_n(x) = \cos(n \arccos x)$. Using the three-term recurrence $T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)$, compute $T_0(x)$ through $T_4(x)$.

**(b)** The roots of $T_n(x)$ are $x_j = \cos\left(\frac{\pi(j-1/2)}{n}\right)$ for $j = 1, \ldots, n$. The Chebyshev Vandermonde matrix is

$$V_{jk} = T_{k-1}(x_j), \qquad j,k = 1,\ldots,n.$$

It is known that $V^{-1} = \frac{1}{n}\begin{pmatrix} 1 \\ & 2 \\ & & \ddots \\ & & & 2 \end{pmatrix} V^\top$.

Show that $\|V^{-1}\| \leq \sqrt{2/n}$, and explain why this means Chebyshev interpolation is *numerically stable* — the condition number $\kappa(V)$ is bounded independent of $n$.

**Hint:** The matrix $DV^\top$ (with appropriate $D$) is a scaled version of an orthogonal matrix (the DCT).

---

## Question 8 — Chebyshev Polynomials of the Second Kind

**(a)** The Chebyshev polynomials of the second kind are $U_n(x) = 2^n x^n + O(x^{n-1})$, orthogonal with respect to $\sqrt{1-x^2}$ on $[-1,1]$. Show that

$$U_n(\cos\theta) = \frac{\sin(n+1)\theta}{\sin\theta}$$

and that they satisfy the three-term recurrence

$$xU_0(x) = \frac{U_1(x)}{2}, \qquad xU_n(x) = \frac{U_{n-1}(x) + U_{n+1}(x)}{2}, \quad n > 0.$$

You may use the trigonometric identity $2\cos\alpha \sin\beta = \sin(\alpha+\beta) - \sin(\alpha-\beta)$.

**(b)** Show that $T_{n+1}'(x) = (n+1)U_n(x)$.

**Hint:** Differentiate $T_{n+1}(\cos\theta) = \cos(n+1)\theta$ implicitly with respect to $x$, using $x = \cos\theta$.

**(c)** Using part (b), or otherwise, prove that the polynomials $C_n(x) := U_{n+1}'(x)$ are orthogonal with respect to the inner product

$$\langle f, g \rangle = \int_{-1}^{1} f(x)g(x)(1-x^2)^{3/2}\,dx.$$

**Hint:** Use integration by parts. You may use the fact (from part (a) or the revision sheet Q7a) that $\int_{-1}^1 U_n(x) p(x) \sqrt{1-x^2}\,dx = 0$ for any polynomial $p$ of degree $m < n$.

---

## Question 9 — Orthogonal Polynomial Construction (Stieltjes)

Consider orthogonal polynomials $q_n(x)$ with respect to the weight $w(x) = 1 - x^2$ on $[-1,1]$, normalised so that $q_n(x) = k_n x^n + O(x^{n-1})$ with

$$k_n = \frac{(2(n+1))!}{2^{n+1}(n+1)!n!}.$$

**(a)** Show that $w(x)$ is an even function and hence that the recurrence coefficients $a_n = 0$ for all $n$. That is, the three-term recurrence takes the form

$$xq_n(x) = b_{n-1}q_{n-1}(x) + b_n q_{n+1}(x).$$

**(b)** Construct $q_0(x)$ and $q_1(x)$, and hence find $q_2(x)$. You may use without proof the formulae

$$\int_{-1}^1 (1-x^2)\,dx = \frac{4}{3}, \qquad \int_{-1}^1 x^2(1-x^2)\,dx = \frac{4}{15}.$$

**(c)** Write down the $2 \times 2$ Jacobi matrix $J_2$ for these polynomials. Find its eigenvalues, and hence give the nodes for the 2-point Gauss quadrature rule for $w(x) = 1 - x^2$. Find the corresponding weights.

---

## Question 10 — Floating Point Rounding

Consider IEEE 16-bit floating point numbers $F_{16} = F_{15,5,10}$.

**(a)** What are the following numbers when rounded to the nearest float in $F_{16}$?

$$2, \qquad 2 + 2^{-10}, \qquad 2 + 2^{-9} + 2^{-10}$$

**(b)** What is the spacing between consecutive floats in the interval $[2, 4)$? How does this compare to the spacing in $[1, 2)$?

**(c)** Show that for all $x$ in the normalised range,

$$\text{fl}(x) = x(1 + \delta_x), \qquad |\delta_x| \leq \frac{\epsilon_m}{2}.$$

---

## Solutions

### Solution 1

By Taylor's theorem with remainder:

$$f(x + 3h) = f(x) + 3hf'(x) + \frac{(3h)^2}{2}f''(t_1) = f(x) + 3hf'(x) + \frac{9h^2}{2}f''(t_1)$$

$$f(x - 2h) = f(x) - 2hf'(x) + \frac{(2h)^2}{2}f''(t_2) = f(x) - 2hf'(x) + 2h^2 f''(t_2)$$

for some $t_1 \in [x, x+3h]$ and $t_2 \in [x-2h, x]$. Subtracting:

$$f(x+3h) - f(x-2h) = 5hf'(x) + \frac{9h^2}{2}f''(t_1) - 2h^2 f''(t_2)$$

Dividing by $5h$:

$$\frac{f(x+3h) - f(x-2h)}{5h} = f'(x) + \frac{9h}{10}f''(t_1) - \frac{2h}{5}f''(t_2)$$

The error is:

$$|\delta| = \left|\frac{9h}{10}f''(t_1) - \frac{2h}{5}f''(t_2)\right| \leq \frac{9h}{10}M + \frac{2h}{5}M = \frac{9 + 4}{10}Mh = \frac{13Mh}{10}. \qquad \square$$

---

### Solution 2

**(a)** For the upper bidiagonal matrix, the multiplication proceeds row by row. For $k = 1$:

$$b_1 = (u_{11} \otimes x_1) \oplus (u_{12} \otimes x_2)$$

For $k = 2$:

$$b_2 = (u_{22} \otimes x_2) \oplus (u_{23} \otimes x_3)$$

For $k = 3$: $b_3 = u_{33} \otimes x_3$ (only one term, no addition needed).

More generally for an $n \times n$ upper bidiagonal: $b_k = (u_{kk} \otimes x_k) \oplus (u_{k,k+1} \otimes x_{k+1})$ for $k < n$, and $b_n = u_{nn} \otimes x_n$.

**(b)** Expanding for a general row $k < n$:

$$b_k = (u_{kk}x_k(1+\delta_1) + u_{k,k+1}x_{k+1}(1+\delta_2))(1+\delta_3)$$

where $|\delta_i| \leq \epsilon_m/2$. This gives:

$$b_k = u_{kk}x_k(1+\delta_1)(1+\delta_3) + u_{k,k+1}x_{k+1}(1+\delta_2)(1+\delta_3)$$

The error for this row is:

$$\epsilon_k = b_k - (u_{kk}x_k + u_{k,k+1}x_{k+1})$$
$$= u_{kk}x_k(\delta_1 + \delta_3 + \delta_1\delta_3) + u_{k,k+1}x_{k+1}(\delta_2 + \delta_3 + \delta_2\delta_3)$$

Using $|\delta_i + \delta_j + \delta_i\delta_j| \leq \epsilon_m/2 + \epsilon_m/2 + \epsilon_m^2/4 \leq \epsilon_m + \epsilon_m/4 \leq \frac{3}{2}\epsilon_m$ (since $\epsilon_m < 1$):

$$|\epsilon_k| \leq \frac{3}{2}\epsilon_m(|u_{kk}||x_k| + |u_{k,k+1}||x_{k+1}|)$$

For $k = n$: $b_n = u_{nn}x_n(1+\delta_1)$, so $|\epsilon_n| \leq \frac{\epsilon_m}{2}|u_{nn}||x_n| \leq \frac{3}{2}\epsilon_m|u_{nn}||x_n|$.

Summing:

$$\|\boldsymbol{\epsilon}\|_1 = \sum_k |\epsilon_k| \leq \frac{3}{2}\epsilon_m \sum_k \sum_j |u_{kj}||x_j| \leq \frac{3}{2}\epsilon_m \|U\|_1 \|\mathbf{x}\|_1. \qquad \square$$

---

### Solution 3

**(a)** Require $(a + b\epsilon)(c + d\epsilon) = 1$, i.e., $ac + (bc + ad)\epsilon = 1$. So $ac = 1$ giving $c = 1/a$, and $bc + ad = 0$ giving $d = -b/a^2$. Thus:

$$(a + b\epsilon)^{-1} = \frac{1}{a} - \frac{b}{a^2}\epsilon.$$

This is consistent with $(1/x)' = -1/x^2$.

**(b)(i)** Expanding and using the nilpotent rules:

$(x + a\epsilon_1 + b\epsilon_2 + c\epsilon_3)(y + d\epsilon_1 + e\epsilon_2 + f\epsilon_3)$

$= xy + (xd + ay)\epsilon_1 + (xe + by + ad)\epsilon_2 + (xf + cy + ae + bd)\epsilon_3$

where we used $\epsilon_1^2 = \epsilon_2$, $\epsilon_1\epsilon_2 = \epsilon_3$, and all other quadratic/higher terms vanish.

**(b)(ii)** By induction. For $p(x) = x^n$, the base cases $n = 0, 1$ are clear. For the inductive step, multiply $(x + a\epsilon_1 + b\epsilon_2 + c\epsilon_3)$ by $x^{n-1} + (n-1)a x^{n-2}\epsilon_1 + ((n-1)bx^{n-2} + \frac{(n-1)(n-2)}{2}a^2 x^{n-3})\epsilon_2 + (\cdots)\epsilon_3$ using the product formula from (i). The $\epsilon_3$ coefficient collects terms giving $cnx^{n-1} + abn(n-1)x^{n-2} + \frac{n(n-1)(n-2)}{6}a^3 x^{n-3}$, which is $cp'(x) + abp''(x) + \frac{a^3}{6}p'''(x)$ evaluated at $x^n$. Linearity extends to general polynomials.

**(b)(iii)** Choose $a = 1$, $b = 0$, $c = 0$. Then the $\epsilon_3$ coefficient gives $\frac{1}{6}p'''(x)$, so $p'''(x) = 6 \times [\epsilon_3\text{-coefficient}]$.

---

### Solution 4

Step 1: $\alpha_1 = A[1,1] = 4$ and $\mathbf{v} = \frac{1}{4}[2, -2]^\top = [1/2, -1/2]^\top$.

$$A_2 = \begin{bmatrix} 5 & 1 \\ 1 & 6 \end{bmatrix} - \frac{1}{4}\begin{bmatrix} 2 \\ -2 \end{bmatrix}\begin{bmatrix} 2 & -2 \end{bmatrix} = \begin{bmatrix} 5 & 1 \\ 1 & 6 \end{bmatrix} - \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix} = \begin{bmatrix} 4 & 2 \\ 2 & 5 \end{bmatrix}$$

Step 2: $\alpha_2 = 4$ and $\mathbf{v}_2 = [2/4] = [1/2]$.

$$A_3 = [5] - \frac{1}{4}[2][2] = 5 - 1 = [4]$$

Since $\alpha_3 = 4 > 0$ (and all $\alpha_i > 0$), $A$ is SPD. The Cholesky factor is:

$$L = \begin{bmatrix} 2 & & \\ 1 & 2 & \\ -1 & 1 & 2 \end{bmatrix}$$

since $\ell_{11} = \sqrt{4} = 2$, first column below: $[2/2, -2/2]^\top = [1, -1]^\top$; $\ell_{22} = \sqrt{4} = 2$, second column below: $[2/2] = [1]$; $\ell_{33} = \sqrt{4} = 2$.

Verify: $LL^\top = A$. $\checkmark$

---

### Solution 5

First column of $A$: $\mathbf{x} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$, $\|\mathbf{x}\| = \sqrt{5}$.

$\mathbf{y} = \mathbf{x} - \|\mathbf{x}\|\mathbf{e}_1 = \begin{pmatrix} 1 - \sqrt{5} \\ 2 \end{pmatrix}$

$\|\mathbf{y}\|^2 = (1-\sqrt{5})^2 + 4 = 1 - 2\sqrt{5} + 5 + 4 = 10 - 2\sqrt{5}$

$\mathbf{w} = \frac{\mathbf{y}}{\|\mathbf{y}\|}$

$Q = I - 2\mathbf{w}\mathbf{w}^\top = I - \frac{2\mathbf{y}\mathbf{y}^\top}{\|\mathbf{y}\|^2} = I - \frac{2}{10-2\sqrt{5}}\begin{pmatrix} (1-\sqrt{5})^2 & 2(1-\sqrt{5}) \\ 2(1-\sqrt{5}) & 4 \end{pmatrix}$

$= I - \frac{1}{5-\sqrt{5}}\begin{pmatrix} 6-2\sqrt{5} & 2-2\sqrt{5} \\ 2-2\sqrt{5} & 4 \end{pmatrix}$

Simplifying (multiply numerator and denominator by $5+\sqrt{5}$):

$$Q = \frac{1}{\sqrt{5}}\begin{bmatrix} -1 & -2 \\ -2 & 1 \end{bmatrix}$$

(Note the sign choice: $Q\mathbf{x} = -\|\mathbf{x}\|\mathbf{e}_1$ with the default sign convention.)

Then $R = QA$:

$$R = \frac{1}{\sqrt{5}}\begin{bmatrix} -1 & -2 \\ -2 & 1 \end{bmatrix}\begin{bmatrix} 1 & 2 \\ 2 & 1 \end{bmatrix} = \frac{1}{\sqrt{5}}\begin{bmatrix} -5 & -4 \\ 0 & -3 \end{bmatrix} = \begin{bmatrix} -\sqrt{5} & -4/\sqrt{5} \\ 0 & -3/\sqrt{5} \end{bmatrix}$$

Since $Q$ is its own inverse ($Q^2 = I$ for reflections), $A = QR$.

If positive diagonal is required: flip signs by noting $A = (Q)(-R) \cdot (-I)$... or simply note $A = Q^\top R$ gives $Q^\top = Q$ since $Q$ is symmetric. The factorisation is $A = QR$ as stated.

---

### Solution 6

Write $\cos^2\theta = \frac{1 + \cos 2\theta}{2} = \frac{1}{2} + \frac{e^{2i\theta} + e^{-2i\theta}}{4}$.

So the exact Fourier coefficients are:

$$\hat{f}_0 = \frac{1}{2}, \qquad \hat{f}_2 = \hat{f}_{-2} = \frac{1}{4}, \qquad \hat{f}_k = 0 \text{ otherwise.}$$

For the discrete coefficients, we use $\hat{f}_k^n = \sum_{m=-\infty}^{\infty} \hat{f}_{k+mn}$ (aliasing formula):

- $n = 1$: $\hat{f}_k^1 = \hat{f}_0 + \hat{f}_2 + \hat{f}_{-2} = 1$ for all $k$.
- $n = 2$: $\hat{f}_{2k}^2 = \hat{f}_0 + \hat{f}_2 + \hat{f}_{-2} = 1$, $\hat{f}_{2k+1}^2 = 0$.
- $n = 3$: $\hat{f}_{3k}^3 = \hat{f}_0 = 1/2$, $\hat{f}_{3k+1}^3 = \hat{f}_{-2} = 1/4$, $\hat{f}_{3k+2}^3 = \hat{f}_2 = 1/4$.
- $n = 4$: $\hat{f}_{4k}^4 = 1/2$, $\hat{f}_{4k+2}^4 = \hat{f}_2 + \hat{f}_{-2} = 1/2$, $\hat{f}_{4k+1}^4 = \hat{f}_{4k+3}^4 = 0$.

For $n \geq 5$: $\hat{f}_k^n = \hat{f}_k$ (no aliasing), so $\hat{f}_0^n = 1/2$, $\hat{f}_{\pm 2}^n = 1/4$, and all others zero (modulo $n$).

---

### Solution 7

**(a)**

$$T_0(x) = 1, \quad T_1(x) = x, \quad T_2(x) = 2x^2 - 1, \quad T_3(x) = 4x^3 - 3x, \quad T_4(x) = 8x^4 - 8x^2 + 1.$$

**(b)** Since $V^{-1} = \frac{1}{n}DV^\top$ where $D = \text{diag}(1, 2, \ldots, 2)$, we can write $V^{-1} = D \cdot C_n$ where $C_n = \frac{1}{n}DV^\top$ is a scaled DCT matrix. The DCT $C_n$ is orthogonal, so $\|C_n\| = 1$. Then:

$$\|V^{-1}\| = \|D C_n \cdot D^{-1} \cdot D\| \leq \|D\| \|C_n\| = \sqrt{2/n} \cdot 1 = \sqrt{2/n}.$$

More precisely, $V^{-1} = \frac{1}{n}DV^\top$ and the matrix $\sqrt{D}V^\top / \sqrt{n}$ is orthogonal (the DCT), so $V^{-1} = \frac{1}{\sqrt{n}}\sqrt{D} \cdot (\sqrt{D}V^\top/\sqrt{n})$, giving $\|V^{-1}\| \leq \|\sqrt{D}/\sqrt{n}\| = \sqrt{2/n}$.

For the condition number: $\kappa(V) = \|V\|\|V^{-1}\|$. Since $V = C_n^{-1}D^{-1} = C_n^\top D^{-1}$, we have $\|V\| \leq \|D^{-1}\|\|C_n^\top\| = \sqrt{n} \cdot 1$... but more carefully, $\kappa(V) \leq \sqrt{2}$. The bounded condition number means errors in the function values are amplified by at most a constant factor when computing Chebyshev coefficients, regardless of $n$. This is in stark contrast to monomial Vandermonde matrices whose condition numbers grow exponentially.

---

### Solution 8

**(a)** Let $x = \cos\theta$. Then $U_n(\cos\theta) = \frac{\sin(n+1)\theta}{\sin\theta}$. Check: $U_0(\cos\theta) = \frac{\sin\theta}{\sin\theta} = 1$. For the leading coefficient: from the recurrence, $U_0 = 1$, $U_1 = 2x$, inductively $U_n = 2^n x^n + \cdots$. ✓

Three-term recurrence: using $2\cos\alpha\sin\beta = \sin(\alpha+\beta) - \sin(\alpha-\beta)$ with $\alpha = \theta$, $\beta = (n+1)\theta$:

$$x U_n(\cos\theta) = \cos\theta \cdot \frac{\sin(n+1)\theta}{\sin\theta} = \frac{\sin(n+2)\theta + \sin(n\theta)}{2\sin\theta} = \frac{U_{n+1}(\cos\theta) + U_{n-1}(\cos\theta)}{2}.$$

For $n = 0$: $xU_0 = x = \frac{\sin 2\theta}{2\sin\theta} = U_1/2$. ✓

**(b)** We have $T_{n+1}(x) = \cos(n+1)\theta$ with $x = \cos\theta$. Differentiating with respect to $x$:

$$T_{n+1}'(x) = \frac{d}{dx}\cos(n+1)\theta = -(n+1)\sin(n+1)\theta \cdot \frac{d\theta}{dx}.$$

Since $x = \cos\theta$, we get $dx = -\sin\theta\, d\theta$, so $\frac{d\theta}{dx} = \frac{-1}{\sin\theta}$.

Therefore:

$$T_{n+1}'(x) = (n+1)\frac{\sin(n+1)\theta}{\sin\theta} = (n+1)U_n(x). \qquad \square$$

**(c)** We need to show $\int_{-1}^1 C_n(x)C_m(x)(1-x^2)^{3/2}\,dx = 0$ for $n \neq m$, where $C_n = U_{n+1}'$.

Integration by parts:

$$\int_{-1}^1 C_n C_m (1-x^2)^{3/2}\,dx = \left[U_{n+1} C_m (1-x^2)^{3/2}\right]_{-1}^1 - \int_{-1}^1 U_{n+1} \frac{d}{dx}\left[C_m(1-x^2)^{3/2}\right]dx$$

The boundary term vanishes since $(1-x^2)^{3/2} = 0$ at $x = \pm 1$. Now:

$$\frac{d}{dx}\left[C_m(1-x^2)^{3/2}\right] = \left[C_m'(1-x^2) - 3xC_m\right]\sqrt{1-x^2}$$

The expression in brackets is a polynomial of degree $m + 1$. Since $m < n$ implies $m + 1 \leq n < n + 1$, and $U_{n+1}$ is orthogonal to all polynomials of degree $< n+1$ with respect to $\sqrt{1-x^2}$ (from the revision sheet Q7a), the integral vanishes. The case $m > n$ follows by swapping the integration by parts. $\square$

---

### Solution 9

**(a)** Since $w(-x) = 1-(-x)^2 = 1-x^2 = w(x)$, the weight is even. For even weights on $[-1,1]$, the recurrence coefficients $a_n = \frac{\langle xq_n, q_n\rangle}{\langle q_n, q_n\rangle}$. But $xq_n(x)^2 w(x)$ is odd (product of odd function $x$ with even functions $q_n^2$ and $w$), so $\langle xq_n, q_n\rangle = \int_{-1}^1 x q_n(x)^2 w(x)\,dx = 0$. Hence $a_n = 0$.

**(b)** $q_0(x) = k_0 = \frac{2!}{2 \cdot 1! \cdot 0!} = 1$.

$\|q_0\|^2 = \int_{-1}^1 (1-x^2)\,dx = 4/3$.

The recurrence gives $xq_0(x) = b_0 q_1(x)$, so $q_1(x) = x/b_0$.

We need $q_1$ to have leading coefficient $k_1 = \frac{4!}{2^2 \cdot 2! \cdot 1!} = \frac{24}{8} = 3$, so $1/b_0 = 3$, giving $b_0 = 1/3$.

Check: $\|q_1\|^2 = 9\int_{-1}^1 x^2(1-x^2)\,dx = 9 \cdot 4/15 = 12/5$.

For $q_2$: $xq_1(x) = b_0 q_0(x) + b_1 q_2(x)$, so $3x^2 = \frac{1}{3} + b_1 q_2(x)$.

We need $q_2$ with leading coefficient $k_2 = \frac{6!}{2^3 \cdot 3! \cdot 2!} = \frac{720}{48} = 15$, so $b_1 = 3/(15) = 1/5$.

Hence $q_2(x) = 5(3x^2 - 1/3) = 15x^2 - 5/3$.

**(c)** The Jacobi matrix is:

$$J_2 = \begin{bmatrix} 0 & 1/3 \\ 1/3 & 0 \end{bmatrix}$$

Eigenvalues: $\lambda^2 = 1/9$, so $\lambda = \pm 1/3$. The Gauss quadrature nodes are $x_1 = -1/3$ and $x_2 = 1/3$.

For weights: the eigenvectors are $\frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ \pm 1 \end{pmatrix}$, so $w_j = \|w\|_1 (e_1^\top Q e_j)^2 = \frac{4}{3} \cdot \frac{1}{2} = \frac{2}{3}$ for both $j = 1, 2$.

The 2-point rule is: $\int_{-1}^1 f(x)(1-x^2)\,dx \approx \frac{2}{3}f(-1/3) + \frac{2}{3}f(1/3)$.

---

### Solution 10

**(a)** In $F_{16}$, machine epsilon is $\epsilon_m = 2^{-10}$. Near 2, the exponent is $q - 15$ where $2 = 2^1$, so $q = 16$. The spacing between consecutive floats in $[2,4)$ is $2^{1-10} = 2^{-9}$.

- $\text{fl}(2) = 2$ — exact.
- $2 + 2^{-10} = 2 \times (1.0000000001)_2$ — but $2^{-10}$ is half the spacing $2^{-9}$, so this is exactly the midpoint between $2$ and $2 + 2^{-9}$. Round to even: $2 = 2 \times (1.0000000000)_2$ has last bit 0, so $\text{fl}(2 + 2^{-10}) = 2$.
- $2 + 2^{-9} + 2^{-10}$: this is the midpoint between $2 + 2^{-9}$ and $2 + 2^{-8}$. The number $2 + 2^{-9} = 2 \times (1.000000001)_2$ has last bit 1; $2 + 2^{-8} = 2 \times (1.00000001\mathbf{0})_2$ has last bit 0. Round to even: $\text{fl}(2 + 2^{-9} + 2^{-10}) = 2 + 2^{-8}$.

**(b)** In $[2,4)$, the exponent satisfies $2^{q-15} = 2$, so the spacing is $2 \cdot \epsilon_m = 2^{1-10} = 2^{-9}$. In $[1,2)$, the spacing is $\epsilon_m = 2^{-10}$. The spacing doubles when crossing from $[1,2)$ to $[2,4)$.

**(c)** Standard proof: write $x = 2^{q-\sigma}(1.b_1\ldots b_S b_{S+1}\ldots)_2$. The nearest floats are $x^- = 2^{q-\sigma}(1.b_1\ldots b_S)_2$ and $x^+ = x^- + 2^{q-\sigma-S}$. The gap is $2^{q-\sigma-S}$, so $|x - \text{fl}(x)| \leq 2^{q-\sigma-S-1}$. Dividing by $|x| \geq 2^{q-\sigma}$ gives $|\delta_x| \leq 2^{-S-1} = \epsilon_m/2$. $\square$