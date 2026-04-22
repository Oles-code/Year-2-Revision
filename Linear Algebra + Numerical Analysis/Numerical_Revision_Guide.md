# MATH50003 Numerical Analysis — Revision Guide

*Based on the lecture notes, 2023–2025 exam papers, and problem sheets.*

---

## How to Use This Guide

The numerical analysis section of the exam typically consists of **3 questions** (Q4, Q5, Q6), each worth 20 marks. From past papers, the topics break down very consistently:

- **Q4:** Floating point arithmetic, divided differences error bounds, dual numbers (every year)
- **Q5:** Matrix factorisations (LU/PLU/Cholesky) + reflections/rotations + backward error analysis (every year)
- **Q6:** Fourier series/DFT + orthogonal polynomials + Gauss quadrature (every year)

The exam solutions label each part as **seen** (bookwork/similar to notes), **sim. seen** (similar method to notes/sheets), **meth seen** (method seen but applied differently), or **unseen** (genuinely new). The hardest marks come from error bound derivations and the orthogonal polynomials question.

---

## Part I: Calculus on a Computer

### 1.1 Rectangular Rule

The right-sided rectangular rule approximates an integral by:

$$\int_a^b f(x)\,dx \approx h \sum_{j=1}^{n} f(x_j), \qquad x_j = a + jh, \quad h = \frac{b-a}{n}$$

**Error on one panel** (proved by integration by parts):

$$\int_a^b f(x)\,dx = (b-a)f(b) + \delta, \qquad |\delta| \leq M(b-a)^2$$

where $M = \sup_{a \leq x \leq b} |f'(x)|$.

**Overall error:** Summing over $n$ panels gives $O(1/n)$ convergence — the error is bounded by $\frac{M(b-a)^2}{n}$.

**Proof technique (important to know):** The key trick is writing $\int_a^b f(x)\,dx = \int_a^b (x-a)' f(x)\,dx$ and integrating by parts. This avoids Taylor series entirely.

**Trapezium rule** (from problem sheet 1): Converges at rate $O(h^2) = O(1/n^2)$, which is a quadratic improvement. For periodic functions, the trapezium rule converges *extremely* fast (faster than any polynomial rate), which is why it reappears in the Fourier series chapter.

---

### 1.2 Divided Differences

Approximating derivatives:

$$f'(x) \approx \frac{f(x+h) - f(x)}{h}$$

**Error bound (via Taylor series):**

$$f'(x) = \frac{f(x+h) - f(x)}{h} + \delta, \qquad |\delta| \leq \frac{Mh}{2}$$

where $M = \sup_{x \leq t \leq x+h} |f''(t)|$.

**Central differences** (from sheet 1, examined 2024): $f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}$, which has error $O(h^2)$ — you need $f$ to be *thrice*-differentiable and use the third derivative in the bound:

$$|\delta| \leq \frac{h^2}{6} \sup_{x-h \leq t \leq x+h} |f'''(t)|$$

**The mystery:** In practice, divided differences *diverge* as $h \to 0$ on a computer. This is because of floating point rounding errors — resolved in Section II.

**Second derivative approximation** (from sheet 1):
$$f''(x) \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^2}$$

---

### 1.3 Dual Numbers

**Definition:** Dual numbers $\mathbb{D} = \{a + b\epsilon : a,b \in \mathbb{R}, \; \epsilon^2 = 0\}$.

Multiplication: $(a + b\epsilon)(c + d\epsilon) = ac + (bc + ad)\epsilon$.

Compare with complex numbers: same structure but $i^2 = -1$ vs $\epsilon^2 = 0$.

**Key theorem (polynomials):** For any polynomial $p$,
$$p(a + b\epsilon) = p(a) + bp'(a)\epsilon$$

*Proof by induction:* The cases $n=0,1$ are immediate. For $n > 1$:
$(a+b\epsilon)^n = (a+b\epsilon)(a+b\epsilon)^{n-1} = (a+b\epsilon)(a^{n-1} + (n-1)ba^{n-2}\epsilon) = a^n + bna^{n-1}\epsilon$.
Linearity then extends to general polynomials.

**Dual extensions of standard functions:** For a differentiable function $f$, define:
$$f(a + b\epsilon) := f(a) + bf'(a)\epsilon$$

So for example:
- $\exp(a + b\epsilon) = \exp(a) + b\exp(a)\epsilon$
- $\sin(a + b\epsilon) = \sin(a) + b\cos(a)\epsilon$
- $\cos(a + b\epsilon) = \cos(a) - b\sin(a)\epsilon$
- $\log(a + b\epsilon) = \log(a) + \frac{b}{a}\epsilon$

**Why this works:** Addition and multiplication of dual extensions recover the product rule (Lemma 2), and composition recovers the chain rule (Lemma 3). So evaluating *any* function built from basic operations on $a + \epsilon$ automatically gives $f(a) + f'(a)\epsilon$.

**Practical usage:** To compute $f'(a)$, evaluate $f(a + \epsilon)$ and read off the $\epsilon$-coefficient.

**Exam variants — know these well:**

1. **2nd order dual numbers** (2023 Q4c): $x + a\epsilon_1 + b\epsilon_2$ with $\epsilon_1^2 = \epsilon_2$, $\epsilon_2^2 = \epsilon_1\epsilon_2 = 0$. Then $p(x + a\epsilon_1 + b\epsilon_2) = p(x) + ap'(x)\epsilon_1 + (bp'(x) + \frac{a^2}{2}p''(x))\epsilon_2$. Choose $b = 0, a = 1$ to get $p''(x)$ from the $\epsilon_2$ coefficient.

2. **Double-dual numbers** (2024 Q4c): $a + b\epsilon + c\delta + d\epsilon\delta$ with $\epsilon^2 = \delta^2 = 0$ and $\epsilon\delta$ commutes. Product formula: $(a+b\epsilon+c\delta+d\epsilon\delta)(e+f\epsilon+g\delta+h\epsilon\delta) = ae + (af+be)\epsilon + (ag+ce)\delta + (ah+bg+cf+de)\epsilon\delta$. For polynomials: $p(a+b\epsilon+c\delta+d\epsilon\delta) = p(a) + bp'(a)\epsilon + cp'(a)\delta + (dp'(a) + bcp''(a))\epsilon\delta$. Choose $b=c=1, d=0$ to get $p''(a)$ from the $\epsilon\delta$ coefficient.

---

### 1.4 Newton's Method

Given $f(r) = 0$, iterate: $x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}$.

**Error bound (quadratic convergence):** If $\epsilon_k = r - x_k$, then
$$|\epsilon_{k+1}| \leq M|\epsilon_k|^2, \qquad M = \frac{1}{2}\sup_{x \in B} |f''(x)| \cdot \sup_{x \in B} \left|\frac{1}{f'(x)}\right|$$

*Proof:* Use Taylor's theorem: $0 = f(r) = f(x_k) + f'(x_k)\epsilon_k + \frac{f''(t)}{2}\epsilon_k^2$. Rearranging: $\epsilon_{k+1} = \epsilon_k + \frac{f(x_k)}{f'(x_k)} = -\frac{f''(t)}{2f'(x_k)}\epsilon_k^2$.

**Convergence:** If $|ε_0| < M^{-1}$ then $|ε_k| \leq \frac{1}{M}(M|ε_0|)^{2^k} \to 0$. The number of correct digits roughly *doubles* each iteration.

**Degenerate case** (from sheet 2): When $f'(r) = 0$ too, convergence slows to linear ($|\epsilon_{k+1}| \leq c|\epsilon_k|$ for some $c < 1$).

---

## Part II: Representing Numbers

### 2.1 Integers and Reals in Binary

**Binary format:** $\pm(B_p \ldots B_1 B_0)_2 = \pm \sum_{k=0}^{p} B_k 2^k$.

**Unsigned 8-bit (UInt8):** Represents $0$ to $255$. The bits $b_7 b_6 \ldots b_0$ represent $\sum_{k=0}^7 b_k 2^k$.

**Signed 8-bit (Int8) — two's complement:** Represents $-128$ to $127$. The value is $-b_7 \cdot 2^7 + \sum_{k=0}^6 b_k 2^k$.

*Example (examined 2023):* The bits `11111110`:
- **UInt8:** $128 + 64 + 32 + 16 + 8 + 4 + 2 = 254$
- **Int8:** $-128 + 64 + 32 + 16 + 8 + 4 + 2 = -2$

### 2.2 Floating Point Numbers

**IEEE format:** A floating point number in $F_{\sigma,Q,S}$ is stored as: **1 sign bit | Q exponent bits | S significand bits**.

The value is: $(-1)^s \times 2^{q - \sigma} \times (1.b_1 b_2 \ldots b_S)_2$

where $q$ is the stored exponent (an unsigned integer from $1$ to $2^Q - 2$), and $\sigma = 2^{Q-1} - 1$ is the exponent bias.

**Key constants:**
- Machine epsilon: $\epsilon_m = 2^{-S}$ (gap between 1 and the next float)
- Smallest positive normal: $2^{1-\sigma}$
- Largest normal: $2^{2^Q - 2 - \sigma}(2 - \epsilon_m)$

**IEEE 16-bit ($F_{16}$):** $\sigma = 15$, $Q = 5$, $S = 10$.

*Example (examined 2023):* Bits `1 00010 1100000000`:
- Sign: $s = 1$ (negative)
- Exponent: $q = (00010)_2 = 2$, so $q - \sigma = 2 - 15 = -13$
- Significand: $(1.1100000000)_2 = 1 + 1/2 + 1/4 = 7/4$
- Value: $-2^{-13} \times 7/4 = -7/2^{15} = -7/32768$

**Special values:** $q = 0$ gives sub-normals (including $\pm 0$); $q = 2^Q - 1$ gives $\pm\infty$ or NaN.

*Example (examined 2024):* Rounding $1, 1+2^{-11}, 1+2^{-10}+2^{-11}$ to $F_{16}$:
- $1 = (1.0000000000)_2$ — exact float
- $1 + 2^{-11} = (1.00000000001)_2$ — the 11th bit is beyond $S=10$ significant bits, so it's at the halfway point. By "round to even" rule (since $b_{10} = 0$), it rounds *down* to $1$.
- $1 + 2^{-10} + 2^{-11} = (1.0000000011)_2$ — this fits exactly in 10 significant bits.

---

### 2.3 Floating Point Arithmetic

**Core principle:** IEEE operations are *exact up to rounding*:
$$x \oplus y := \mathrm{fl}(x + y), \quad x \ominus y := \mathrm{fl}(x - y), \quad x \otimes y := \mathrm{fl}(x \times y), \quad x \oslash y := \mathrm{fl}(x / y)$$

**WARNING:** These operations are *not associative*! $(x \oplus y) \oplus z \neq x \oplus (y \oplus z)$ in general.

**Round bound (Proposition 2):** If $x$ is in the normalised range, then
$$\mathrm{fl}(x) = x(1 + \delta_x), \qquad |\delta_x| \leq \frac{\epsilon_m}{2} \text{ (nearest rounding)}$$

*Proof sketch:* Write $x = 2^{q-\sigma}(1.b_1 \ldots b_S b_{S+1} \ldots)_2$. The nearest floats are $x_- = 2^{q-\sigma}(1.b_1 \ldots b_S)_2$ and $x_+ = x_- + 2^{q-\sigma-S}$. The maximum rounding error is half the gap: $|x - \mathrm{fl}(x)| \leq 2^{q-\sigma-S-1}$. Dividing by $|x| \geq 2^{q-\sigma}$ gives $|\delta_x| \leq 2^{-S-1} = \epsilon_m/2$.

This immediately implies: $x \oplus y = (x+y)(1+\delta_1)$ where $|\delta_1| \leq \epsilon_m/2$.

**Idealised floating point $F_{\infty,S}$:** We extend to all exponents (no overflow/underflow), making the round bound valid for all reals. This simplifies error analysis.

### Bounding Errors in Practice — The Exam Technique

This is the most commonly examined skill. The pattern is always the same:

1. **Replace each constant** with $\mathrm{fl}(\text{value}) = \text{value} \times (1 + \delta_i)$
2. **Replace each operation** $\oplus, \otimes$ etc. with the exact result times $(1 + \delta_j)$
3. **Expand** and collect terms
4. **Bound** using $|\delta_i| \leq \epsilon_m/2$
5. **Simplify** by bounding $\epsilon_m^2$ terms: $|\delta_i \delta_j| \leq \epsilon_m^2/4 \leq \epsilon_m/4$
6. **Round up** constants to integers for cleaner bounds

*Worked example from notes:* $(1.1 + 1.2) \times 1.3 = 2.99$ on a computer:

$(\mathrm{fl}(1.1) \oplus \mathrm{fl}(1.2)) \otimes \mathrm{fl}(1.3) = (1.1(1+\delta_1) + 1.2(1+\delta_2))(1+\delta_3) \times 1.3(1+\delta_4)(1+\delta_5)$

Expanding and bounding carefully gives absolute error $|\delta| \leq 23\epsilon_m$.

### Divided Differences with Floating Point (Theorem 4)

**This is examined every single year.** Given $f^{FP}(x) = f(x) + \delta_x^f$ with $|\delta_x^f| \leq c\epsilon_m$:

$$\frac{f^{FP}(x+h) \ominus f^{FP}(x)}{h} = f'(x) + \delta^{FD}_{x,h}$$

where

$$|\delta^{FD}_{x,h}| \leq \frac{|f'(x)|}{2}\epsilon_m + Mh + \frac{4c\epsilon_m}{h}$$

**The three terms tell a story:**
1. $\frac{|f'(x)|}{2}\epsilon_m$ — fixed small error from rounding the division
2. $Mh$ — Taylor series truncation, decreases as $h \to 0$
3. $\frac{4c\epsilon_m}{h}$ — floating point error amplification, *grows* as $h \to 0$

**Heuristic:** Choose $h \propto \sqrt{\epsilon_m}$ to balance terms 2 and 3. For double precision, $\sqrt{\epsilon_m} \approx 1.5 \times 10^{-8}$.

*Proof technique:* Use Taylor expansion on $f(x+h)$, then bound $|1+\delta_1| \leq 2$ (pessimistic but clean).

**Exam variants seen:** 2023 asked for a *central* difference error bound with floating point; 2024 asked for central differences too. The structure is always the same — just the Taylor expansion step changes (you get $f'''$ instead of $f''$, and the algebra adjusts slightly).

---

### 2.4 Interval Arithmetic

For intervals $X = [a,b]$, $Y = [c,d]$ with $0 < a \leq b$, $0 < c \leq d$:
- $X + Y = [a+c, b+d]$
- $XY = [ac, bd]$
- $X/n = [a/n, b/n]$

**Floating point interval arithmetic:** Round endpoints in the "safe" direction:
- $[a,b] \oplus [c,d] = [\mathrm{fl}^{\text{down}}(a+c), \mathrm{fl}^{\text{up}}(b+d)]$

This guarantees that the true result is *contained* in the computed interval.

**Application:** Computing $e$ rigorously using Taylor series with error bounds. Combine interval arithmetic on the partial sum with an explicit bound on the Taylor remainder $|R_n| \leq \frac{3}{(n+1)!}$.

This topic appears on problem sheets but has not been directly examined in the final exams (2023–2025). Focus on understanding the concept rather than drilling computations.

---

## Part III: Numerical Linear Algebra

### 3.1 Structured Matrices

**Dense:** $m \times n$ matrix-vector multiplication costs $O(mn)$ operations.

**Triangular:** Forward/back-substitution solves $Lx = b$ or $Ux = b$ in $O(n^2)$ operations.

**Banded:** A matrix with bandwidths $(l, u)$ has $a_{kj} = 0$ when $j - k > u$ or $k - j > l$.
- **Bidiagonal:** $l = 1, u = 0$ (or vice versa). Multiplication and solve in $O(n)$.
- **Tridiagonal:** $l = u = 1$. Multiplication and solve in $O(n)$.

**Floating point error in matrix-vector multiplication** (from sheet 4 — very relevant to exam Q5): If $b = \text{mul}(A, x)$ is matrix-vector multiplication in floating point, then $Ax = b + \delta$ where $\|\delta\|_\infty \leq \frac{3}{2}\epsilon_m \|A\|_\infty \|x\|_\infty$ (for bidiagonal; coefficients change for other structures). The infinity norm of a matrix is $\|A\|_\infty = \max_k \sum_j |a_{kj}|$.

---

### 3.2 LU Factorisation

**What it is:** Write $A = LU$ where $L$ is lower triangular (with 1s on the diagonal) and $U$ is upper triangular. This is just Gaussian elimination recorded as matrices.

**The recursive construction:** Write $A = \begin{pmatrix} \alpha & \mathbf{w}^\top \\ \mathbf{v} & K \end{pmatrix}$, then:

$$A = \underbrace{\begin{pmatrix} 1 \\ \mathbf{v}/\alpha & I \end{pmatrix}}_{L_1} \begin{pmatrix} \alpha & \mathbf{w}^\top \\ & K - \mathbf{v}\mathbf{w}^\top/\alpha \end{pmatrix}$$

Continue recursively on $A_2 = K - \mathbf{v}\mathbf{w}^\top/\alpha$.

**How to compute LU by hand:** This is the bread-and-butter exam technique.

1. Take the first column of $A$: $\alpha_1 = a_{11}$, $\mathbf{v}_1 = $ rest of column 1
2. First column of $L$: entries are $\mathbf{v}_1/\alpha_1$ below the 1
3. First row of $U$: is the first row of $A$
4. Compute $A_2 = K - \mathbf{v}_1 \mathbf{w}_1^\top/\alpha_1$ (the Schur complement)
5. Repeat on $A_2$

*Example (from 2025 exam solution):*
$$A = \begin{pmatrix} 1 \\ 0 & 1 \\ 3 & 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 2 & 0 \\ & 5 & 1 \\ & & -5+2 \end{pmatrix}$$

---

### 3.3 PLU Factorisation

When a diagonal entry is zero (or for stability, when a larger entry exists below it), we **pivot** — swap rows using a permutation matrix.

**Result:** Every invertible matrix has a PLU factorisation: $A = P^\top L U$, where $P$ is a permutation matrix.

**Permutation matrices:** $P_\sigma$ permutes the entries of a vector according to $\sigma$. Key property: $P_\sigma^{-1} = P_\sigma^\top$ (permutation matrices are orthogonal).

**How to compute PLU by hand:**
1. Find the largest (in magnitude) entry in the first column
2. Swap that row to the top (record in $P$)
3. Proceed with LU on the permuted matrix
4. At each step, pivot again if needed

**Exam tip (2023 examiner comments):** Many students confused $P$ with $P^{-1}$. Remember: $A = P^\top L U$, so $PA = LU$.

---

### 3.4 Cholesky Factorisation

**For symmetric positive definite (SPD) matrices only:** $A = LL^\top$ where $L$ is lower triangular with positive diagonal entries.

**Key theorem:** A matrix is SPD if and only if it has a Cholesky factorisation.

**Properties of positive definite matrices:**
- All diagonal entries are positive ($a_{kk} = e_k^\top A e_k > 0$)
- $V^\top A V$ is positive definite if $V$ is non-singular
- Any principal submatrix is positive definite

**How to compute Cholesky by hand:**

Write $A = \begin{pmatrix} \alpha & \mathbf{v}^\top \\ \mathbf{v} & K \end{pmatrix}$. Then $L_1 = \begin{pmatrix} \sqrt{\alpha} \\ \mathbf{v}/\sqrt{\alpha} & I \end{pmatrix}$ and recursively factorise $A_2 = K - \mathbf{v}\mathbf{v}^\top/\alpha$.

*Examined every year (2023 Q5 Cholesky by hand, 2024 Q5b Cholesky by hand, 2025 Q5 LU by hand). This is essentially free marks — drill it.*

---

### 3.5 Orthogonal and Unitary Matrices

**Definition:** $Q \in O(n)$ if $Q^\top Q = I$; $Q \in U(n)$ if $Q^* Q = I$.

**Key properties:**
1. Norm-preserving: $\|Qx\| = \|x\|$
2. Eigenvalues have $|\lambda| = 1$
3. $\det Q = \pm 1$ (real case)
4. Trivially invertible: $Q^{-1} = Q^\top$ (or $Q^*$)
5. Stable for numerical computation (errors don't amplify)

### Rotations

A 2D rotation by angle $\theta$:
$$Q_\theta = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$$

has $\det Q_\theta = 1$. To rotate a vector $\mathbf{x} = \begin{pmatrix} a \\ b \end{pmatrix}$ to align with $e_1$ (i.e., $Q\mathbf{x} = \|\mathbf{x}\|e_1$), choose $\cos\theta = a/\|\mathbf{x}\|$, $\sin\theta = b/\|\mathbf{x}\|$:

$$Q = \frac{1}{\sqrt{a^2+b^2}}\begin{pmatrix} a & b \\ -b & a \end{pmatrix}$$

The condition for this to be a rotation (not a reflection) is $\det Q = 1$.

**Examined 2025 Q5b:** Given a matrix, identify which rotation matrix $Q$ has $\det Q = 1$.

### Reflections

**Definition:** $Q_v = I - 2vv^\top$ where $\|v\| = 1$. This reflects through the hyperplane orthogonal to $v$.

Properties: $Q_v$ is symmetric ($Q_v^\top = Q_v$) and orthogonal ($Q_v^2 = I$).

**Key fact** (examined 2024 Q5c): Any diagonal matrix $D_k$ with one $-1$ on the diagonal and all other entries $+1$ can be written as a reflection. Choose $v = e_k$, then $I - 2e_k e_k^\top = D_k$.

### Householder Reflections

Given $\mathbf{x} \in \mathbb{R}^n$, the Householder reflection maps $\mathbf{x}$ to $\pm\|\mathbf{x}\|e_1$:

$$Q_x^H = I - 2\mathbf{w}\mathbf{w}^\top, \qquad \mathbf{y} = \mp\|\mathbf{x}\|e_1 + \mathbf{x}, \quad \mathbf{w} = \frac{\mathbf{y}}{\|\mathbf{y}\|}$$

**Default choice of sign:** $Q_x^H \mathbf{x} = -\text{sign}(x_1)\|\mathbf{x}\|e_1$ (opposite sign to $x_1$ for stability).

**Variant (2023 Q5b):** Choose $v$ so that $Q_v x = \|x\|e_n$ (map to *last* basis vector instead of first). Then $\mathbf{y} = \mathbf{x} - \|\mathbf{x}\|e_n$.

**Key application:** Householder reflections are used to prove QR factorisation exists (by induction).

---

### 3.6 QR Factorisation

**Theorem:** Every $A \in \mathbb{C}^{m \times n}$ has a QR factorisation: $A = QR$ where $Q \in U(m)$ and $R$ is right-triangular (upper triangular if square).

**Reduced QR:** $A = \hat{Q}\hat{R}$ where $\hat{Q} \in \mathbb{C}^{m \times n}$ has orthonormal columns and $\hat{R} \in \mathbb{C}^{n \times n}$ is upper triangular.

**Two methods to compute QR:**

1. **Gram–Schmidt** (triangular orthogonalisation): Orthogonalise columns of $A$ one by one. Gives the *reduced* QR.

2. **Householder reflections** (orthogonal triangularisation): Apply Householder reflections to introduce zeros column by column. The proof is by induction (same structure as LU/Cholesky proofs): apply $Q_1 = Q_{a_1}^H$ to get $Q_1 A = \begin{pmatrix} \alpha_1 & \mathbf{w}_1^\top \\ & A_2 \end{pmatrix}$, then recurse on $A_2$.

**Note from the notes:** Computing QR by hand via Householder is very hard to do with nice numbers. The notes say this is "not examinable" for Householder QR by hand, but the *proof technique* (induction via Householder) is examinable and has been tested (2023 Q5b: prove $A = QL$ factorisation; 2025 Q5c: prove LQ factorisation).

### QR and Least Squares

Given $A \in \mathbb{C}^{m \times n}$ with $m \geq n$ and full column rank, the vector $\mathbf{x} = \hat{R}^{-1}\hat{Q}^*\mathbf{b}$ minimises $\|A\mathbf{x} - \mathbf{b}\|$.

*Proof:* $\|A\mathbf{x} - \mathbf{b}\| = \|QR\mathbf{x} - \mathbf{b}\| = \|R\mathbf{x} - Q^*\mathbf{b}\|$ (since $Q$ is norm-preserving). The bottom $m-n$ rows of $R$ are zero, giving a fixed contribution. Minimise by setting the top $n$ rows to zero: $\hat{R}\mathbf{x} = \hat{Q}^*\mathbf{b}$.

---

### 3.7 Backward Error Analysis

**This appears on every exam Q5.** The idea: instead of asking "how far is our computed answer from the true answer?" (forward error), ask "for what *perturbed problem* does our computation give the *exact* answer?" (backward error).

**Pattern (examined 2023 Q5c, 2024 Q5a):** For a bidiagonal system $Lx = b$ solved with floating point forward substitution:

$x_1 = \mathrm{fl}(b_1) = b_1(1+\delta_1)$

$x_k = (b_k \ominus (\ell_{k-1} \otimes x_{k-1})) = b_k(1+\delta_k') - \ell_{k-1}x_{k-1}(1+\delta_k'')$

Then we can write $\text{bidiagforwardsub}(L, b) = (L + \delta L)^{-1}(b + \delta b)$ where the perturbations satisfy bounds like $|\delta b_k| \leq \frac{\epsilon_m}{2}|b_k|$ and $|\delta \ell_k| \leq 2\epsilon_m |\ell_k|$.

**For matrix-vector multiplication** (2024 Q5a): $Ux = \text{bidiagmul}(U,x) + \delta$ where $\|\delta\|_\infty \leq \frac{3}{2}\epsilon_m \|U\|_\infty \|x\|_\infty$.

---

## Part IV: Linear Algebra Applications

### 4.1 Polynomial Interpolation

Given $n$ points $(x_j, f_j)$, find a polynomial $p$ of degree $< n$ passing through all points: $p(x_j) = f_j$.

This is equivalent to solving the **Vandermonde system**:
$$V\mathbf{c} = \mathbf{f}, \qquad V_{jk} = x_j^{k-1}$$

**Warning:** Interpolation at evenly spaced points with monomials is unstable for large $n$.

### 4.2 Polynomial Regression (Least Squares)

When $m > n$ (more data points than polynomial degree), solve $\min \|V\mathbf{c} - \mathbf{f}\|$ via the QR factorisation of $V$.

### 4.3 Singular Value Decomposition (SVD)

**Definition:** $A = U\Sigma V^*$ where $U, V$ are unitary and $\Sigma$ is diagonal with non-negative entries $\sigma_1 \geq \sigma_2 \geq \ldots \geq \sigma_r > 0$ (singular values).

**Existence:** Follows from diagonalising the Gram matrix $A^*A = Q\Lambda Q^*$. The eigenvalues $\lambda_k \geq 0$ give $\sigma_k = \sqrt{\lambda_k}$, and $V = Q$, $U$ columns are $Av_k/\sigma_k$.

**Key property:** $\ker(A) = \ker(A^*A)$ (Gram matrix kernel proposition).

**Matrix 2-norm:** $\|A\|_2 = \sigma_1$ (the largest singular value).

**Best rank-$k$ approximation (Eckart–Young theorem):** $A_k = \sum_{j=1}^k \sigma_j u_j v_j^*$ minimises $\|A - B\|$ over all rank-$k$ matrices $B$. The error is $\|A - A_k\| = \sigma_{k+1}$.

**Pseudoinverse** (from sheet 7): $A^+ = V\Sigma^{-1}U^*$ gives the least-squares solution $x = A^+ b$.

---

## Part V: Numerical Fourier Series

### 5.1 Fourier Expansions

A function $f$ has a Fourier expansion:
$$f(\theta) = \sum_{k=-\infty}^{\infty} \hat{f}_k e^{ik\theta}, \qquad \hat{f}_k = \frac{1}{2\pi}\int_0^{2\pi} e^{-ik\theta} f(\theta)\,d\theta$$

**Convergence:** If $f$ is $p$ times differentiable and periodic, then $|\hat{f}_k| = O(|k|^{-p})$. If $f$ is smooth and periodic, coefficients decay faster than any polynomial.

**Discrete Fourier coefficients** (via the trapezium rule):
$$\hat{f}_k^n = \frac{1}{n}\sum_{j=0}^{n-1} f(\theta_j) e^{-ik\theta_j}, \qquad \theta_j = \frac{2\pi j}{n}$$

**Aliasing:** $\hat{f}_k^n = \sum_{p=-\infty}^{\infty} \hat{f}_{k+pn}$ (the discrete coefficients "fold in" contributions from frequencies shifted by multiples of $n$). In particular, $\hat{f}_k^n = \hat{f}_{k+pn}^n$ for all integers $p$.

**Discrete orthogonality (key lemma):**
$$\sum_{j=0}^{n-1} e^{ik\theta_j} = \begin{cases} n & k = 0, \pm n, \pm 2n, \ldots \\ 0 & \text{otherwise} \end{cases}$$

*Proof:* For $k = pn$, every term is 1. Otherwise, use the geometric series: $\sum_{j=0}^{n-1} \omega^{kj} = \frac{\omega^{kn} - 1}{\omega^k - 1} = 0$ since $\omega^n = 1$ but $\omega^k \neq 1$.

*Example (examined 2023 Q6a):* For $f(\theta) = \cos 2\theta = \frac{1}{2}e^{-2i\theta} + \frac{1}{2}e^{2i\theta}$, so $\hat{f}_{-2} = \hat{f}_2 = 1/2$ and all others are zero. The discrete coefficients $\hat{f}_k^n$ follow from aliasing.

---

### 5.2 Discrete Fourier Transform (DFT)

**DFT matrix:** With $\omega = e^{2\pi i/n}$:

$$Q_n = \frac{1}{\sqrt{n}} \begin{pmatrix} 1 & 1 & 1 & \cdots & 1 \\ 1 & \omega^{-1} & \omega^{-2} & \cdots & \omega^{-(n-1)} \\ 1 & \omega^{-2} & \omega^{-4} & \cdots & \omega^{-2(n-1)} \\ \vdots & & & & \vdots \\ 1 & \omega^{-(n-1)} & \omega^{-2(n-1)} & \cdots & \omega^{-(n-1)^2} \end{pmatrix}$$

So $\hat{\mathbf{f}}^n = \frac{1}{\sqrt{n}} Q_n \mathbf{f}^n$ and $\mathbf{f}^n = \sqrt{n} Q_n^* \hat{\mathbf{f}}^n$.

**$Q_n$ is unitary** (examined 2025 Q6a): $Q_n Q_n^* = I$. Proof uses the discrete orthogonality lemma: $(Q_n Q_n^*)_{k\ell} = \frac{1}{n}\sum_{j=0}^{n-1} \omega^{(k-\ell)j}$, which equals 1 if $k = \ell$ and 0 otherwise.

**Trigonometric interpolation:** $f_n(\theta) = \sum_{k=0}^{n-1} \hat{f}_k^n e^{ik\theta}$ interpolates $f$ at $\theta_j$: $f_n(\theta_j) = f(\theta_j)$.

**Connection to polynomial interpolation** (examined 2024 Q6a): Interpolating $\exp(z)$ at $z = 1, i, -1, -i$ (i.e., $z = e^{i\theta_j}$ for $\theta_j = 0, \pi/2, \pi, 3\pi/2$) can be done via the DFT with $\omega = i$.

---

### 5.3 Orthogonal Polynomials

**Definition:** Polynomials $\{q_n\}$ are orthogonal with respect to a weight $w(x)$ on $[a,b]$ if:
$$\langle q_m, q_n \rangle = \int_a^b q_m(x) q_n(x) w(x)\,dx = 0 \quad \text{for } m \neq n$$

**Three-term recurrence:** Orthogonal polynomials always satisfy:
$$xq_n(x) = a_n q_n(x) + b_n q_{n+1}(x) + b_{n-1} q_{n-1}(x)$$

where $a_n = \frac{\langle xq_n, q_n \rangle}{\|q_n\|^2}$ and $b_n$ depends on normalisation.

**Jacobi matrix:** The three-term recurrence coefficients form a symmetric tridiagonal matrix:
$$X = \begin{pmatrix} a_0 & b_0 \\ b_0 & a_1 & b_1 \\ & b_1 & a_2 & \ddots \\ & & \ddots & \ddots \end{pmatrix}$$

**Symmetry property** (examined 2023 Q6c): If $w(x) = w(-x)$ (even weight), then $a_k = 0$ for all $k$ (even polynomials have even coefficients, odd have odd). This means the Jacobi matrix has zero diagonal.

### Key Families

1. **Legendre:** $P_n(x) = \frac{(2n)!}{2^n (n!)^2} x^n + O(x^{n-1})$, weight $w(x) = 1$ on $[-1,1]$

2. **Chebyshev (2nd kind):** $U_n(x) = 2^n x^n + O(x^{n-1})$, weight $w(x) = \sqrt{1-x^2}$ on $[-1,1]$.
   Key identity: $U_n(\cos\theta) = \frac{\sin(n+1)\theta}{\sin\theta}$

3. **Ultraspherical $C_n$:** Various weights of the form $(1-x^2)^\alpha$.
   2023: $w(x) = (1-x^2)^{3/2}$, normalisation $C_n(x) = 2^n n(n+1) x^n + \ldots$
   2024: weight $\sqrt{1-x^2}$ (Chebyshev)
   2025: $w(x) = 1 - x^2$, normalisation $C_n(x) = \frac{(2(n+1))!}{2^{n+1}(n+1)!n!}x^n + \ldots$

### Constructing Orthogonal Polynomials (Stieltjes procedure)

This is the main exam technique for Q6. Given a weight and normalisation:

1. Start with $C_0(x) = 1$ (or given normalisation)
2. Compute $\|C_0\|^2 = \int_{-1}^1 w(x)\,dx$
3. Use three-term recurrence: $xC_n(x) = a_n C_n(x) + $ terms, where $a_n = \frac{\langle xC_n, C_n \rangle}{\|C_n\|^2}$
4. If $w$ is even, $a_n = 0$ by symmetry (big shortcut!)
5. Solve for $C_{n+1}$ using the required leading coefficient

### Derivatives and Orthogonal Polynomials

A recurring exam pattern: show that the derivative of one family of orthogonal polynomials gives another family. The technique is:

1. **Check normalisation:** Verify the leading coefficient matches after differentiation
2. **Check orthogonality:** Use integration by parts. If $f_m$ has degree $< n$:
$$\int_{-1}^1 P'_{n+1}(x) f_m(x) w(x)\,dx = -\int_{-1}^1 P_{n+1}(x) [f_m(x)w(x)]'\,dx = 0$$
since $[f_m w]'$ has degree $\leq m+1 < n+1$ and $P_{n+1}$ is orthogonal to all lower degree polynomials w.r.t. weight 1. (The boundary terms vanish when $w$ vanishes at $\pm 1$.)

*Examined:* 2023 Q6b(ii) ($U'_{n+1} = 2C_n$), 2025 Q6b(ii) ($C_n = P'_{n+1}$).

---

### 5.4 Gauss Quadrature

**$n$-point Gauss quadrature** for weight $w$: Choose nodes $x_1, \ldots, x_n$ and weights $w_1, \ldots, w_n$ so that:
$$\int_{-1}^1 f(x) w(x)\,dx \approx \sum_{j=1}^n w_j f(x_j)$$

is *exact* for all polynomials of degree $\leq 2n-1$.

**The nodes** are the eigenvalues of the $n \times n$ truncated Jacobi matrix $J_n$ (equivalently, the roots of $q_n$).

**The weights** are $w_j = \|w\|_1 (e_1^\top Q e_j)^2$ where $J_n = Q \text{diag}(x_1, \ldots, x_n) Q^\top$ and $\|w\|_1 = \int w(x)\,dx$.

**Symmetry** (examined 2023 Q6c(ii)): If $w$ is even, then the Gauss quadrature points come in symmetric pairs: if $x_k$ is a node, so is $-x_k$.

**Computing 2-point Gauss quadrature** (examined 2025 Q6b(iii)):
1. Find the roots of $q_2(x)$ using the recurrence
2. The 2×2 Jacobi matrix $J_2 = \begin{pmatrix} 0 & b_0 \\ b_0 & 0 \end{pmatrix}$ (if $w$ is even) has eigenvalues $\pm b_0$
3. Compute weights from the eigenvectors

---

## Exam Technique Summary

### What comes up every year

| Topic | Typical marks | Difficulty |
|-------|------|------|
| Integer/float representation | 3–4 | Easy (bookwork) |
| Divided difference error bound (with FP) | 5 | Medium (method seen) |
| Dual number variant | 8–12 | Medium-Hard (some unseen) |
| LU or Cholesky by hand | 6 | Easy (drill this) |
| Rotation/reflection construction | 3–5 | Medium |
| Householder/QR proof by induction | 4–6 | Hard (unseen variants) |
| Backward error analysis | 6 | Hard |
| DFT definition + unitarity proof | 4–6 | Easy-Medium (bookwork) |
| Orthogonal polynomial construction | 4–5 | Medium (Stieltjes) |
| Derivative = new family proof | 4–5 | Medium-Hard |
| Gauss quadrature computation | 4–5 | Medium |

### Common pitfalls from examiner comments

1. **Floating point bounds (2024 comments):** Many students made leaps to fit the given answer. The bounds must be *clearly justified* — don't skip steps when $\epsilon_m^2$ terms appear.

2. **Householder method (2024 comments):** Students were unsure what was expected. Know the formula for the reflection that maps $x$ to $\|x\|e_1$ (or $e_n$).

3. **Confusing $P$ and $P^{-1}$** (2023 comments): In $A = P^\top LU$, $P$ permutes *rows* of $A$.

4. **Taylor's theorem precision (2024 comments):** State it precisely with remainder term. Don't just wave hands.

---

## Key Problem Sheet Questions to Review

Questions marked *(optional)* on the course GitHub are not directly examined but may give good practice for "unseen" style questions.

### Sheet 1 — Rectangular Rule & Divided Differences
- **Q1:** Left-sided rectangular rule error bound — same proof technique as lectures, good warmup
- **Q2(a):** One-panel trapezium rule error bound ($O(h^3)$) — the integration-by-parts-twice proof technique is important
- **Q2(b):** Global trapezium rule error bound ($O(h^2)$) — summing panel errors, mirrors the lecture proof
- **Q3:** Left-sided divided difference error bound — almost identical to right-sided (lectures)
- **Q4:** Central differences convergence rate ($O(h^2)$) — **directly examined 2024 Q4b**. Uses third-order Taylor expansion
- **Q5:** Second derivative approximation error bound *(optional)*

### Sheet 2 — Dual Numbers & Newton's Method
- **Q1:** Computing polynomial derivatives via dual number arithmetic by hand — core skill, practise until automatic
- **Q2:** Dual extensions of $x^{100}+1$, $1/x$, $\tan x$ — knowing the extensions for standard functions
- **Q3:** Differentiating compositions with dual numbers (e.g. $\exp(\exp(x)\cos(x) + \sin(x))$) — *(Q3(b) optional)*
- **Q5(a):** 2D dual number product formula — **directly examined 2024 Q4c(i)**
- **Q5(b):** 2D dual numbers compute gradients of polynomials — **directly examined 2024 Q4c(ii)**
- **Q5(c):** Computing a specific gradient via 2D duals — good practice
- **Q6:** Newton convergence for degenerate case ($f'(r) = 0$) — linear convergence bound $|\epsilon_{k+1}| \leq \tilde{M}|\epsilon_k|$

### Sheet 3 — Floating Point Arithmetic
- **Q1:** $\pi$ in binary — basic skill
- **Q2:** $F_{32}$ representations for specific numbers — **examined every year Q4a(ii)**
- **Q3:** Float spacing ($m(y) - y$ for different $y$) — builds intuition for machine epsilon
- **Q4:** Exact floating point computations (when there's no rounding error)
- **Q5:** $1/5$ in binary and rounding in $F_{16}$ — good practice for the rounding rules (round-to-even)
- **Q6:** Bounding FP expressions: $(fl(1.1) \otimes fl(1.2)) \oplus fl(1.3)$ and $(fl(1.1) \ominus 1) \oslash fl(0.1)$ — **core exam technique**, directly mirrors lecture Example 9
- **Q7:** Central differences with floating point error bound — **directly examined 2023 Q4b, 2024 Q4b** *(optional but highly exam-relevant)*

### Sheet 4 — Interval Arithmetic & Structured Matrices
- **Q1:** Proving interval arithmetic formulae ($X/n$, $XY$) and generalising to negative values
- **Q2(a):** Interval arithmetic for computing $\sin 1$ with rigorous bounds *(optional)*
- **Q2(b):** Further interval computations — good pen-and-paper practice
- **Q2(c):** Extended interval computations *(optional)*
- **Q3:** Floating point error in bidiagonal matrix-vector multiplication — **directly examined 2024 Q5a**. Shows $\|\delta\|_\infty \leq \frac{3}{2}\epsilon_m \|A\|_\infty \|x\|_\infty$

### Sheet 5 — LU, PLU & Cholesky Factorisations
- **Q1:** Backward error for banded matrices *(optional)* — but the technique is **directly examined 2023 Q5c**
- **Q2:** LU factorisation by hand — **examined every year Q5a**. Drill until effortless
- **Q3:** PLU factorisation by hand — know the pivoting procedure cold
- **Q4:** Cholesky factorisation by hand — **examined every year** (2024 Q5b, 2025 Q5a). Free marks if practised
- **Q5:** Reverse Cholesky ($A = U U^\top$) — good extension, tests understanding of the proof structure

### Sheet 6 — Orthogonal Matrices & QR Factorisation
- **Q1:** Constructing simple rotations and reflections *(optional)* — but rotations **examined 2025 Q5b**
- **Q2:** Properties of orthogonal matrices: (i) norm-preserving, (ii) eigenvalues on unit circle, (iii) $\det = \pm 1$, (iv) normal, (v) $Q = I$ iff all eigenvalues are 1 — all important for proofs
- **Q3:** Householder reflection construction — **directly examined 2023 Q5b, 2024 Q5c(i), 2025 Q5c(i)**
- **Q4:** QR factorisation by hand via Householder *(optional — hard to get nice numbers)*
- **Q5:** QR uniqueness (when $R$ has positive diagonal)

### Sheet 7 — Interpolation, Regression & SVD
- **Q1(a):** Polynomial interpolation via Lagrange basis — good background
- **Q1(b):** Least squares via QR by hand — the full Householder→solve pipeline
- **Q2:** Vandermonde matrix properties and determinant — connects to interpolation theory
- **Q3:** Interpolatory quadrature rules — **directly examined 2025 Q6b(iii)** (Gauss quadrature is a special case)
- **Q4:** SVD and pseudoinverse *(optional)*
- **Q5:** Vandermonde determinant *(optional)*

### Sheet 8 — Fourier Series & DFT (from previous year, covered in project chat)
- **Q1:** Computing discrete Fourier coefficients via aliasing formula — **directly examined 2023 Q6a**
- **Q2:** Decay of Fourier coefficients via integration by parts — important for understanding convergence
- **Q3:** DFT matrix properties: (a) unitarity proof, (b) bounding sums by integrals, (c) interpolation — **(a) directly examined 2025 Q6a(ii)**
- **Q4:** Discrete Cosine Transform — extends the DFT theory
- **Q5:** Interpolation via DFT — **directly examined 2024 Q6a**

### Sheet 9 — Orthogonal Polynomials & Gauss Quadrature
- **Q1:** Constructing orthogonal polynomials via Gram-Schmidt/Stieltjes — **directly examined every year Q6b**
- **Q2:** Three-term recurrence and Jacobi matrix — key theoretical framework
- **Q3:** Symmetry of Jacobi matrix when weight is even ($a_k = 0$) — **examined 2023 Q6c(i)**
- **Q4:** Gauss quadrature nodes as eigenvalues of Jacobi matrix — **examined 2023 Q6c(ii), 2024 Q6b(ii-iii), 2025 Q6b(iii)**
- **Q5:** Derivatives of orthogonal polynomial families — **examined 2023 Q6b(ii) ($U'_{n+1} = 2C_n$), 2025 Q6b(ii) ($C_n = P'_{n+1}$)**

### Priority Ranking for Revision

If time is limited, focus on these in order:

1. **Sheet 5 Q2–Q4** (LU/PLU/Cholesky by hand) — guaranteed easy marks every year
2. **Sheet 3 Q6** (bounding FP expressions) — the core technique for Q4
3. **Sheet 6 Q3** (Householder reflections) — appears every year in Q5
4. **Sheet 9 Q1, Q4** (orthogonal polynomial construction + Gauss quadrature) — the hardest Q6 marks
5. **Sheet 2 Q1, Q5** (dual numbers + 2D variant) — the trickiest part of Q4
6. **Sheet 1 Q4** (central differences error) — often examined directly
7. **Sheet 8 Q1, Q3(a)** (Fourier coefficients + DFT unitarity) — bookwork marks in Q6