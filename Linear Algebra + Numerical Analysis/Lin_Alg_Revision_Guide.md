# MATH50003 Linear Algebra — Exam Revision Guide

## How This Guide Is Organised

This guide covers the linear algebra topics most commonly examined in the MATH50003 final, based on analysis of the 2023, 2024, and 2025 papers. Topics are ordered roughly by how they build on each other. Each section flags whether definitions, proofs, or computations are most likely to be tested.

---

## 1. Similarity, Eigenvalues, and Diagonalisation

**Exam relevance:** Definitions appear as bookwork (category A marks). Diagonalisation criteria underpin almost every question.

### Key Definitions

**Similarity:** Matrices A, B are *similar* (written A ~ B) if there exists invertible P with B = P⁻¹AP. Similar matrices share: determinant, characteristic polynomial, eigenvalues, rank, trace, minimal polynomial.

**Algebraic multiplicity** a(λ): the power of (x − λ) dividing the characteristic polynomial.

**Geometric multiplicity** g(λ): dimension of the eigenspace E_λ = ker(T − λI).

**Key inequality:** g(λ) ≤ a(λ) always. (Proof uses block matrix structure — see notes §3.)

### Diagonalisation Criterion (Theorem 3.2)

T is diagonalisable ⟺ g(λᵢ) = a(λᵢ) for every eigenvalue λᵢ ⟺ Σ g(λᵢ) = dim V.

**Why this works:** If all geometric multiplicities equal algebraic multiplicities, then the union of bases of the eigenspaces gives a full basis of eigenvectors.

---

## 2. Direct Sums and T-Invariant Subspaces

**Exam relevance:** Definitions asked in 2023, 2024. The concept of T-invariant subspaces is tested every year — both in proofs and in concrete examples.

### Direct Sums

V = V₁ ⊕ ··· ⊕ Vₖ means every v ∈ V can be *uniquely* written as v = v₁ + ··· + vₖ with vᵢ ∈ Vᵢ.

**How to check (k = 2):** V = V₁ ⊕ V₂ ⟺ V₁ ∩ V₂ = {0} and dim V₁ + dim V₂ = dim V.

**General check:** V = V₁ ⊕ ··· ⊕ Vₖ ⟺ dim V = Σ dim Vᵢ and the union of bases of the Vᵢ forms a basis of V.

### T-Invariant Subspaces

W is *T-invariant* if T(W) ⊆ W. Key consequence: if V = V₁ ⊕ ··· ⊕ Vₖ with each Vᵢ T-invariant, then choosing a basis Bᵢ of each Vᵢ and setting B = B₁ ∪ ··· ∪ Bₖ gives

[T]_B = [T_{V₁}]_{B₁} ⊕ ··· ⊕ [T_{Vₖ}]_{Bₖ}

(a block-diagonal matrix). This is the fundamental mechanism behind *all* canonical form theorems.

### Quotient Spaces and the Quotient Map

If W is T-invariant with basis B_W extended to a basis B of V, and T̄ : V/W → V/W is defined by T̄(W + v) = W + T(v), then:

[T]_B = ( X  Z )
        ( 0  Y )

where X = [T_W]_{B_W} and Y = [T̄]_{B̄}. Crucially, c_T(x) = c_{T_W}(x) · c_{T̄}(x).

**Why this matters for exams:** The 2023 and 2024 papers ask you to determine which dimensions of T-invariant subspaces exist. Use: the min poly of T_W divides m_T(x), and c_{T_W}(x) divides c_T(x). If m_T(x) is irreducible, then V has no non-trivial T-invariant subspaces.

---

## 3. The Cayley-Hamilton Theorem

**Exam relevance:** Statement asked frequently (category A). The proof technique (using T-invariant subspaces and companion matrices) underpins the entire course.

### Statement

If T : V → V has characteristic polynomial p(x), then p(T) = 0.

### Proof Sketch

The proof is by induction on dim V.

**Case A — Non-trivial T-invariant subspace W exists:** Write [T]_B in block form. Since p(x) = p_X(x)·p_Y(x) where p_X, p_Y are the char polys of the diagonal blocks, and by induction p_X(X) = 0 and p_Y(Y) = 0, we get p(A) = 0 by block multiplication.

**Case B — No non-trivial T-invariant subspace:** For any v ≠ 0, the set {v, T(v), ..., T^{n-1}(v)} must be a basis of V (otherwise Sp(v, ..., T^{j-1}(v)) would be a proper T-invariant subspace). The matrix [T]_B is then a *companion matrix*, and the result follows from the definition.

---

## 4. Polynomials: Irreducibility and Factorisation

**Exam relevance:** Factorising polynomials over F₂ and F₃ comes up *every year*. Common source of errors in exams (markers repeatedly note students confusing "no roots" with "irreducible").

### Key Facts

- **Irreducible** means it cannot be written as a product of polynomials of smaller degree over that field.
- **Critical warning:** A polynomial with no roots is NOT necessarily irreducible! For example, x⁴ + x² + 1 = (x² + x + 1)² over F₂. You must check for factors of *all* degrees, not just linear ones.

### Irreducibles over F₂ of small degree

| Degree | Irreducible polynomials |
|--------|------------------------|
| 1      | x, x+1                |
| 2      | x²+x+1                |
| 3      | x³+x+1, x³+x²+1      |
| 4      | x⁴+x+1, x⁴+x³+1, x⁴+x³+x²+x+1 |

### Irreducibility over Q

For monic polynomials with integer coefficients: if α ∈ Q is a root, then α ∈ Z (and α divides the constant term). To show a polynomial is irreducible over Q, check it has no integer roots and (by Gauss's Lemma) cannot factor into lower-degree polynomials with integer coefficients.

### The Euclidean Algorithm and gcd

Given f, g ∈ F[x], repeatedly apply the division algorithm to get gcd(f,g). The key consequence is Bezout's identity: there exist r, s ∈ F[x] with gcd(f,g) = r·f + s·g. This is used in the Primary Decomposition proof.

---

## 5. The Minimal Polynomial

**Exam relevance:** Computing minimal polynomials is tested every year. The relationship between m_T(x) and c_T(x) is fundamental.

### Definition and Uniqueness

m_T(x) is the unique monic polynomial of smallest degree such that m_T(T) = 0.

**Key property:** p(T) = 0 ⟺ m_T(x) | p(x). (Proof: Euclidean algorithm.)

### Relationship with Characteristic Polynomial

1. m_T(x) | c_T(x) (by Cayley-Hamilton).
2. m_T(x) and c_T(x) have *exactly the same irreducible factors* (Theorem 9.3).
3. For a companion matrix C(p(x)): both the minimal and characteristic polynomials equal p(x).

### How to Compute m_T(x)

Given c_T(x) = ∏ fᵢ(x)^{aᵢ}, we know m_T(x) = ∏ fᵢ(x)^{nᵢ} with 1 ≤ nᵢ ≤ aᵢ. Test: compute fᵢ(T)^j for various j to find the smallest nᵢ such that the product annihilates V.

**Exam tip (from markers):** When A is a companion matrix, its minimal polynomial equals its characteristic polynomial. Recognising companion matrices saves significant computation time (noted as a common missed shortcut in 2024 markers' comments).

---

## 6. Primary Decomposition

**Exam relevance:** Used in nearly every JCF/RCF computation. The 2025 Q2(d) explicitly requires understanding the dimensions of primary components.

### Statement (Theorem 10.1)

Let m_T(x) = ∏ fᵢ(x)^{nᵢ} where f₁, ..., fₖ are distinct irreducible polynomials. Define Vᵢ = ker(fᵢ(T)^{nᵢ}). Then:

1. V = V₁ ⊕ ··· ⊕ Vₖ
2. Each Vᵢ is T-invariant
3. The restriction T_{Vᵢ} has minimal polynomial fᵢ(x)^{nᵢ}

### Proof Idea (for k = 2)

The core argument uses Bezout: if g₁, g₂ are coprime with g₁(T)g₂(T) = 0, then ∃ s₁, s₂ with s₁g₁ + s₂g₂ = 1, so s₁(T)g₁(T) + s₂(T)g₂(T) = I. This lets you decompose any v = v₁ + v₂ with vᵢ ∈ ker gᵢ(T), and show the sum is direct.

### Special Case: All Irreducibles Linear

When F = C (or whenever c_T(x) splits into linear factors), we get the *generalised eigenspace decomposition*:

Vᵢ = ker(T − λᵢI)^{nᵢ}

This is the starting point for Jordan Canonical Form.

### Key Application: Diagonalisation (Corollary 10.2)

T is diagonalisable ⟺ m_T(x) = ∏(x − λᵢ), a product of *distinct* linear factors.

**Intuition:** If all the powers nᵢ are 1, then each Vᵢ = ker(T − λᵢI) is just the eigenspace, and V is the direct sum of eigenspaces.

---

## 7. Jordan Canonical Form (JCF)

**Exam relevance:** JCF questions appear every year. Expect: definitions (A marks), computations (B marks), and harder proof/counting questions (C/D marks).

### Definitions

The **Jordan block** J_n(λ) is the n×n matrix with λ on the diagonal, 1's on the superdiagonal, and 0's elsewhere.

### Properties of Jordan Blocks (Proposition 11.1)

For J = J_n(λ):
- c_J(x) = m_J(x) = (x − λ)ⁿ
- λ is the only eigenvalue, with a(λ) = n and g(λ) = 1
- (J − λI)ⁱ has rank n − i (for i ≤ n) and sends eₙ → eₙ₋ᵢ

### The JCF Theorem (Theorem 11.3)

If A is n×n over F and c_A(x) is a product of linear factors, then A is similar to a *unique* (up to block ordering) block-diagonal matrix:

J = J_{n₁}(λ₁) ⊕ J_{n₂}(λ₂) ⊕ ··· ⊕ J_{nₖ}(λₖ)

### How to Compute the JCF

For each eigenvalue λ:

1. **Number of λ-blocks** = g(λ) = dim E_λ = dim ker(A − λI)
2. **Sum of λ-block sizes** = a(λ)
3. **Largest λ-block size** = highest power of (x−λ) dividing m_A(x)
4. **If still ambiguous:** compute rank(A − λI)², rank(A − λI)³, etc.

**The nullity method (from uniqueness proof):** Let nᵢ = null(A − λI)ⁱ. Then the number of blocks of size exactly j is:

aⱼ = (nⱼ − nⱼ₋₁) − (nⱼ₊₁ − nⱼ) = 2nⱼ − nⱼ₋₁ − nⱼ₊₁

(with appropriate boundary terms). This is the fully general method.

### Uniqueness Proof Sketch

For one eigenvalue λ, write J = J₁(λ)^{a₁} ⊕ ··· ⊕ J_r(λ)^{a_r}. The nullities nᵢ = null(A − λI)ⁱ satisfy a system of equations that uniquely determines all the aⱼ. For multiple eigenvalues, handle each separately since the λ-blocks of J and the blocks for other eigenvalues are independent.

### Finding a Jordan Basis (The Algorithm)

Given a nilpotent map S (i.e., S = T − λI for the one-eigenvalue case):

1. Compute the chain of kernels N₁ ⊂ N₂ ⊂ ··· ⊂ N_r = V where Nᵢ = ker(Sⁱ).
2. Choose a basis for V mod N_{r−1} — these are cyclic vectors for the largest blocks.
3. Apply S to get vectors in N_{r−1}, extend to a basis of N_{r−1} mod N_{r−2}, etc.
4. The Jordan basis consists of each cyclic vector v followed by S(v), S²(v), ..., reversed for each block.

**Exam example (2025 Q5a):** To find the LU factorisation, you're essentially doing Gaussian elimination. For the JCF: compute N₁ = ker(A − I), then work modulo N₁ to find cyclic vectors.

---

## 8. Companion Matrices and the Rational Canonical Form (RCF)

**Exam relevance:** RCF is tested every year. The 2024 and 2025 papers both ask for the RCF statement, companion matrix definition, and computations over F₂.

### Companion Matrix

For a monic polynomial p(x) = xⁿ + aₙ₋₁x^{n-1} + ··· + a₀:

C(p(x)) = ( 0  0  ···  0  −a₀  )
           ( 1  0  ···  0  −a₁  )
           ( 0  1  ···  0  −a₂  )
           ( ···                  )
           ( 0  0  ···  1  −aₙ₋₁ )

Both the characteristic and minimal polynomials of C(p(x)) equal p(x).

### Cyclic Subspaces

For v ≠ 0, the cyclic subspace Z(v, T) = Sp(v, T(v), T²(v), ...) has:
- Basis {v, T(v), ..., T^{k-1}(v)} where k = dim Z(v, T)
- The restriction T_v has matrix [T_v] = C(m_v), the companion matrix of the T-annihilator of v

### The RCF Theorem (Theorem 12.5)

Let m_T(x) = ∏ fᵢ(x)^{kᵢ} with fᵢ distinct irreducibles. Then T is similar to:

C(f₁^{k₁₁}) ⊕ ··· ⊕ C(f₁^{k₁r₁}) ⊕ ··· ⊕ C(fₜ^{kₜ₁}) ⊕ ··· ⊕ C(fₜ^{kₜrₜ})

where kᵢ = kᵢ₁ ≥ kᵢ₂ ≥ ··· ≥ kᵢrᵢ, all uniquely determined.

### Computing the RCF: The Rank Method

For each irreducible factor fᵢ of degree d, the RCF has blocks C(fᵢ), C(fᵢ²), ..., C(fᵢ^{kᵢ}). Call the numbers of each size a₁, a₂, ..., aₖ. Then:

- rank(fᵢ(T)^{k−1}) restricted to the fᵢ-primary component = aₖ · d
- rank(fᵢ(T)^{k−2}) restricted = aₖ₋₁ · d + 2aₖ · d

and so on. So from the ranks, you can read off all the aⱼ.

**Crucial detail (from 2025 Q2d):** When computing rank(fᵢ(A)) on the primary component for fᵢ, remember that fᵢ(T) acts invertibly on the other primary components (by Proposition 12.7(iii)). So:

rank(fᵢ(A)) = rank on fᵢ-component + dim(other components)

### Proving c_{C(f)} = f (2025 Q2b)

Expand det(xI − C(f)) along the first row. By induction, the (1,1)-minor gives g(x) = x^{n-1} + aₙ₋₁x^{n-2} + ··· + a₁, and the expansion yields c(x) = xg(x) + (−1)^{n-1}a₀(−1)^{n-1} = f(x).

---

## 9. Counting Similarity Classes

**Exam relevance:** Asked every year as a harder question (C/D marks). The 2023, 2024, and 2025 papers all have counting questions.

### General Strategy

1. **Identify possible characteristic polynomials** (products of irreducibles with correct total degree).
2. **For each char poly, determine possible minimal polynomials** (same irreducible factors, powers between 1 and the char poly power).
3. **For each min poly, count compatible RCFs/JCFs** using the rank constraints.

### Example: Counting with Constraints (as in 2025 Q1d)

If A³ = A, then m_A(x) | x³ − x. Over F₂: x³ − x = x(x+1)² (since x³ + x = x(x² + 1) = x(x+1)² in F₂). So the JCF is J₁(0)^a ⊕ J₂(1)^b ⊕ J₁(1)^c with a + 2b + c = n. Count all valid (a,b,c).

### Example: Matrices with No Eigenvalues over F₂ (2023 Q2d)

A 4×4 matrix with no eigenvalue in F₂ has RCF either C(f) for an irreducible degree-4 polynomial, or C(f₁) ⊕ C(f₂) for irreducible quadratics. Over F₂ there is only one irreducible quadratic (x²+x+1) and three irreducible quartics. Total: 4 similarity classes.

---

## 10. Inner Product Spaces

**Exam relevance:** Definitions asked every year (A marks). The 2025 Q3 is entirely about inner products.

### Definition

An inner product on V (over R or C) satisfies:
1. (λ₁v₁ + λ₂v₂, w) = λ₁(v₁, w) + λ₂(v₂, w) [left-linearity]
2. (w, v) = conjugate of (v, w) [conjugate symmetry]
3. (v, v) > 0 for v ≠ 0 [positive definiteness]

**Important subtlety:** Over C, the inner product is *not* right-linear: (v, λw) = λ̄(v,w).

### Verifying Inner Products (Common Exam Question)

To verify (A, B) = trace(AᵀB) is an inner product on Mₙ(R):
- **Left-linearity:** trace is linear, transpose is linear — follows immediately.
- **Symmetry:** (B, A) = trace(BᵀA) = trace((AᵀB)ᵀ) = trace(AᵀB) = (A, B).
- **Positive definiteness:** (A, A) = trace(AᵀA) = Σᵢ,ⱼ aᵢⱼ², which is > 0 when A ≠ 0.

(This specific example appeared in 2025 Q3c.)

### Orthogonality and Gram-Schmidt

V = W ⊕ W⊥ for any subspace W of an inner product space. The Gram-Schmidt process converts any basis to an orthonormal one.

### Orthogonal Projection

π_W(v) is the closest point in W to v. If {v₁,...,vᵣ} is an orthonormal basis of W:

π_W(v) = Σ (v, vᵢ) vᵢ

---

## 11. The Adjoint and the Spectral Theorem

**Exam relevance:** Definitions asked every year. The Spectral Theorem is stated in 2023, 2024, and 2025. Proofs involving self-adjoint maps are common (B/D marks).

### The Adjoint Map

T* : V → V is the unique linear map satisfying (T(u), v) = (u, T*(v)) for all u, v.

**Matrix version:** If E is an orthonormal basis and A = [T]_E, then [T*]_E = Āᵀ (conjugate transpose). So T is self-adjoint iff A = Āᵀ (symmetric if real, Hermitian if complex).

### Key Lemma (15.5) — Needed for Spectral Theorem Proof

If T is self-adjoint:
1. **Eigenvalues are real.** Proof: λ(v,v) = (Tv,v) = (v,Tv) = λ̄(v,v), so λ = λ̄.
2. **Eigenvectors for distinct eigenvalues are orthogonal.** Proof: λ(u,v) = (Tu,v) = (u,Tv) = μ(u,v), so (λ−μ)(u,v) = 0.
3. **If W is T-invariant, so is W⊥.** Proof: for x ∈ W⊥ and w ∈ W: (w, Tx) = (Tw, x) = 0 since Tw ∈ W.

### The Spectral Theorem (Theorem 15.3)

If T : V → V is self-adjoint, then V has an **orthonormal basis of T-eigenvectors**.

**Proof by induction:** T has a real eigenvalue λ (using that eigenvalues are real). Let u₁ be a unit eigenvector, W = Sp(u₁). Then W⊥ is T-invariant (Lemma 15.5(3)), and T restricted to W⊥ is still self-adjoint. By induction, W⊥ has an orthonormal basis of eigenvectors. Combining with u₁ gives the full basis.

**Matrix corollary:** Every real symmetric matrix A can be written as A = PDP⁻¹ where D is diagonal and P is orthogonal. Every complex Hermitian matrix is unitarily diagonalisable.

---

## 12. Bilinear and Quadratic Forms

**Exam relevance:** The 2024 paper has a full question on quadratic forms over Q. Definitions are commonly asked.

### Bilinear Forms

A bilinear form on V is a map V × V → F that is linear in both arguments. It is *symmetric* if (v,u) = (u,v), and *skew-symmetric* if (v,u) = −(u,v).

**Non-degenerate** means: (u,v) = 0 for all v ⟹ u = 0.

**Matrix of a bilinear form:** If B = {v₁,...,vₙ} is a basis, then f_B = (aᵢⱼ) where aᵢⱼ = (vᵢ, vⱼ). Under change of basis with matrix P: f_{B₂} = PᵀfB₁P (**congruence**, not similarity).

### Quadratic Forms

Q(v) = (v, v) where (,) is a symmetric bilinear form. Recover the bilinear form via:

(x, y) = ½(Q(x+y) − Q(x) − Q(y))

(requires char(F) ≠ 2).

### Diagonalisation of Symmetric Forms (Theorem 16.6)

Every non-degenerate symmetric bilinear form over a field with char ≠ 2 has an **orthogonal basis** (i.e., a basis where the matrix is diagonal).

**Proof:** Find v₁ with (v₁,v₁) ≠ 0 (exists since form is non-degenerate and char ≠ 2). Then V = Sp(v₁) ⊕ v₁⊥, and recurse on v₁⊥.

### Equivalence of Quadratic Forms

Q and Q' are equivalent if Q(x) = Q'(Px) for some invertible P. Equivalently, their matrices are congruent.

### Classification over R (Sylvester's Law of Inertia)

Every real non-degenerate quadratic form is equivalent to x₁² + ··· + x_p² − x_{p+1}² − ··· − x_n² for unique p, q with p + q = n. The pair (p, q) is called the **signature**.

### Classification over Q

Over Q, the classification is much harder. The 2024 Q3(c) asks whether Q(x) = x₁² + 2x₁x₂ − x₂² = k has solutions in Q² for various k. Strategy:
1. Diagonalise: complete the square to get Q' ~ λ₁x₁² + λ₂x₂².
2. For Q'(x) = 0: check whether −λ₂/λ₁ is a square in Q.
3. For Q'(x) = k: try to find an explicit solution, or use modular arithmetic to show no solution exists (e.g., consider the equation mod 3 or mod 8).

---

## 13. Unitary and Orthogonal Matrices

**Exam relevance:** Connected to the Spectral Theorem. The 2024 paper asks about expressing orthogonal matrices as products of reflections.

### Definitions

- **Orthogonal:** PᵀP = I (real matrices preserving the dot product)
- **Unitary:** P̄ᵀP = I (complex matrices preserving the inner product)

Both have |det P| = 1, and all eigenvalues have absolute value 1.

**Key fact:** The change-of-basis matrix between two orthonormal bases is unitary/orthogonal.

---

## Summary of What Gets Examined and How

### Every Year (expect these):
- **Define:** Jordan block, JCF theorem statement, companion matrix, RCF theorem statement, inner product, adjoint, Spectral Theorem statement, T-invariant subspace
- **Compute:** JCF of a given matrix, RCF of a matrix over F₂ or F₃, eigenvalues/minimal polynomial, verify an inner product
- **Prove:** Short arguments using Cayley-Hamilton, primary decomposition, spectral theorem properties

### Common Harder Questions (C/D marks):
- Count similarity classes given constraints (every year)
- Determine existence of T-invariant subspaces of specific dimensions
- Induction proofs (e.g., LQ factorisation in 2025, QL in 2023)
- Orthonormal eigenvector bases for specific operators
- Quadratic form problems over Q (2024)

### Common Pitfalls (from markers' comments):
- **2025:** Confusing minimal polynomial with characteristic polynomial in RCF context; not stating constraints on exponents
- **2024:** Not recognising companion matrices; failing to factorise polynomials properly over F₂
- **2023:** Assuming "no roots" implies irreducible (x⁴+x²+1 over F₂ is NOT irreducible!)
- **All years:** Sloppy proofs — state results from lectures explicitly when using them

---

## Key Problem Sheet Questions to Revisit


The following problem sheet questions cover techniques that appear directly in exam questions:

### Foundations (Sheets 1, 3)
- **Sheet 1 Q7:** Proving the characteristic polynomial of a companion matrix equals p(x) — this proof was asked in 2025 Q2(b).
- **Sheet 3 Q4:** Direct proof of Cayley-Hamilton for upper triangular matrices.
- **Sheet 3 Q8:** Finding irreducible polynomials over F₂ (degree 4) — essential for RCF questions.
- **Sheet 3 Q9:** Testing irreducibility over Q.

### Minimal Polynomial and Primary Decomposition (Sheet 4)
- **Sheet 4 Q1(b):** Proving m_{A₁⊕···⊕Aₖ} = lcm(m_{A₁},...,m_{Aₖ}) — used constantly.
- **Sheet 4 Q3(a):** Proving the minimal polynomial of C(p(x)) is p(x) — asked/used every year.
- **Sheet 4 Q5:** Computing a Primary Decomposition explicitly.

### JCF Computations and Theory (Sheets 5, 6)
- **Sheet 5 Q1:** Listing all possible JCFs for given char/min polys — directly exam-style.
- **Sheet 5 Q3:** Determining similarity of upper triangular matrices by computing JCFs.
- **Sheet 5 Q4:** Finding JCF from rank conditions — this is exactly the exam technique.
- **Sheet 5 Q6:** Proving J_n(λ) is similar to its transpose.
- **Sheet 6 Q2:** Finding JCF *and* the change-of-basis matrix P (Jordan basis computation).

### RCF and Counting (Sheet 7)
- **Sheet 7 Q1:** Listing all RCFs with given min poly — directly exam-style.
- **Sheet 7 Q4:** Proving C(fg) ~ C(f) ⊕ C(g) for coprime f, g.
- **Sheet 7 Q6:** Counting conjugacy classes in GL(3, F₃) and GL(4, F₂) — exam-style.

### Inner Products and Spectral Theorem (Sheets 8, 9)
- **Sheet 8 Q5:** Proving Cauchy-Schwarz consequences, Pythagoras, orthogonal set independence.
- **Sheet 8 Q8:** Proving existence of an upper triangular matrix with respect to an ONB (complex case).
- **Sheet 9 Q1:** Orthogonal projection computations — concretely useful.
- **Sheet 9 Q2:** Properties of the adjoint (especially ker(T*) = Im(T)⊥).
- **Sheet 9 Q5:** Positive self-adjoint maps and square roots — a style of question seen in exams.

### Bilinear/Quadratic Forms (Sheet 10)
- **Sheet 10 Q5:** Working with symmetric bilinear forms over Fₚ and finding orthogonal bases.
- **Sheet 10 Q6:** Diagonalising quadratic forms over Q and solving Q(x) = k — directly like 2024 Q3.
- **Sheet 10 Q7:** Proving surjectivity of non-degenerate quadratic forms with an isotropic vector.
