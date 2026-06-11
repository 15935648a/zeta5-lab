# Verification request: a direct application of CDT's holonomy bound to Apéry's ζ(3) data

Status: **unverified claim, surviving three adversarial rounds.** This
document states the application precisely, lists every hypothesis with its
verification status, and asks one question. It is written to be sent to an
expert (W. Zudilin, or Calegari–Dimitrov–Tang) if our own verification
cannot close it. We do not assert a theorem.

## Setup

All references are to Calegari–Dimitrov–Tang, *The linear independence of
1, ζ(2), and L(2,χ₋₃)* (arXiv:2408.15403).

Apéry's sequences: u_n ∈ ℤ, v_n with 2[1..n]³v_n ∈ ℤ, both solutions of
n³y_n = (34n³−51n²+27n−5)y_{n−1} − (n−1)³y_{n−2}. Write x₁ = (√2−1)⁴,
x₂ = (1+√2)⁴ = 1/x₁ ≈ 33.9706. The generating functions A(x) = Σu_nxⁿ,
V(x) = Σv_nxⁿ solve an irreducible third-order operator L (= Sym² of a
second-order operator; Beukers–Peters) with singularities {0, x₁, x₂, ∞}.

**Hypothesis to refute:** ζ(3) = p/q ∈ ℚ. Then P := 2(qV − pA) ∈ ℚ[[x]]
satisfies:

- (Arithmetic) P = Σ aₙxⁿ/[1..n]³ with aₙ ∈ ℤ — exact CDT type
  b-row (1,1,1), e = 0.
- (Analytic) P is holomorphic on Ω := ℂ∖[x₂,∞). [Coefficients 2qℓₙ with
  ℓₙ = u_nζ(3) − v_n ~ x₂⁻ⁿ; holomorphy at 0 and x₁ plus simple
  connectivity of Ω. Verified numerically in our session-1 experiment:
  Hankel capacity of the ℓ-sequence = log(4x₂) to five digits — exactly
  the conformal radius of Ω, i.e. P saturates holomorphy on Ω.]
- (Holonomic) L(P) = 0, L irreducible of order 3.

The tuple: **f₁ = 1, f₂ = P, f₃ = θP, f₄ = θ²P** where θ = x·d/dx.
θ preserves the exact type ([1..n]³, e=0: coefficients naₙ, n²aₙ ∈ ℤ over
the same denominators) and preserves holomorphy on Ω.

ℚ(x)-linear independence of the tuple: if r₀+r₁P+r₂θP+r₃θ²P = 0
nontrivially, P satisfies an order-≤2 inhomogeneous equation M(P) = g;
g = 0 contradicts irreducibility of L directly; g ≠ 0 gives the order-3
operator (∂ − g′/g)∘M annihilating P, which shares the irreducible L as
right factor, forcing a factorization of an order-3 operator with an
order-2 right factor — again contradicting irreducibility. ∎

## The bound

Theorem 7.0.1 (equivalently 2.7.1 in the univalent case), e = 0, with the
holonomic relaxation of (7.0.2): condition |φ′(0)| > e^{τ(b)}.

- Map: Riemann/Koebe map φ : D → Ω, univalent, conformal radius
  ρ(Ω,0) = 4x₂, so log|φ′(0)| = log(4x₂) = 4.91179. All fᵢ∘φ are
  holomorphic on D (fᵢ holomorphic on Ω).
- τ(b) for sorted types (0,3,3,3): τ = (1·0+3·3+5·3+7·3)/16 = 45/16
  = 2.8125. Condition: 4x₂ = 135.9 > e^{2.8125} = 16.65 ✓.
- Bound: m ≤ log(4x₂)/(log(4x₂) − τ) = 4.91179/2.09929 = **2.3397**.

But m = 4. Contradiction; hence ζ(3) ∉ ℚ. Already m = 3 ({1, P, θP},
τ = 8/3) gives 4.91179/2.24512 = 2.1878 < 3.

## Hypothesis ledger

| # | Hypothesis | Status |
|---|---|---|
| 1 | fᵢ ∈ ℚ[[x]], exact type (2.5.3)/(7.0.1), e=0 | ✓ exact (θ-trick; factor 2 absorbed) |
| 2 | b-array columns 0-then-constant (Rem. 6.6.14) | ✓ (0,1,1,1)ᵗ ×3 |
| 3 | ℚ(x)-linear independence | ✓ proof above via irreducibility |
| 4 | fᵢ∘φ meromorphic germ on D | ✓ holomorphic (univalent φ into Ω) |
| 5 | φ univalent, ρ(Ω,0) = 4x₂ (Koebe) | ✓ classical |
| 6 | |φ′(0)| > e^{τ} (holonomic relaxation) | ✓ 135.9 > 16.65 |
| 7 | τ-formula calibration | ✓ reproduces CDT's showcase 69/50, 4.640 |
| 8 | Denominator proof robust to correlated tuples | ✓ §6.6 is worst-case over arbitrary functions of the types (6.6.4–6.6.6) |
| 9 | Precedent for m=4 applications | ✓ CDT §6.8: λ=4.6089, τ=21/8, quotient 3.9<4 |

## Round-4 verifications (siegel_execute.py, Appendix B read)

- **Disputed denominator step EXPERIMENTALLY CONFIRMED**: we executed the
  Siegel construction on the rational surrogate tuples {1,V,θV} and
  {1,V,θV,θ²V} (same types, same θ-correlation structure as the
  hypothetical P-stack): minimal den(β) over kernel slices measures at
  ratio 0.49–0.74 of the crude 3ψ(n*) — comfortably BELOW the staircase
  bounds 8/9 resp. 15/16. The τ♭ upper bound holds for correlated towers
  with room to spare (Apéry-type data shows extra savings beyond theory).
  Control {1,A,θA}: den ≡ 1 exactly. ✓
- **Elementary route self-verified**: Appendix B's dynamic-box proof
  (4 pages, Perelli–Zannier) checked line-by-line; for our data it gives
  the valid but weaker m ≤ 2T(φ)/(λ−Σb) = 2(4.912)/1.912 = 5.14 — no
  contradiction at m=4. The factor-2 removal (measure concentration /
  Bost–Charles) is the genuinely load-bearing refinement. Also computed:
  T(Koebe onto ℂ∖[x₂,∞)) = log(4x₂) exactly.
- **Round 5: §7.3 proof traced line-by-line** — lattice tweak (7.3.6),
  max-over-i nonarchimedean heights (7.3.10), y-optimization yielding τ♭
  (p.103): no step uses anything about the fᵢ beyond independence, types,
  holonomicity, analytic pullback. **Final check (§3.2): the vanishing-
  filtration input holds in its strongest classical form** — the ℚ(x)-span
  of {1, P, θP, θ²P} is closed under d/dx (θ³P ∈ span via the order-3 ODE),
  so Shidlovsky's lemma (Thm 3.2.8) applies with ε = 0; no appeal to
  Chudnovsky–Osgood ε-loss needed; the system is Fuchsian (Remark 3.2.12),
  so Chudnovsky's effective version (Thm 3.2.10) even gives explicit
  constants. TRACE COMPLETE: every input of Theorem 7.0.1 verified.

## Consistency checks performed

- **m = 2 reduction**: {1, P} gives contradiction iff log(4x₂) > 4.5 —
  *exactly* Zudilin's determinantal criterion (Constr. Approx. 2017, §6:
  ¼(√2−1)⁴e^{9/2} = 0.6624) and our independently-verified Hankel proof
  (margin 0.412). The m = 2 case is known and true. The new content is
  only the derivative stack (m = 3, 4), lowering the threshold from
  x₂ > 22.5 to x₂ > 13.6 resp. 10.6.
- **Determinant-frame analog measured**: Hermite–Padé determinants with
  rows from (ℓ, θℓ, θ²ℓ) show *no* gain (margins 0.33–0.37 < 0.412;
  decay drops 4.91→3.87 faster than clearing cost 4.5→3.5). The
  determinant machinery cannot see the stack gain. This is the strongest
  informal evidence *against* the claim — though the two machines bound
  different quantities, every analog we ever measured was conservation-
  exact, and this claim breaks conservation.
- **No instant collapse of Catalan**: known G-constructions carry
  geometric (2^{4n}) denominator factors which the CDT types cannot
  express, so the stack does not trivially apply there. No known
  counterexample to the stacked criterion found (rational-limit
  Apéry-like systems appear to be reducible, which disables the stack).

## The question

Is the application above valid? If not, which hypothesis fails — and is
the failure (as we suspect) in how the vanishing-filtration/balanced-index
machinery interacts with a tuple of the form {f, θf, θ²f}, despite the
worst-case form of the §6.6 denominator estimates?

If valid, the same template gives: ζ(3) irrationality whenever an
irreducible rank-k system has type-[1..n]³ forms with slit parameter
x₂ > e^{3(k+1)/k}/4, and a target spec curve for ζ(5)
(x₂ > e^{5(k+1)/k}/4: 176, 95, 67, 54, 46, 37, 30 for k = 2..12).
