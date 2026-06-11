# zeta5-lab — Unipotent Holonomy Bound (UHB) Program

Goal: a dimension bound for arithmetically-constrained sections of *unipotent
extensions* on modular curves, strong enough to feed the conditional
construction below and contradict ζ(5) ∈ ℚ (intermediate target: Catalan G,
and linear independence of 1, ζ(3), ζ(5)).

## Why this route (established by the 2026-06 determinant campaign)

The scalar Hankel/capacity framework is *rigid*: four levers measured, each
closed by an empirical law (see experiments/ and memory notes):
capacity saturates; refinement bonus = log(4/#slits), exactly conserved;
pure ℓ-Hankel is the optimal determinant shape (integer blocks cannot decay,
touching supports give no Angelesco rescue); form-level purification costs
~2.5–3 per eliminated constant. Standing wins: second proof of ζ(3)
(margin 0.412), ζ(2) margin 0.792. G deficits −6.3 / −8.5.
Conclusion: progress requires new *arithmetic objects* (lower effective δ,
richer function supply), i.e. the modular/unipotent world — not new shapes.

## Objects

For even weight k, level 1 (later level N): the Eichler integral
  F_k(τ) = Σ_{n≥1} σ_{k−1}(n) n^{1−k} q^n = Σ_m m^{1−k} q^m/(1−q^m)
is a holomorphic section of a 2-step unipotent extension: it transforms with
weight 2−k up to a polynomial cocycle p_γ(τ) of degree ≤ k−2 whose
coefficients are the critical values of Λ(s) = (2π)^{−s}Γ(s)ζ(s)ζ(s+1−k).

**Load-bearing arithmetic fact (Stone B, verified numerically here):**
for k = 6 the cocycle coefficients are
  { ±ζ(5)/2 (slots y⁰, y⁴), ±π⁵/540 (slots y¹, y³), 0 (slot y²), ζ(6)/2π (y⁻¹) }.
The y² slot — where ζ(3)π²-type contamination would sit — is killed exactly
by the trivial zero ζ(−2) = 0. This is the modular incarnation of the purity
of Ext¹(Q(0), Q(5)): the mixing wall that defeats every integral construction
is *absent* in this packaging. Contrast k = 4: the y²-end slots carry ζ(3)/2.

## Conditional construction scheme

Assume ζ(5) = p/q. Then the k=6 cocycle becomes rational-up-to-π-powers, and
G₀ := corrected F₆ acquires (i) coefficients σ₅(n)/n⁵ with cleared denominator
growth δ ≤ 5 (adelically possibly less — measure), (ii) genuine quasi-modular
transformation with *rational* defect. The supply of derived objects:
Hecke translates T_p G₀ (same δ, p+1 per prime), multiplications by the
modular forms ring, level raisings. A UHB must bound such families by
(δ_eff + archimedean growth)/λ(domain) and lose to the count.

## Milestones

- **M1 (done 2026-06-11)**: cocycle data verified by self-derived Mellin
  computation + numerical fit (experiments/eisenstein_cocycle.py);
  purity slot confirmed empty at 1e−12.
- **M2 (RESOLVED 2026-06-11, from CDT §2)**: the λ-side has a HARD CEILING.
  CDT's bound is m ≤ λ/(λ−τ) with λ = log ρ(Ω,0) or log|φ'(0)|; Carathéodory
  caps |φ'(0)| ≤ 16 for maps to ℂ∖{1}, attained by the modular λ-map; more
  punctures (higher level) only lower it. Ceiling λ ≤ log 16 = 2.7726.
  Ledger: τ(log) = 1 ✓, τ(ζ(2)/L(2,χ)) = 2 ✓ barely (CDT's frontier),
  τ(ζ(3)) = 3 ✗, τ(ζ(5)) = 5 ✗. Direct UHB for ζ(5) is structurally dead;
  the analytic side cannot be improved. Also: τ(b) staircase = our 1.5δ
  formalized ((2i−1)/m² weights; cf. their (2.7.9) constant 2/3).
- **M3+ CANDIDATE RESULT (2026-06-11, cdt_direct_zeta3.py — UNVERIFIED,
  adversarial check required)**: applying CDT Thm 2.5.1 literally with the
  x₂-scaled bivalent map (2.11.2), Σ⁰={0,x₁}, Σ¹={x₂,∞}, to the tuple
  {1, P, P'} (P = Apéry's conditional function, type [1..n]³, slit-plane
  holomorphic — verified session 1): BC = log(8x₂)+4G/π = 6.771,
  log φ'(0) = 5.605, τ = 8/3 ⇒ bound 2.30 < 3 ⇒ contradiction ⇒ ζ(3) ∉ ℚ
  in three lines. m=2 lands knife-edge (2.018 vs 2). TOO STRONG TO TRUST:
  CDT know Apéry intimately and don't state this. Flaw-hunt list:
  (a) derivative-stacking vs exact (2.5.3) types (prime-power shifts);
  (b) hidden closed-disc/boundary hypotheses on φ (their showcase uses the
  same pole-on-boundary map, so unlikely); (c) single-valuedness of P∘φ
  (argued via monodromy factoring through π₁(ℂ∖{x₂}) since P holomorphic
  at 0, x₁ — spelled out, looks airtight); (d) ℚ(x)-independence of
  {1,P,P'} (only needs P non-first-order — trivial); (e) misreading of
  τ/BC (calibrated against their showcase 69/50 and 4.640 ✓).
  If valid: ζ(5) target spec curve computed (rank k, decay x₂):
  k=4 ⇒ x₂ ≥ 67.5; k=12 ⇒ x₂ ≥ 30; mixing with ζ(3) handled by
  independence framing.
  ADVERSARIAL ROUND 1 (done): m=2 univalent (Koebe onto ℂ∖[x₂,∞),
  λ = log 4x₂ = 4.912) IS Zudilin's criterion IS our Hankel proof —
  contradiction ⟺ log(4x₂) > 4.5, margin 0.412; triangle closed, digits
  match; the bivalent knife-edge was the same fact through a weaker
  instrument. The ONLY surviving new claim is the derivative stack
  (m=3: x₂ > 13.6; m=4: x₂ > 10.6). Legality evidence FOR: Remark 2.7.9's
  explicit m/(m+1) improvement; AGAINST: CDT silence on the implied ζ(3)
  reproof. ROUND 2 (next): (a) read §5/§6.6/§7 denominator-arithmetic
  proofs; (b) empirical: Hermite–Padé determinants with derivative rows
  (ℓ, ℓ', ℓ'') on Apéry data — does the measured margin exceed 0.412?
- **M3 (refined 2026-06-11, staircase_ledger.py)**: the naive staircase
  loophole (2 type-0 companions + P ⇒ need φ'(0) > 12.18 < 16, ladder
  reaches 13.45) SELF-DESTRUCTS: Borel–Pólya (CDT 2.7.10) forbids
  transcendental type-0 functions on λ>0 domains; companions carry their
  own floors (≥ (2/3)λ by 2.7.9), and re-inserting floors kills all small-m
  tuples. Structural lesson: the bound family is self-consistent —
  contradictions live ONLY in the gap between analytically-forced type
  FLOORS and hypothetical types. ζ(3) question becomes: is the floor for
  the weight-4 conditional function P (4-point configuration, covered
  δ = (√2−1)⁸ ≈ 0.000867) above 3? Simple staircase floor: 1.85 (far);
  CDT's Theorem A needed §6–8 fine bounds + §5 integrations to lift floors
  past the line at τ=2. NEXT: extract §6–8 fine-bound constants, compute
  the ζ(3)-configuration floor, measure the deficit to 3.
- **M4**: toy UHB: rank-2 (one extension step), lowest weight, prove the
  dimension bound rigorously in a model case; calibrate on a KNOWN statement
  (re-derive irrationality of ζ(3) or log 2 in unipotent packaging).
- **M5**: assemble for weight 6 / target inequality; iterate levels.

## Risks (honest)

- M2 may show λ structurally < 5 at every level → program dies (cheap, good).
- M4 is genuinely new mathematics; the archimedean estimates for sections
  with polynomial defects may erode λ (the defect grows at cusps).
- Known-territory caveat: CDT's arithmetic holonomy bounds and their
  L(2,χ₋₃) work are adjacent; literature check deferred deliberately
  (self-derivation first), required before claiming novelty.

## Method rules (from campaign lessons)

1. Measure the would-be theorem's constants before proving it.
2. Integer-clear all elimination coefficients (solver-denominator trap).
3. Watch for exact collapses under row/functional operations (overlap trap).
4. Every claim gets a numerical verification at ≥10 significant digits.
