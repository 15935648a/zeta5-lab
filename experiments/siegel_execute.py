"""
Route 3: numerically EXECUTE the Siegel construction on Apery data and
measure whether the staircase denominator rate tau-flat(b) is realized
for theta-tower tuples.

Surrogate: V (Apery's rational sequence, exact type [1..n]^3) replaces the
hypothetical P. Tuples:
  T3 = {1, V, thetaV}          types (0,3,3)   staircase ratio 8/9  = 0.8889
  T4 = {1, V, thetaV, th2V}    types (0,3,3,3) staircase ratio 15/16 = 0.9375
  TA = {1, A, thetaA}          control, types (0,0,0): ratio ~ 0
For each degree D and slack s: find all (Q_1..Q_m), deg < D, with
F = sum Q_i f_i vanishing to order >= mD-1-s; search the kernel slice for
the vector minimizing den(beta), beta = leading coefficient of F.
Report  r := log den(beta) / (3*psi(n*))  where psi(n) = log lcm(1..n).
  crude rate  : r -> 1
  staircase   : r -> 1 - 1/m^2
"""
from fractions import Fraction
from math import lcm, log
from itertools import product as iproduct

def apery(NMAX):
    u = [Fraction(1), Fraction(5)]
    v = [Fraction(0), Fraction(6)]
    for n in range(2, NMAX + 1):
        a = 34 * n**3 - 51 * n**2 + 27 * n - 5
        b = (n - 1) ** 3
        u.append((a * u[-1] - b * u[-2]) / n**3)
        v.append((a * v[-1] - b * v[-2]) / n**3)
    return u, v

NSEQ = 75
uA, vV = apery(NSEQ)
dl = [1]
for k in range(1, NSEQ + 1):
    dl.append(lcm(dl[-1], k))
PSI = [0.0] + [log(x) for x in dl[1:]]

def coeffs(seq, dd, n):
    # coefficient of x^n in theta^dd applied to gen.func of seq
    return seq[n] * (n ** dd) if dd > 0 else seq[n]

def kernel_basis(M):
    """Exact kernel of rational matrix M (rows = constraints)."""
    rows = [r[:] for r in M]
    ncol = len(M[0])
    piv = {}
    rr = 0
    for c in range(ncol):
        sel = None
        for r in range(rr, len(rows)):
            if rows[r][c] != 0:
                sel = r
                break
        if sel is None:
            continue
        rows[rr], rows[sel] = rows[sel], rows[rr]
        pv = rows[rr][c]
        rows[rr] = [x / pv for x in rows[rr]]
        for r in range(len(rows)):
            if r != rr and rows[r][c] != 0:
                f = rows[r][c]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[rr])]
        piv[c] = rr
        rr += 1
        if rr == len(rows):
            break
    basis = []
    for c in range(ncol):
        if c in piv:
            continue
        w = [Fraction(0)] * ncol
        w[c] = Fraction(1)
        for pc, pr in piv.items():
            w[pc] = -rows[pr][c]
        basis.append(w)
    return basis

def to_primitive(w):
    den = 1
    for x in w:
        den = lcm(den, x.denominator)
    iv = [int(x * den) for x in w]
    from math import gcd
    g = 0
    for x in iv:
        g = gcd(g, abs(x))
    return [x // (g if g else 1) for x in iv]

def run_tuple(name, seqs, Ds, slacks, ratio_pred):
    m = len(seqs)
    print("\n-- %s (m=%d, staircase ratio %.4f) --" % (name, m, ratio_pred))
    print("   D   s   n*   log den(beta)   3*psi(n*)   ratio")
    for D in Ds:
        N = m * D
        for s in slacks:
            ntar = N - 1 - s
            if ntar + 12 > NSEQ:
                continue
            # constraint rows k=0..ntar-1; unknowns (i,j): Q_i coeff of x^j
            M = []
            for k in range(ntar):
                row = []
                for (seq, dd) in seqs:
                    for j in range(D):
                        row.append(coeffs(seq, dd, k - j) if k - j >= 0 else Fraction(0))
                M.append(row)
            basis = kernel_basis(M)
            if not basis:
                continue
            prim = [to_primitive(w) for w in basis]
            best = None
            rng = range(-2, 3) if len(prim) > 1 else range(1, 2)
            for cvec in iproduct(rng, repeat=len(prim)):
                if all(c == 0 for c in cvec):
                    continue
                w = [sum(c * b[t] for c, b in zip(cvec, prim)) for t in range(N)]
                if all(x == 0 for x in w):
                    continue
                # leading coefficient of F
                nstar, beta = None, None
                for k in range(ntar, ntar + 12):
                    val = Fraction(0)
                    t = 0
                    for (seq, dd) in seqs:
                        for j in range(D):
                            if w[t] and k - j >= 0:
                                val += w[t] * coeffs(seq, dd, k - j)
                            t += 1
                    if val != 0:
                        nstar, beta = k, val
                        break
                if beta is None:
                    continue
                ld = log(beta.denominator) if beta.denominator > 1 else 0.0
                key = (ld, nstar)
                if best is None or ld / max(1e-9, 3 * PSI[nstar]) < \
                   best[0] / max(1e-9, 3 * PSI[best[1]]):
                    best = (ld, nstar)
            if best:
                ld, nstar = best
                den3psi = 3 * PSI[nstar]
                print("  %2d  %2d  %3d   %10.3f   %9.3f    %.4f"
                      % (D, s, nstar, ld, den3psi, ld / den3psi if den3psi else 0))

V3 = [(vV, 0), (vV, 1)]
run_tuple("T3 = {1, V, thetaV}", [([Fraction(1)] + [Fraction(0)] * NSEQ, 0),
                                  (vV, 0), (vV, 1)],
          [4, 6, 8, 10, 12, 14, 16], [0, 2], 8.0 / 9)
run_tuple("T4 = {1, V, thetaV, theta2V}",
          [([Fraction(1)] + [Fraction(0)] * NSEQ, 0),
           (vV, 0), (vV, 1), (vV, 2)],
          [3, 5, 7, 9, 11, 13], [0, 2], 15.0 / 16)
run_tuple("TA control = {1, A, thetaA}",
          [([Fraction(1)] + [Fraction(0)] * NSEQ, 0),
           (uA, 0), (uA, 1)],
          [6, 10, 14], [0], 0.0)
