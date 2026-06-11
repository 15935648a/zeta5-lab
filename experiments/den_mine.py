"""
Denominator mine: quantify and LOCALIZE the extra savings in den(beta)
for Siegel constructions over the Apery theta-towers.

For each (tuple, D, slack): find min-den kernel vector, then factor
den(beta) prime-by-prime and compare with the cap 3*floor(log_p n*):
savings binned by prime size (p <= n/4, n/4..n/2, n/2..3n/4, 3n/4..n).
Output: tau_eff trend + implied improved thresholds for the stacked
criterion  lambda > tau_eff * m/(m-1).
"""
from fractions import Fraction
from math import lcm, log, exp
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

NSEQ = 90
uA, vV = apery(NSEQ)
dl = [1]
for k in range(1, NSEQ + 1):
    dl.append(lcm(dl[-1], k))
PSI = [0.0] + [log(x) for x in dl[1:]]

def primes_upto(n):
    s = list(range(n + 1))
    P = []
    for p in range(2, n + 1):
        if s[p] == p:
            P.append(p)
            for q in range(p * p, n + 1, p):
                if s[q] == q:
                    s[q] = p
    return P

def coeffs(seq, dd, n):
    return seq[n] * (n ** dd) if dd > 0 else seq[n]

def kernel_basis(M):
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
    from math import gcd
    den = 1
    for x in w:
        den = lcm(den, x.denominator)
    iv = [int(x * den) for x in w]
    g = 0
    for x in iv:
        g = gcd(g, abs(x))
    return [x // (g if g else 1) for x in iv]

def mine(name, seqs, Ds, slacks):
    m = len(seqs)
    print("\n== %s (m=%d; staircase ratio %.4f) ==" % (name, m, 1 - 1.0 / m**2))
    print("   D   s   n*    tau_eff   ratio    savings by prime bin (nats):")
    print("                                    [<=n/4] (n/4,n/2] (n/2,3n/4] (3n/4,n]")
    results = []
    for D in Ds:
        N = m * D
        for s in slacks:
            ntar = N - 1 - s
            if ntar + 12 > NSEQ:
                continue
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
            R = 2 if len(prim) <= 3 else 1
            best = None
            for cvec in iproduct(range(-R, R + 1), repeat=len(prim)):
                if all(c == 0 for c in cvec):
                    continue
                w = [sum(c * b[t] for c, b in zip(cvec, prim)) for t in range(N)]
                if all(x == 0 for x in w):
                    continue
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
                if best is None or ld / max(1, 3 * PSI[nstar]) < \
                   best[0] / max(1, 3 * PSI[best[1]]):
                    best = (ld, nstar, beta.denominator)
            if not best:
                continue
            ld, nstar, den = best
            # prime profile of savings
            bins = [0.0, 0.0, 0.0, 0.0]
            for p in primes_upto(nstar):
                cap = 0
                q = p
                while q <= nstar:
                    cap += 1
                    q *= p
                vp = 0
                dd2 = den
                while dd2 % p == 0:
                    vp += 1
                    dd2 //= p
                save = (3 * cap - vp) * log(p)
                idx = min(3, int(4 * p / (nstar + 1)))
                bins[idx] += save
            ratio = ld / (3 * PSI[nstar])
            results.append((nstar, ld))
            print("  %2d  %2d  %3d    %6.3f   %.4f    %7.1f %8.1f %9.1f %8.1f"
                  % (D, s, nstar, ld / nstar, ratio, *bins))
    if len(results) >= 4:
        # tail fit of tau_eff
        tail = sorted(results)[-4:]
        te = sum(ld / n for n, ld in tail) / len(tail)
        print("  tail tau_eff ~ %.3f  (theory cap tau-flat*sigma-normalized: %.3f, crude 3*psi(n)/n)"
              % (te, (1 - 1.0 / m**2) * 3))
        lam_req = te * m / (m - 1)
        print("  implied threshold IF provable: lambda > %.3f  <=>  x2 > %.2f (Koebe)"
              % (lam_req, exp(lam_req) / 4))

one = [Fraction(1)] + [Fraction(0)] * NSEQ
mine("T3 = {1, V, thetaV}", [(one, 0), (vV, 0), (vV, 1)],
     [6, 10, 14, 18, 22], [0, 2, 4])
mine("T4 = {1, V, thetaV, theta2V}", [(one, 0), (vV, 0), (vV, 1), (vV, 2)],
     [5, 8, 11, 14, 17], [0, 2, 4])
