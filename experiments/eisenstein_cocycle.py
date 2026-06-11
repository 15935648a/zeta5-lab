"""
Stone B / M1: verify the Eisenstein Eichler cocycle data for weights 6 and 4.

F_k(iy) = sum_{m>=1} m^{1-k} e^{-2pi m y}/(1 - e^{-2pi m y})

Predicted (self-derived via Mellin contour shift, residues of
Lambda(s) = (2pi)^{-s} Gamma(s) zeta(s) zeta(s+1-k)):

k=6:  D6(y) := F6(iy) - y^4 F6(i/y)
            =  (zeta(6)/2pi) y^{-1} - zeta(5)/2 + (pi^5/540) y
               + 0*y^2 - (pi^5/540) y^3 + (zeta(5)/2) y^4
      The y^2 slot (the zeta(3)pi^2 slot) is killed by zeta(-2)=0:  PURITY.

k=4:  D4(y) := F4(iy) + y^2 F4(i/y)
            =  (zeta(4)/2pi) y^{-1} - zeta(3)/2 + (pi^3/36) y - (zeta(3)/2) y^2
      Here zeta(3) OCCUPIES its end slots (contrast).

Fit numerically at sample points, compare coefficient ratios to predictions.
"""
from math import exp, pi, log

def zeta(s, N=60):
    # Euler-Maclaurin: sum_{n<N} n^-s + N^{1-s}/(s-1) + N^-s/2 + s*N^{-s-1}/12
    tot = sum(n ** -float(s) for n in range(1, N))
    tot += N ** (1 - s) / (s - 1) + 0.5 * N ** -float(s) + s * N ** (-s - 1) / 12.0
    tot -= s * (s + 1) * (s + 2) * N ** (-s - 3) / 720.0
    return tot

Z3, Z5 = zeta(3), zeta(5)
Z4, Z6 = pi ** 4 / 90, pi ** 6 / 945

def F(k, y, M=80):
    tot = 0.0
    for m in range(1, M + 1):
        e = exp(-2 * pi * m * y)
        if e == 0.0:
            break
        tot += m ** (1 - k) * e / (1 - e)
    return tot

def solve(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [a - f * bb for a, bb in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]

def fit(k, powers, Dfun, pts):
    A = [[y ** p for p in powers] for y in pts]
    b = [Dfun(y) for y in pts]
    return solve(A, b)

def run(k, powers, pred, pts, checks):
    best = None
    for sgn in (-1, +1):
        Df = lambda y: F(k, y) + sgn * y ** (k - 2) * F(k, 1 / y) * (1j ** (k - 2)).real
        Df = lambda y, s=sgn: F(k, y) + s * y ** (k - 2) * F(k, 1 / y)
        c = fit(k, powers, Df, pts)
        resid = max(abs(Df(y) - sum(cc * y ** p for cc, p in zip(c, powers)))
                    for y in checks)
        if best is None or resid < best[0]:
            best = (resid, sgn, c, Df)
    resid, sgn, c, Df = best
    print(" reflection sign: %+d   off-sample residual: %.2e" % (sgn, resid))
    print(" power   fitted            predicted         ratio")
    for p, cc, pr in zip(powers, c, pred):
        r = "%+.12f" % (cc / pr) if pr != 0 else ("ZERO-SLOT |c|=%.2e" % abs(cc))
        print(" y^%+d   %+.12f   %+.12f   %s" % (p, cc, pr, r))

print("== weight 6 ==")
run(6, [-1, 0, 1, 2, 3, 4, 5],
    [pi ** 5 / 1890, -Z5 / 2, pi ** 5 / 540, 0.0, -pi ** 5 / 540, Z5 / 2,
     -pi ** 5 / 1890],
    [0.62, 0.75, 0.88, 1.0, 1.14, 1.33, 1.58], (0.7, 0.95, 1.45))

print("\n== weight 4 (contrast: zeta(3) occupies its slots) ==")
run(4, [-1, 0, 1, 2, 3],
    [pi ** 3 / 180, -Z3 / 2, pi ** 3 / 36, -Z3 / 2, pi ** 3 / 180],
    [0.66, 0.82, 1.0, 1.21, 1.48], (0.74, 0.94, 1.35))

print("\nkey claims: (1) k=6 y^2 slot empty (purity of Ext^1(Q(0),Q(5)));")
print("            (2) k=6 ends carry zeta(5)/2; (3) k=4 ends carry zeta(3)/2.")
