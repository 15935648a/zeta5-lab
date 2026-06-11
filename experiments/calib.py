"""
Calibration: determinant capacity of the Apery slit plane vs arithmetic cost.

Pure-stdlib (big ints + Fraction). No external packages.

Stage A: Apery sequences u_n, v_n; verify integrality of u, denominator
         structure 2*d_n^3*v_n in Z, and decay rate of l_n = u_n*zeta3 - v_n.
Stage B: Hankel determinants H_m of (l_n) in fixed point; measure
         a_m = log|H_m|/m^2  ->  expected -log(4*x2), x2=(1+sqrt2)^4.
Stage C: exact denominator D_m of the Hankel determinant viewed as a
         polynomial in t (t = placeholder for zeta3), via exact integer
         determinants at t=0..m and Lagrange interpolation.
         Compare log(D_m)/m^2 with the column-clearing bound log(S_m)/m^2.
"""
from fractions import Fraction
from math import lcm, log, sqrt

# ---------- Apery recurrence ----------
def apery(NMAX):
    u = [Fraction(1), Fraction(5)]
    v = [Fraction(0), Fraction(6)]
    for n in range(2, NMAX + 1):
        a = 34 * n**3 - 51 * n**2 + 27 * n - 5
        b = (n - 1) ** 3
        u.append((a * u[-1] - b * u[-2]) / n**3)
        v.append((a * v[-1] - b * v[-2]) / n**3)
    return u, v

NSEQ = 64
u, v = apery(NSEQ)
assert all(x.denominator == 1 for x in u), "u_n not integral?!"
U = [int(x) for x in u]

d = [1]
for k in range(1, NSEQ + 1):
    d.append(lcm(d[-1], k))

print("== Stage A: structure checks ==")
bad = [n for n in range(NSEQ + 1) if (2 * d[n] ** 3 * v[n]).denominator != 1]
print("2*d_n^3*v_n integral for all n<=%d : %s" % (NSEQ, not bad))
if bad:
    print("  failures at n =", bad[:10])

# high-precision zeta(3) from the convergent v_N/u_N (error ~ x2^(-2N))
NZ = 1500
uz, vz = apery(NZ)
zeta3 = vz[NZ] / uz[NZ]          # Fraction, error < 1e-4500
del uz, vz

P = 4500                          # fixed-point decimal digits
SC = 10 ** P
x2 = (1 + sqrt(2)) ** 4
log_x2 = 4 * log(1 + sqrt(2))
log_4x2 = log(4) + log_x2

# l_n in fixed point (exact rational -> floor scale)
Lfx = []
for n in range(NSEQ + 1):
    ln = u[n] * zeta3 - v[n]      # Fraction
    Lfx.append((ln.numerator * SC) // ln.denominator)

# decay-rate check
print("\n n   log(l_n)/n   (target -log x2 = %.6f)" % (-log_x2))
for n in (10, 20, 30, 40):
    ln_log = log(Lfx[n]) - P * log(10)
    print("%3d   %.6f" % (n, ln_log / n))
# second difference estimate of the rate (removes the n^{-3/2} prefactor)
r1 = (log(Lfx[40]) - log(Lfx[30])) / 10
print("slope n=30..40: %.6f" % r1)

# ---------- exact integer determinant (Bareiss) ----------
def det_int(M):
    M = [row[:] for row in M]
    n = len(M)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
            M[i][k] = 0
        prev = M[k][k]
    return sign * M[-1][-1]

# ---------- Stage B: Hankel determinants of l_n ----------
print("\n== Stage B: Hankel determinant capacity ==")
print("targets: classical -log x2 = %.5f ; slit-plane -log(4 x2) = %.5f"
      % (-log_x2, -log_4x2))
logH = {}
MMAX = 21
print(" m    log|H_m|/m^2    2nd-diff/2    H_m>0")
for m in range(1, MMAX + 1):
    M = [[Lfx[i + j] for j in range(m)] for i in range(m)]
    D = det_int(M)                       # = H_m * 10^(m*P), perturbed
    pos = D > 0
    lH = log(abs(D)) - m * P * log(10)
    logH[m] = lH
    sd = ""
    if m >= 3:
        sd = "%.5f" % ((logH[m] - 2 * logH[m - 1] + logH[m - 2]) / 2)
    print("%2d    %10.5f    %10s    %s" % (m, lH / m**2, sd, pos))

# ---------- Stage C: exact denominators of H_m(t) ----------
def polymul(a, b):
    r = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r

def lagrange(xs, ys):
    n = len(xs)
    out = [Fraction(0)] * n
    for i in range(n):
        num = [Fraction(1)]
        den = Fraction(1)
        for j in range(n):
            if j == i:
                continue
            num = polymul(num, [Fraction(-xs[j]), Fraction(1)])
            den *= xs[i] - xs[j]
        w = ys[i] / den
        for k in range(len(num)):
            out[k] += w * num[k]
    return out

print("\n== Stage C: exact denominator of H_m as polynomial in t ==")
print("              (cost = log(D_m)/m^2 ; clearing bound = log(S_m)/m^2 ; asympt 4.5)")
print(" m    cost(D_m)    bound(S_m)    margin = 4.9118 - cost")
for m in range(2, 14):
    s = [2 * d[m - 1 + j] ** 3 for j in range(m)]
    S = 1
    for x in s:
        S *= x
    ys = []
    xs = list(range(m + 1))
    for t in xs:
        Mt = [[int((Fraction(U[i + j] * t) - v[i + j]) * s[j]) for j in range(m)]
              for i in range(m)]
        ys.append(Fraction(det_int(Mt), S))
    coeffs = lagrange([Fraction(x) for x in xs], ys)
    # sanity: leading coeff = det of integer Hankel of U -> denominator 1
    assert coeffs[m].denominator == 1, "leading coeff not integral?!"
    Dm = 1
    for c in coeffs:
        Dm = lcm(Dm, c.denominator)
    cost = log(Dm) / m**2
    bnd = log(S) / m**2
    print("%2d    %8.5f     %8.5f     %8.5f" % (m, cost, bnd, log_4x2 - cost))

print("\nreference ledger:")
print("  classical route : need delta < log x2      : 3      vs %.5f  (margin %.5f)"
      % (log_x2, log_x2 - 3))
print("  naive Hankel    : need cost  < log(4 x2)   : 4.5    vs %.5f  (margin %.5f)"
      % (log_4x2, log_4x2 - 4.5))
