"""
Round-2 empirical test: do DERIVATIVE ROWS beat the Hankel conservation law?

Sequences: l(d)_n = (n+1)...(n+d) * l_{n+d}  (coefficients of P^(d)), d=0,1,2.
Patterns (m x m, columns j=0..m-1):
  CONTROL : rows l(0)_{i+j}, i=0..m-1            (pure Hankel; margin -> 0.412)
  THIRDS  : d=0,1,2 each with shifts i=0..m/3-1  (kappa ~ 3.5)
  HALVES  : d=0,1   each with shifts i=0..m/2-1  (kappa ~ 3.75)
Measure rho_m = -log|D_m|/m^2 (2nd differences), exact kappa_m from actual
lcm clearing (each row (d,i): 2*d_{i+m-1+d}^3), margin = rho - kappa.

Conservation hypothesis : margins all -> 0.412.
Stack-legality hypothesis: margins grow as kappa drops.
"""
from fractions import Fraction
from math import lcm, log

def apery(NMAX):
    u = [Fraction(1), Fraction(5)]
    v = [Fraction(0), Fraction(6)]
    for n in range(2, NMAX + 1):
        a = 34 * n**3 - 51 * n**2 + 27 * n - 5
        b = (n - 1) ** 3
        u.append((a * u[-1] - b * u[-2]) / n**3)
        v.append((a * v[-1] - b * v[-2]) / n**3)
    return u, v

NSEQ = 50
u, v = apery(NSEQ)
d = [1]
for k in range(1, NSEQ + 1):
    d.append(lcm(d[-1], k))

NZ = 1500
uz, vz = apery(NZ)
zeta3 = vz[NZ] / uz[NZ]
del uz, vz

P = 3200
S = 10 ** P
Lfx = []
for n in range(NSEQ + 1):
    ln = u[n] * zeta3 - v[n]
    Lfx.append((ln.numerator * S) // ln.denominator)

def lder(dd, n):
    """l(d)_n * 10^P as exact int multiple of Lfx."""
    f = 1
    for t in range(1, dd + 1):
        f *= (n + t)
    return f * Lfx[n + dd]

def _flog(x):
    bl = x.bit_length()
    if bl <= 900:
        return log(x)
    sh = bl - 60
    return log(x >> sh) + sh * log(2)

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

def run(name, rowspec_fn, ms):
    print("\n-- %s --" % name)
    print("  m    rho_m     2nd-diff    kappa_m   margin_m")
    lh = {}
    prev = []
    for m in ms:
        rows = rowspec_fn(m)
        Mx = [[lder(dd, i + j) for j in range(m)] for (dd, i) in rows]
        D = det_int(Mx)
        if D == 0:
            print("  m=%d det=0!" % m)
            continue
        lg = _flog(abs(D)) - m * P * log(10)
        lh[m] = (lg, m)
        kap = sum(3 * _flog(d[i + m - 1 + dd]) + log(2) for (dd, i) in rows) / m**2
        sd = ""
        prev.append((m, lg))
        if len(prev) >= 3:
            (m0, l0), (m1, l1), (m2, l2) = prev[-3], prev[-2], prev[-1]
            h1, h2 = m1 - m0, m2 - m1
            if h1 == h2:
                sd = "%8.4f" % (-(l2 - 2 * l1 + l0) / (2 * h1 * h1))
        print("  %2d  %8.4f  %9s  %8.4f  %8.4f" % (m, -lg / m**2, sd, kap, -lg / m**2 - kap))

run("CONTROL (pure Hankel)", lambda m: [(0, i) for i in range(m)],
    [6, 9, 12, 15, 18, 21])
run("THIRDS (d=0,1,2; i<m/3)",
    lambda m: [(dd, i) for dd in (0, 1, 2) for i in range(m // 3)],
    [6, 9, 12, 15, 18, 21])
run("HALVES (d=0,1; i<m/2)",
    lambda m: [(dd, i) for dd in (0, 1) for i in range(m // 2)],
    [6, 10, 14, 18])

print("\nreference: conservation predicts margin -> 0.412 in every pattern;")
print("stack-legality predicts margins ~ (rho stays high while kappa drops).")
