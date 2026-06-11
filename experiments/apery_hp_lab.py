"""
Mixed Hermite-Pade determinant laboratory on Apery's zeta(3) data.

Rows drawn from two pools, columns j = 0..m-1:
  u-rows (shift i): entries u_{i+j}   -- integers, no arithmetic cost
  l-rows (shift i): entries l_{i+j}   -- tiny, each row costs ~ d_{i+m-1}^3

Omega_m = det(mixed matrix). Expanding all l-rows (l = u*z3 - v) makes
Omega a degree-k polynomial form in zeta(3) with rational coefficients;
rationality of zeta(3)=p/q is contradicted iff
   rho := -lim log|Omega|/m^2   >   kappa := lim log(den)/m^2
(q^k contributes only O(m), ignorable).

Shapes (m even, k = #l-rows):
  S1 control  : alpha=1   (pure l-Hankel; must reproduce 4.912 / 4.5)
  S2 alpha=3/4, S3 alpha=1/2, S4 alpha=1/4  (l shifts 0..k-1, u shifts 0..m-k-1)
  S5 alpha=1/2 disjoint (l shifts 0..k-1, u shifts k..m-1)

Plus exact denominator measurement for S3 via t-interpolation (small m).
"""
from fractions import Fraction
from math import lcm, log

# ---------- Apery data ----------
def apery(NMAX):
    u = [Fraction(1), Fraction(5)]
    v = [Fraction(0), Fraction(6)]
    for n in range(2, NMAX + 1):
        a = 34 * n**3 - 51 * n**2 + 27 * n - 5
        b = (n - 1) ** 3
        u.append((a * u[-1] - b * u[-2]) / n**3)
        v.append((a * v[-1] - b * v[-2]) / n**3)
    return u, v

NSEQ = 46
u, v = apery(NSEQ)
U = [int(x) for x in u]
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
Ufx = []
for n in range(NSEQ + 1):
    ln = u[n] * zeta3 - v[n]
    Lfx.append((ln.numerator * S) // ln.denominator)
    Ufx.append(U[n] * S)

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

def shape_rows(name, m):
    if name == "S1":
        return [], list(range(m))
    if name == "S2":
        k = (3 * m) // 4
    elif name in ("S3", "S5"):
        k = m // 2
    elif name == "S4":
        k = m // 4
    if name == "S5":
        return list(range(k, m)), list(range(k))
    return list(range(m - k)), list(range(k))

def kappa_bound(lshifts, m):
    tot = 0.0
    for i in lshifts:
        tot += 3 * _flog(d[i + m - 1]) + log(2)
    return tot / m**2

print("== mixed-determinant scan (P=%d) ==" % P)
print("references: control target rho=4.912 kappa=4.5 margin=0.412")
import time
for name in ("S1", "S2", "S3", "S4", "S5"):
    print("\n-- shape %s --" % name)
    print("  m   k   rho_m     2nd-diff   kappa_m   margin_m   sign")
    lh = {}
    t0 = time.time()
    for m in range(4, 23, 2):
        ush, lsh = shape_rows(name, m)
        rows = [[Ufx[i + j] for j in range(m)] for i in ush] + \
               [[Lfx[i + j] for j in range(m)] for i in lsh]
        D = det_int(rows)
        if D == 0:
            print("  m=%d: det = 0 !" % m)
            continue
        lg = _flog(abs(D)) - m * P * log(10)
        lh[m] = lg
        sd = ""
        if m - 2 in lh and m - 4 in lh:
            sd = "%8.4f" % (-(lh[m] - 2 * lh[m - 2] + lh[m - 4]) / 8)
        kb = kappa_bound(lsh, m)
        print(" %2d  %2d  %8.4f  %8s  %8.4f  %8.4f   %+d"
              % (m, len(lsh), -lg / m**2, sd, kb, -lg / m**2 - kb, 1 if D > 0 else -1))
    print("  (%.1fs)" % (time.time() - t0))

# ---------- exact denominators for S3 via t-interpolation ----------
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
        for k2 in range(len(num)):
            out[k2] += w * num[k2]
    return out

print("\n== S3 exact denominators (t-interpolation) ==")
print("  m   k   kappa_exact   kappa_bound")
for m in range(4, 15, 2):
    ush, lsh = shape_rows("S3", m)
    k = len(lsh)
    scal = [2 * d[i + m - 1] ** 3 for i in lsh]
    Sc = 1
    for x in scal:
        Sc *= x
    xs = list(range(k + 1))
    ys = []
    for t in xs:
        rows = [[U[i + j] for j in range(m)] for i in ush]
        for idx, i in enumerate(lsh):
            rows.append([int((Fraction(U[i + j] * t) - v[i + j]) * scal[idx])
                         for j in range(m)])
        ys.append(Fraction(det_int(rows), Sc))
    coeffs = lagrange([Fraction(x) for x in xs], ys)
    Dm = 1
    for c in coeffs:
        Dm = lcm(Dm, c.denominator)
    print(" %2d  %2d   %9.4f    %9.4f" % (m, k, log(Dm) / m**2, kappa_bound(lsh, m)))
