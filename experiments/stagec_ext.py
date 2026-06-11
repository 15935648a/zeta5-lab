"""Extend Stage C to m=18: is log(S_m/D_m) linear in m (savings -> 0 per m^2)
or quadratic (true arithmetic dividend)?"""
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

NSEQ = 40
u, v = apery(NSEQ)
U = [int(x) for x in u]
d = [1]
for k in range(1, NSEQ + 1):
    d.append(lcm(d[-1], k))

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

print(" m   cost=logD/m^2   bound=logS/m^2   gap=log(S/D)   gap increment")
prev_gap = None
gaps = {}
for m in range(2, 19):
    s = [2 * d[m - 1 + j] ** 3 for j in range(m)]
    S = 1
    for x in s:
        S *= x
    xs = list(range(m + 1))
    ys = []
    for t in xs:
        Mt = [[int((Fraction(U[i + j] * t) - v[i + j]) * s[j]) for j in range(m)]
              for i in range(m)]
        ys.append(Fraction(det_int(Mt), S))
    coeffs = lagrange([Fraction(x) for x in xs], ys)
    assert coeffs[m].denominator == 1
    Dm = 1
    for c in coeffs:
        Dm = lcm(Dm, c.denominator)
    gap = log(S) - log(Dm)
    gaps[m] = gap
    inc = "" if prev_gap is None else "%.3f" % (gap - prev_gap)
    prev_gap = gap
    print("%2d     %8.5f       %8.5f       %8.3f      %s"
          % (m, log(Dm) / m**2, log(S) / m**2, gap, inc))

# least-squares fit gap ~ a*m + b and gap ~ c*m^2 + a*m + b on m>=8
ms = [m for m in gaps if m >= 8]
ys = [gaps[m] for m in ms]
n = len(ms)
sx = sum(ms); sy = sum(ys); sxx = sum(m * m for m in ms); sxy = sum(m * y for m, y in zip(ms, ys))
a = (n * sxy - sx * sy) / (n * sxx - sx * sx)
b = (sy - a * sx) / n
print("\nlinear fit (m>=8): gap = %.3f*m + %.3f" % (a, b))
resid = [y - (a * m + b) for m, y in zip(ms, ys)]
print("residuals:", " ".join("%.2f" % r for r in resid))
print("=> per-m^2 savings at m=18: %.4f ; extrapolated m=40: %.4f ; m->inf: 0 if linear"
      % (gaps[18] / 324, (a * 40 + b) / 1600))
