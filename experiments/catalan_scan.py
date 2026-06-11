"""
Self-designed Apery-style linear forms in (1, G), G = Catalan.

  l_n = sum_{k>=0} (-1)^k R_n(k),  u = k+1/2,
  R_n(u) = prod_{j=0}^{W-1} (u^2-(j+1/2)^2)  /  [ u^2 * prod_{j=1}^{n} (u^2-j^2)^b ]

R is EVEN in u => odd-order pole constants (pi, pi^3) cancel by parity;
double poles give only beta(2)=G tails; pole order <=3 keeps beta(4) out.
So  l_n = A_n*G - B_n  exactly, with A_n, B_n in Q computed by exact
partial fractions (Laurent series at each pole, pure Fraction arithmetic).

Measured per design: decay rate c = -log|l_n|/n, denominator growth
delta = log(den)/n, classical margin c - delta, Hankel margin
log4 + c - 1.5*delta  (window test).
"""
from fractions import Fraction
from math import lcm, log, isqrt, pi as PI_F

# ---------------- high-precision constants (fixed point, pure int) -------
P = 420
S = 10 ** P

def machin_pi():
    def atan_inv(x):
        t = S // x
        x2 = x * x
        total = 0
        k = 0
        while t:
            total += t // (2 * k + 1) if k % 2 == 0 else -(t // (2 * k + 1))
            t //= x2
            k += 1
        return total
    return 16 * atan_inv(5) - 4 * atan_inv(239)

PI = machin_pi()
SQRT3 = isqrt(3 * 10 ** (2 * P))                 # sqrt(3)*10^P
def log_2_plus_sqrt3():
    # = 2*atanh(1/sqrt3) = (2/sqrt3) * sum 3^-k/(2k+1)
    t = 2 * S * S // SQRT3
    total = 0
    k = 0
    while t:
        total += t // (2 * k + 1)
        t //= 3
        k += 1
    return total
L23 = log_2_plus_sqrt3()

def catalan_fx():
    # G = (3/8) sum_{n>=0} 1/(C(2n,n)(2n+1)^2) + (pi/8) log(2+sqrt3)
    total = 0
    c = 1                                        # C(2n,n)
    n = 0
    while True:
        term = S // (c * (2 * n + 1) ** 2)
        if term == 0:
            break
        total += term
        c = c * (2 * n + 1) * (2 * n + 2) // ((n + 1) * (n + 1))
        n += 1
    return (3 * total + (PI * L23) // S) // 8
G_FX = catalan_fx()
assert str(G_FX)[:18] == "915965594177219015", str(G_FX)[:20]
print("G =", "0." + str(G_FX)[:30], "... (validated 18 digits)")

# ---------------- truncated series over Fraction --------------------------
def smul(a, b, L):
    out = [Fraction(0)] * L
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if i + j >= L:
                break
            out[i + j] += x * y
    return out

def sinv(a, L):
    assert a[0] != 0
    out = [Fraction(0)] * L
    out[0] = 1 / a[0]
    for k in range(1, L):
        acc = Fraction(0)
        for i in range(1, k + 1):
            acc += a[i] * out[k - i]
        out[k] = -acc / a[0]
    return out

# ---------------- tails V_r(p) --------------------------------------------
# V_r(p) = sum_{k>=0} (-1)^k / (k+1/2-p)^r  =  rat + sign*K_r
# K_1 = pi/2, K_2 = 4G, K_3 = pi^3/4  (track coefficients separately)
def V(r, p):
    rat = Fraction(0)
    if p <= 0:
        q = -p
        sign = -1 if q % 2 else 1
        # rat part: -sign * sum_{m<q} (-1)^m/(m+1/2)^r
        for m in range(q):
            t = Fraction((-1) ** m * 2 ** r, (2 * m + 1) ** r)
            rat -= sign * t
        return rat, sign
    sign = -1 if p % 2 else 1
    for k in range(p):
        base = Fraction(2, 2 * k + 1 - 2 * p)    # 1/(k+1/2-p)
        t = base ** r
        if k % 2:
            t = -t
        rat += t
    return rat, sign

# ---------------- exact partial fractions for one (n, b, W) ---------------
def linear_form(n, b, W):
    num_roots = []
    for j in range(W):
        h = Fraction(2 * j + 1, 2)
        num_roots += [h, -h]
    den_roots = [(Fraction(0), 2)]
    for i in range(1, n + 1):
        den_roots += [(Fraction(i), b), (Fraction(-i), b)]

    cK = [Fraction(0)] * 4                       # coeffs of K_1..K_3 (1-indexed)
    rat_total = Fraction(0)
    for p in range(-n, n + 1):
        pf = Fraction(p)
        e = 2 if p == 0 else b
        L = e
        ser = [Fraction(1)] + [Fraction(0)] * (L - 1)
        for rho in num_roots:
            ser = smul(ser, [pf - rho, Fraction(1)][:L] + [Fraction(0)] * max(0, L - 2), L)
        dser = [Fraction(1)] + [Fraction(0)] * (L - 1)
        for rho, m in den_roots:
            if rho == pf:
                continue
            fac = [pf - rho, Fraction(1)][:L] + [Fraction(0)] * max(0, L - 2)
            for _ in range(m):
                dser = smul(dser, fac, L)
        ser = smul(ser, sinv(dser, L), L)
        # c_r = ser[e-r],  contribution c_r * V_r(p)
        for r in range(1, e + 1):
            c = ser[e - r]
            if c == 0:
                continue
            rat, sign = V(r, p)
            rat_total += c * rat
            cK[r] += c * sign
    assert cK[1] == 0, ("pi survives!", n, b, W)
    if len(cK) > 3:
        pass
    assert cK[3] == 0 if b >= 3 or True else True
    assert cK[3] == 0, ("pi^3 survives!", n, b, W)
    A = 4 * cK[2]                                # coeff of G
    B = -rat_total
    return A, B                                  # l_n = A*G - B

def _flog(x):
    bl = x.bit_length()
    if bl <= 900:
        return log(x)
    sh = bl - 60
    return log(x >> sh) + sh * log(2)

def ell_log(A, B):
    # log |A*G - B| via fixed point
    val = A * Fraction(G_FX, S) - B
    if val == 0:
        return None, 0
    return _flog(abs(val.numerator)) - _flog(val.denominator), 1 if val > 0 else -1

# ---------------- numeric validation (n small) ----------------------------
def numeric_sum(n, b, W, terms=400000):
    G = G_FX / S
    tot = 0.0
    prev = 0.0
    for k in range(terms):
        u = k + 0.5
        u2 = u * u
        num = 1.0
        for j in range(W):
            num *= (u2 - (j + 0.5) ** 2)
        den = u2
        for j in range(1, n + 1):
            den *= (u2 - j * j) ** b
        t = num / den
        prev = tot
        tot += t if k % 2 == 0 else -t
    return (tot + prev) / 2                      # Cesaro pairing

print("\n== validation (n=2) ==")
for (b, Wf) in [(1, 1), (2, 2), (3, 3)]:
    n = 2
    W = Wf * n
    A, B = linear_form(n, b, W)
    exact = float(A) * (G_FX / S) - float(B)
    approx = numeric_sum(n, b, W)
    rel = abs(exact - approx) / max(1e-300, abs(exact))
    print("b=%d W=%dn: exact=%.10e  numeric=%.10e  rel.err=%.1e" % (b, Wf, exact, approx, rel))
    assert rel < 1e-4, "pipeline mismatch!"

# ---------------- scan ----------------------------------------------------
print("\n== design scan:  l_n = A_n G - B_n ==")
DESIGNS = [(1, 1, 56), (2, 2, 40), (3, 3, 30)]
for b, Wf, NMAX in DESIGNS:
    data = []
    E = 1                                        # lcm of denominators so far
    for n in range(1, NMAX + 1):
        W = Wf * n
        A, B = linear_form(n, b, W)
        lg, sgn = ell_log(A, B)
        if lg is None:
            print("  n=%d: l_n = 0 exactly (skipping)" % n)
            continue
        E = lcm(E, B.denominator, A.denominator)
        data.append((n, lg, _flog(E), sgn, A.denominator))
    print("\n-- design b=%d, W=%dn --" % (b, Wf))
    print("  n   -log|l_n|/n    logE_n/n    sign  logden(A)/n")
    for (n, lg, lE, sgn, dA) in data:
        if n % max(1, NMAX // 10) == 0 or n == NMAX:
            print(" %3d   %9.5f    %9.5f     %+d    %8.4f"
                  % (n, -lg / n, lE / n, sgn, _flog(dA) / n if dA > 1 else 0.0))
    n2, n1 = data[-1], data[-11]
    c = -(n2[1] - n1[1]) / (n2[0] - n1[0])
    dE = (n2[2] - n1[2]) / (n2[0] - n1[0])
    print("  tail slopes: decay c = %.5f , denom delta = %.5f" % (c, dE))
    print("  classical margin  c - delta          = %+.5f" % (c - dE))
    print("  Hankel    margin  log4 + c - 1.5*d   = %+.5f" % (log(4) + c - 1.5 * dE))
