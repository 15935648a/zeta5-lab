"""
Gaussian-integer route for Catalan G.

Main family   : I_n(i) = II [x(1-x)y(1-y)]^n/(1-ixy)^{n+1}
                = alpha*Li2(i) + beta*Lam + gamma  over Q(i),
                Li2(i) = -pi^2/48 + i G,  Lam = -log(1-i) = -(1/2)log2 + i pi/4.
Helpers       : zeta(2) diagonal forms (kill pi^2),
                J_m = int [x(1-x)]^m/(1-ix)^{m+1} dx  (kill log2 and pi).
Elimination   : exact rational linear algebra per n -> pure forms  b + a*G.
Measured      : decay c_E, denominators delta_E, classical margin, plus the
                primitive rates (c_i, c_Lam, coefficient growths) for diagnosis.
"""
from fractions import Fraction
from math import lcm, log, factorial, isqrt

# ---------------- fixed-point constants ----------------
P = 1200
S = 10 ** P

def machin_pi():
    def atan_inv(x):
        t = S // x
        x2 = x * x
        tot = 0
        k = 0
        while t:
            tot += t // (2 * k + 1) if k % 2 == 0 else -(t // (2 * k + 1))
            t //= x2
            k += 1
        return tot
    return 16 * atan_inv(5) - 4 * atan_inv(239)

PI = machin_pi()

def log_int(a, b):
    # log(a/b) for small ints via atanh: log(x) = 2 atanh((x-1)/(x+1))
    num, den = a - b, a + b
    t = S * num // den
    x2n, x2d = num * num, den * den
    tot = 0
    k = 0
    while t:
        tot += t // (2 * k + 1)
        t = t * x2n // x2d
        k += 1
    return 2 * tot

LOG2 = log_int(2, 1)
assert str(LOG2)[:10] == "6931471805"

SQRT3 = isqrt(3 * 10 ** (2 * P))
def log_2_plus_sqrt3():
    t = 2 * S * S // SQRT3
    tot = 0
    k = 0
    while t:
        tot += t // (2 * k + 1)
        t //= 3
        k += 1
    return tot

def catalan_fx():
    tot = 0
    c = 1
    n = 0
    while True:
        term = S // (c * (2 * n + 1) ** 2)
        if term == 0:
            break
        tot += term
        c = c * (2 * n + 1) * (2 * n + 2) // ((n + 1) * (n + 1))
        n += 1
    return (3 * tot + (PI * log_2_plus_sqrt3()) // S) // 8

G_FX = catalan_fx()
assert str(G_FX)[:18] == "915965594177219015"
PI2_FX = PI * PI // S
print("constants ok: pi, log2, G, pi^2 at %d digits" % P)

CONST_FX = {"one": S, "G": G_FX, "pi2": PI2_FX, "log2": LOG2, "pi": PI}
BASIS = ["one", "G", "pi2", "log2", "pi"]

def _flog(x):
    bl = x.bit_length()
    if bl <= 900:
        return log(x)
    sh = bl - 60
    return log(x >> sh) + sh * log(2)

def form_value_log(vec):
    # vec: dict basis->Fraction ; returns log|value|
    total = Fraction(0)
    for b, c in vec.items():
        if c:
            total += c * Fraction(CONST_FX[b], S)
    if total == 0:
        return None
    return _flog(abs(total.numerator)) - _flog(total.denominator)

def vec_den(vec):
    d = 1
    for c in vec.values():
        d = lcm(d, c.denominator)
    return d

# ---------------- series helpers ----------------
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
    out = [Fraction(0)] * L
    out[0] = 1 / a[0]
    for k in range(1, L):
        acc = Fraction(0)
        for i in range(1, k + 1):
            acc += a[i] * out[k - i]
        out[k] = -acc / a[0]
    return out

# ---------------- Q(i) helpers: pairs (re, im) of Fractions ----------------
def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])

def cmulq(a, q):          # complex * rational
    return (a[0] * q, a[1] * q)

def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

IPOW = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)),
        (Fraction(-1), Fraction(0)), (Fraction(0), Fraction(-1))]
def ipow(e):
    return IPOW[e % 4]

# partial sums S1(j) = sum_{m<j} i^m/m, S2(j) = sum i^m/m^2 (over Q(i))
NMAXP = 140
S1 = [(Fraction(0), Fraction(0))]
S2 = [(Fraction(0), Fraction(0))]
for m in range(1, NMAXP + 1):
    S1.append(cadd(S1[-1], cmulq(ipow(m), Fraction(1, m))))
    S2.append(cadd(S2[-1], cmulq(ipow(m), Fraction(1, m * m))))
# pad index 0/1 alignment: S1[j] = sum_{m=1}^{j-1}? fix: rebuild so S*[j] = sum_{m<j}
S1 = [(Fraction(0), Fraction(0))] + S1[:-1]
S2 = [(Fraction(0), Fraction(0))] + S2[:-1]
# now S1[j] = sum_{m=1}^{j-1} i^m/m

# ---------------- main family I_n(i) ----------------
def main_form(n):
    """I_n(i) -> (alpha, beta, gamma) in Q(i); Q(k)=(k+1..k+n)/n! *
       [n!/((k+n+1)..(k+2n+1))]^2, double poles j=n+1..2n+1."""
    alpha = (Fraction(0), Fraction(0))
    beta = (Fraction(0), Fraction(0))
    gamma = (Fraction(0), Fraction(0))
    for j in range(n + 1, 2 * n + 2):
        L = 2
        ser = [Fraction(factorial(n)), Fraction(0)]   # n!^2/n! = n!
        for i in range(1, n + 1):
            ser = smul(ser, [Fraction(i - j), Fraction(1)], L)
        dser = [Fraction(1), Fraction(0)]
        for j2 in range(n + 1, 2 * n + 2):
            if j2 == j:
                continue
            fac = [Fraction(j2 - j), Fraction(1)]
            dser = smul(dser, fac, L)
            dser = smul(dser, fac, L)
        ser = smul(ser, sinv(dser, L), L)
        A_j, B_j = ser[0], ser[1]
        w = ipow(-j)                                   # i^{-j}
        alpha = cadd(alpha, cmulq(w, A_j))
        beta = cadd(beta, cmulq(w, B_j))
        gamma = cadd(gamma, cmul(w, cadd(cmulq(S2[j], -A_j), cmulq(S1[j], -B_j))))
    return alpha, beta, gamma

# ---------------- Lambda helper J_m ----------------
def lam_form(m):
    """J_m = beta~*Lam + gamma~ ; Q(k) = (k+1..k+m)/[(k+m+1)..(k+2m+1)],
       simple poles."""
    beta = (Fraction(0), Fraction(0))
    gamma = (Fraction(0), Fraction(0))
    for j in range(m + 1, 2 * m + 2):
        num = Fraction(1)
        for i in range(1, m + 1):
            num *= (i - j)
        den = Fraction(1)
        for j2 in range(m + 1, 2 * m + 2):
            if j2 != j:
                den *= (j2 - j)
        B_j = num / den
        w = ipow(-j)
        beta = cadd(beta, cmulq(w, B_j))
        gamma = cadd(gamma, cmul(w, cmulq(S1[j], -B_j)))
    return beta, gamma

# ---------------- zeta(2) helper (pure real) ----------------
def z2_form(m):
    A = Fraction(0)
    rat = Fraction(0)
    Bsum = Fraction(0)
    H1v = Fraction(0)
    H2v = Fraction(0)
    H1l = [Fraction(0)]
    H2l = [Fraction(0)]
    for i in range(1, 2 * m + 2):
        H1v += Fraction(1, i)
        H2v += Fraction(1, i * i)
        H1l.append(H1v)
        H2l.append(H2v)
    for j in range(m + 1, 2 * m + 2):
        L = 2
        ser = [Fraction(factorial(m)), Fraction(0)]
        for i in range(1, m + 1):
            ser = smul(ser, [Fraction(i - j), Fraction(1)], L)
        dser = [Fraction(1), Fraction(0)]
        for j2 in range(m + 1, 2 * m + 2):
            if j2 == j:
                continue
            fac = [Fraction(j2 - j), Fraction(1)]
            dser = smul(dser, fac, L)
            dser = smul(dser, fac, L)
        ser = smul(ser, sinv(dser, L), L)
        A += ser[0]
        rat += -ser[0] * H2l[j - 1] - ser[1] * H1l[j - 1]
        Bsum += ser[1]
    assert Bsum == 0
    return A, rat            # l = A*zeta2 + rat, zeta2 = pi^2/6

# ---------------- to basis vectors ----------------
def main_vecs(n):
    al, be, ga = main_form(n)
    re = {"one": ga[0], "G": -al[1], "pi2": -al[0] / 48,
          "log2": -be[0] / 2, "pi": -be[1] / 4}
    im = {"one": ga[1], "G": al[0], "pi2": -al[1] / 48,
          "log2": -be[1] / 2, "pi": be[0] / 4}
    return re, im

def lam_vecs(m):
    be, ga = lam_form(m)
    re = {"one": ga[0], "G": Fraction(0), "pi2": Fraction(0),
          "log2": -be[0] / 2, "pi": -be[1] / 4}
    im = {"one": ga[1], "G": Fraction(0), "pi2": Fraction(0),
          "log2": -be[1] / 2, "pi": be[0] / 4}
    return re, im

def z2_vec(m):
    A, rat = z2_form(m)
    return {"one": rat, "G": Fraction(0), "pi2": A / 6,
            "log2": Fraction(0), "pi": Fraction(0)}

# ---------------- validation ----------------
print("\n== validation (n=3) ==")
al, be, ga = main_form(3)
exact = complex(float(al[0]) + 1j * float(al[1])) * complex(-(PI2_FX / S) / 48, G_FX / S) \
      + complex(float(be[0]) + 1j * float(be[1])) * complex(-(LOG2 / S) / 2, (PI / S) / 4) \
      + complex(float(ga[0]) + 1j * float(ga[1]))
num = 0j
for k in range(200000):
    t = 1.0
    for i in range(1, 4):
        t *= (k + i)
    t /= factorial(3)
    f = factorial(3)
    for j in range(4, 8):
        f /= (k + j)
    t *= f * f
    num += (1j) ** (k % 4) * t
print("I_3(i): exact=%.10e%+.10ej  numeric=%.10e%+.10ej" % (exact.real, exact.imag, num.real, num.imag))
assert abs(exact - num) / abs(exact) < 1e-6

# ---------------- measure primitives ----------------
print("\n== primitive rates ==")
NMAX = 30
prim = {}
for n in (10, 20, 30):
    re, im = main_vecs(n)
    lg = form_value_log(re)
    lgi = form_value_log(im)
    al, be, ga = main_form(n)
    galpha = _flog(abs(al[0].numerator) + abs(al[1].numerator)) - _flog(al[0].denominator)
    print(" n=%2d  -log|ReI|/n=%.4f  -log|ImI|/n=%.4f  g(alpha)/n=%.4f"
          % (n, -lg / n, -lgi / n, galpha / n))
for m in (15, 30, 45):
    if m > NMAXP // 2 - 1:
        continue
    re, im = lam_vecs(m)
    lg = form_value_log(re)
    lgi = form_value_log(im)
    print(" Lam m=%2d  -log|ReJ|/m=%.4f  -log|ImJ|/m=%.4f" % (m, -lg / m, -lgi / m))

# ---------------- exact elimination per n ----------------
print("\n== eliminated pure (1,G) forms ==")
print("  n  (mz2,mL)   -log|E|/n   logden/n   [two basis solutions]")
results = []
E_lcm = 1
for n in range(2, NMAX + 1):
    # helper index matching: decay of main ~ c_i*n; z2 decay 2.406/m ; Lam decay c_L/m
    mz2 = max(1, round(n * 2.80 / 2.406))
    mL = max(1, min(NMAXP // 2 - 1, round(n * 2.80 / 1.49)))
    vecs = []
    re, im = main_vecs(n)
    vecs += [re, im]
    vecs.append(z2_vec(mz2))
    reL, imL = lam_vecs(mL)
    vecs += [reL, imL]
    # solve: sum lam_i vec_i has pi2 = log2 = pi = 0
    # unknowns lam in Q^5; build 3x5 matrix
    M = [[vecs[i][b] for i in range(5)] for b in ("pi2", "log2", "pi")]
    # find 2D nullspace by exact Gaussian elimination
    rows = [r[:] for r in M]
    ncol = 5
    piv = []
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
        piv.append(c)
        rr += 1
        if rr == len(rows):
            break
    free = [c for c in range(ncol) if c not in piv]
    sols = []
    for fc in free:
        lamv = [Fraction(0)] * ncol
        lamv[fc] = Fraction(1)
        for r, pc in enumerate(piv):
            lamv[pc] = -rows[r][fc]
        sols.append(lamv)
    best = None
    for lamv in sols:
        # clear solver denominators: integer lambda, so den(E) comes from
        # the FORMS only; the size of lambda is the honest elimination tax
        dl = 1
        for x in lamv:
            dl = lcm(dl, x.denominator)
        lamv = [x * dl for x in lamv]
        E = {"one": Fraction(0), "G": Fraction(0)}
        for i, lv in enumerate(lamv):
            if lv == 0:
                continue
            for b in ("one", "G"):
                E[b] += lv * vecs[i][b]
        chk = {b: sum(lv * vecs[i][b] for i, lv in enumerate(lamv))
               for b in ("pi2", "log2", "pi")}
        assert all(v == 0 for v in chk.values())
        if E["G"] == 0:
            continue
        # normalize: make integral-ish by clearing denominators? measure as-is
        lg = form_value_log(E)
        dn = vec_den(E)
        # scale-invariant measure: normalize by content of (one,G) coeffs
        from math import gcd
        sc = Fraction(gcd(E["one"].numerator * E["G"].denominator,
                          E["G"].numerator * E["one"].denominator),
                      E["one"].denominator * E["G"].denominator) if E["one"] else None
        cand = (lg, dn, E)
        if best is None or lg < best[0]:
            best = cand
    if best is None:
        continue
    lg, dn, E = best
    dn2 = lcm(E["one"].denominator, E["G"].denominator)
    E_lcm = lcm(E_lcm, dn2)
    results.append((n, lg, _flog(dn2), _flog(E_lcm)))
    if n % 4 == 0 or n >= NMAX - 4:
        print(" %2d  (%2d,%2d)   %8.4f   %8.4f   (E-filtration %.4f)"
              % (n, mz2, mL, -results[-1][1] / n, results[-1][2] / n, results[-1][3] / n))

if len(results) >= 12:
    n2 = results[-1]
    n1 = results[-9]
    span = n2[0] - n1[0]
    cE = -(n2[1] - n1[1]) / span
    dE = (n2[3] - n1[3]) / span
    print("\n tail slopes: c_E = %.4f ; delta_E = %.4f" % (cE, dE))
    print(" classical margin c_E - delta_E = %+.4f" % (cE - dE))
    print(" (odd-lattice route reference deficit: about -6.3)")
