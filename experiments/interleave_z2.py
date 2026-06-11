"""
Experiment 1: Hankel capacity of an interleaved Apery path for zeta(2).

Family: I(n1,n2,m) = II x^n1(1-x)^n1 y^n2(1-y)^n2 / (1-xy)^(m+1) dxdy
      = sum_{k>=0} R(k),
  R(k) = [n1! n2! / m!] * prod_{i=1}^m (k+i)
         / [ prod_{j=n1+1}^{2n1+1}(k+j) * prod_{j=n2+1}^{2n2+1}(k+j) ]
Partial fractions: double poles -> zeta(2) tails, simple poles -> harmonic
rationals (sum of residues = 0).   l = A*zeta(2) + R_rat  exactly.

Path A: t=2s -> (s,s,s) ; t=2s+1 -> (s+1,s,s).

Measured: per-step decay c~, per-step denominator delta~ (even/odd splits),
Hankel rate of the path sequence.  Discriminator:
  slit-conservation law  => rate = c~ + log 2
  fine-graining loophole => rate = c~ + log 4
Also: diagonal-only Hankel (coarse) should show rate = 2c~ + log 4 = 3.7924.
"""
from fractions import Fraction
from math import lcm, log, factorial

# ---------- fixed-point constants ----------
P = 1600
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
Z2_FX = PI * PI // (6 * S)          # zeta(2)*10^P
assert str(Z2_FX)[:12] == "164493406684", str(Z2_FX)[:14]

def _flog(x):
    bl = x.bit_length()
    if bl <= 900:
        return log(x)
    sh = bl - 60
    return log(x >> sh) + sh * log(2)

# ---------- series helpers over Fraction ----------
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

NH = 400
H1 = [Fraction(0)] * (NH + 1)
H2 = [Fraction(0)] * (NH + 1)
for i in range(1, NH + 1):
    H1[i] = H1[i - 1] + Fraction(1, i)
    H2[i] = H2[i - 1] + Fraction(1, i * i)

def lform(n1, n2, m):
    """l(n1,n2,m) = A*zeta(2) + rat (exact)."""
    poles = {}
    for j in range(n1 + 1, 2 * n1 + 2):
        poles[j] = poles.get(j, 0) + 1
    for j in range(n2 + 1, 2 * n2 + 2):
        poles[j] = poles.get(j, 0) + 1
    const = Fraction(factorial(n1) * factorial(n2), factorial(m))
    A = Fraction(0)
    rat = Fraction(0)
    Bsum = Fraction(0)
    for j, e in poles.items():
        L = e
        pad = [Fraction(0)] * max(0, L - 2)
        ser = [const] + [Fraction(0)] * (L - 1)
        for i in range(1, m + 1):
            ser = smul(ser, [Fraction(i - j), Fraction(1)][:L] + pad, L)
        dser = [Fraction(1)] + [Fraction(0)] * (L - 1)
        for j2, e2 in poles.items():
            if j2 == j:
                continue
            fac = [Fraction(j2 - j), Fraction(1)][:L] + pad
            for _ in range(e2):
                dser = smul(dser, fac, L)
        ser = smul(ser, sinv(dser, L), L)
        c1 = ser[e - 1]
        c2 = ser[e - 2] if e >= 2 else Fraction(0)
        A += c2
        rat += -c2 * H2[j - 1] - c1 * H1[j - 1]
        Bsum += c1
    assert Bsum == 0, ("residue sum nonzero", n1, n2, m)
    return A, rat

# ---------- validation ----------
print("== validation ==")
A, r = lform(1, 1, 1)
print("(1,1,1): A=%s (Apery b_1=3, sign free), l=%.8e" % (A, float(A) * (Z2_FX / S) + float(r)))
def numeric_sum(n1, n2, m, terms=200000):
    tot = 0.0
    for k in range(terms):
        t = factorial(n1) * factorial(n2) / factorial(m)
        for i in range(1, m + 1):
            t *= (k + i)
        for j in range(n1 + 1, 2 * n1 + 2):
            t /= (k + j)
        for j in range(n2 + 1, 2 * n2 + 2):
            t /= (k + j)
        tot += t
    return tot
for (n1, n2, m) in [(2, 2, 2), (3, 2, 2), (3, 2, 3)]:
    A, r = lform(n1, n2, m)
    ex = float(A) * (Z2_FX / S) + float(r)
    nu = numeric_sum(n1, n2, m)
    print("(%d,%d,%d): exact=%.10e numeric=%.10e rel=%.1e" % (n1, n2, m, ex, nu, abs(ex - nu) / abs(ex)))

# ---------- path sequence ----------
TMAX = 52
def path_tuple(t):
    s = t // 2
    return (s, s, s) if t % 2 == 0 else (s + 1, s, s)

ells = []          # Fraction value of l~_t  (using fixed-point zeta2)
Z2F = Fraction(Z2_FX, S)
E = 1
logE = []
print("\n== path A:  t, -log|l|/t, logE_t/t, dE increment ==")
prevE = 0.0
for t in range(0, TMAX + 1):
    n1, n2, m = path_tuple(t)
    A, r = lform(n1, n2, m)
    val = A * Z2F + r
    ells.append(val)
    E = lcm(E, A.denominator, r.denominator)
    lE = _flog(E)
    logE.append(lE)
    if t >= TMAX - 12 or t % 8 == 0:
        lg = _flog(abs(val.numerator)) - _flog(val.denominator) if val != 0 else float("nan")
        print(" t=%3d (%2d,%2d,%2d)  %8.5f   %8.5f   +%.4f"
              % (t, n1, n2, m, -lg / max(t, 1), lE / max(t, 1), lE - prevE))
    prevE = lE

# per-step decay from tail
lg_hi = _flog(abs(ells[TMAX].numerator)) - _flog(ells[TMAX].denominator)
lg_lo = _flog(abs(ells[TMAX - 16].numerator)) - _flog(ells[TMAX - 16].denominator)
c_step = -(lg_hi - lg_lo) / 16
# even/odd denominator increments (tail average)
inc_e = [logE[t] - logE[t - 1] for t in range(TMAX - 15, TMAX + 1) if t % 2 == 0]
inc_o = [logE[t] - logE[t - 1] for t in range(TMAX - 15, TMAX + 1) if t % 2 == 1]
d_step = (logE[TMAX] - logE[TMAX - 16]) / 16
print("\nper-step decay  c~ = %.5f   (half of 5*log(phi)=2.40606 -> 1.20303)" % c_step)
print("per-step denom  d~ = %.5f   (even-step avg %.4f, odd-step avg %.4f)"
      % (d_step, sum(inc_e) / len(inc_e), sum(inc_o) / len(inc_o)))

# ---------- Hankel of the path sequence ----------
def det_int(M):
    M = [row[:] for row in M]
    n = len(M)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r2 in range(k + 1, n):
                if M[r2][k] != 0:
                    M[k], M[r2] = M[r2], M[k]
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

fx = [(v.numerator * S) // v.denominator for v in ells]
print("\n== Hankel of path sequence ==")
print("predictions: slit-law c~+log2 = %.5f ; loophole c~+log4 = %.5f"
      % (c_step + log(2), c_step + log(4)))
print(" m    log|H_m|/m^2   2nd-diff/2   nonzero")
lh = {}
MMAX = (TMAX // 2)
for mm in range(1, MMAX + 1):
    Mx = [[fx[i + j] for j in range(mm)] for i in range(mm)]
    D = det_int(Mx)
    if D == 0:
        print(" m=%d: H=0" % mm)
        continue
    l = _flog(abs(D)) - mm * P * log(10)
    lh[mm] = l
    sd = ""
    if mm >= 3 and mm - 1 in lh and mm - 2 in lh:
        sd = "%.5f" % ((lh[mm] - 2 * lh[mm - 1] + lh[mm - 2]) / 2)
    print(" %2d    %10.5f    %10s   %s" % (mm, l / mm**2, sd, D != 0))

# ---------- coarse (diagonal-only) Hankel for comparison ----------
print("\n== Hankel of diagonal subsequence (coarse, prediction 2c~+log4 = %.5f) ==" % (2 * c_step + log(4)))
diag = [ells[2 * s] for s in range(0, TMAX // 2 + 1)]
fxd = [(v.numerator * S) // v.denominator for v in diag]
lh2 = {}
for mm in range(1, len(fxd) // 2 + 1):
    Mx = [[fxd[i + j] for j in range(mm)] for i in range(mm)]
    D = det_int(Mx)
    if D == 0:
        continue
    l = _flog(abs(D)) - mm * P * log(10)
    lh2[mm] = l
    sd = ""
    if mm >= 3 and mm - 1 in lh2 and mm - 2 in lh2:
        sd = "%.5f" % ((lh2[mm] - 2 * lh2[mm - 1] + lh2[mm - 2]) / 2)
    print(" %2d    %10.5f    %10s" % (mm, l / mm**2, sd))

# ---------- ledger ----------
print("\n== ledger (per fine step) ==")
cost = 1.5 * d_step
print(" Hankel cost 1.5*d~          = %.5f" % cost)
print(" capacity if slit-law        = %.5f  -> margin %+0.5f" % (c_step + log(2), c_step + log(2) - cost))
print(" capacity if loophole        = %.5f  -> margin %+0.5f" % (c_step + log(4), c_step + log(4) - cost))
print(" classical diagonal margin (ref) = +0.40606 per coarse step = +0.20303 per fine")
