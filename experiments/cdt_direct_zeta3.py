"""
CANDIDATE (to be adversarially verified): direct application of CDT
Theorem 2.5.1 to Apery's zeta(3) data.

Setup: punctures {0, x1, x2, oo}, x1 = (sqrt2-1)^4, x2 = 1/x1 = (1+sqrt2)^4.
Under zeta(3) = p/q: P = qV - pA in Q[[x]], type [1..n]^3 (2*d_n^3*v_n in Z),
holomorphic on the slit plane C \\ [x2, oo)  (overconvergent at 0 and x1;
verified numerically to capacity precision in calib.py, session 1).
Apery ODE has order 3, irreducible => {P, P', P''} Q(x)-linearly independent;
derivatives preserve [1..bn]-types asymptotically (CDT p.8).

Map: scaled bivalent  phi(z) = x2 * 8(z+z^3)/(1+z)^4   (CDT 2.11.2):
  - misses {x2, oo} on D (Sigma^1), covers Sigma^0 = {0, x1} once,
    preimages in the univalent interval (-1,1)  (CDT Prop 2.9.3 setup,
    same as their showcase 2.11 with Sigma^0 = {0, delta, delta/(delta-1)});
  - phi'(0) = 8*x2;  BC double integral = log(8 x2) + 4G/pi  (Lemma 2.11.7
    + scaling adds log x2 to both numerator and log phi'(0)).

Bound (2.5.4):  m <= BC / (log phi'(0) - tau(b)).
"""
from math import log, sqrt, pi

G = 0.915965594177219015054603514932384110774
x2 = (1 + sqrt(2)) ** 4
BCbase = log(8.0) + 4 * G / pi          # Lemma 2.11.7
BC = log(x2) + BCbase
LP = log(8 * x2)                         # log phi'(0)

def tau(types):
    m = len(types)
    s = sorted(types)
    return sum((2 * i + 1) * s[i] for i in range(m)) / m**2

print("x2 = %.6f, log x2 = %.5f" % (x2, log(x2)))
print("BC = log(8 x2) + 4G/pi = %.5f ; log phi'(0) = %.5f" % (BC, LP))
print("sigma_m condition: phi'(0) = %.1f > e^3 = %.2f : %s"
      % (8 * x2, 2.718281828 ** 3, 8 * x2 > 2.718281828 ** 3))

print("\n m   tuple              tau      bound    verdict")
for m, types, label in [
    (2, [0, 3], "{1, P}        "),
    (3, [0, 3, 3], "{1, P, P'}    "),
    (4, [0, 3, 3, 3], "{1, P, P', P''}"),
]:
    t = tau(types)
    bd = BC / (LP - t)
    verdict = "CONTRADICTION" if m > bd else "allowed (margin %.3f)" % (bd - m)
    print("  %d  %s  %.4f   %.4f   %s" % (m, label, t, bd, verdict))

print("""
Sanity check: m=2 is allowed by a hair (2.018 vs 2) -- the bound 'knows'
Apery's P exists classically. The contradiction at m=3,4 would be a proof
of zeta(3) irrationality. TOO STRONG TO TRUST until adversarially verified.
""")

# ---- ADVERSARIAL VERIFICATION ROUND 1 (univalent form, Thm 2.7.1) ----
# Univalent Koebe map onto Omega = C \ [x2, oo): rho(Omega,0) = 4*x2.
# m <= log rho / (log rho - tau).  Threshold: lambda > tau*m/(m-1).
lam = log(4 * x2)
print("\n== univalent (Thm 2.7.1, Koebe onto slit plane), lambda = log(4 x2) = %.5f ==" % lam)
for m, types, label in [
    (2, [0, 3], "{1, P}        "),
    (3, [0, 3, 3], "{1, P, P'}    "),
    (4, [0, 3, 3, 3], "{1, P, P', P''}"),
]:
    t = tau(types)
    bd = lam / (lam - t)
    thr = t * m / (m - 1)
    print("  m=%d %s tau=%.4f bound=%.4f  threshold lam>%.3f (x2>%.1f)  %s"
          % (m, label, t, bd, thr, 2.718281828 ** thr / 4,
             "CONTRADICTION" if m > bd else "allowed"))
print("""
 KEY IDENTIFICATION: the m=2 univalent case is EXACTLY Zudilin's
 determinantal criterion / our session-1 Hankel proof:
   contradiction iff log(4 x2) > 4.5  <=>  capacity 4.912 > cost 4.5,
 margin 0.412 -- three formalisms, one inequality. VERIFIED CONSISTENT.
 Everything beyond Zudilin hinges on ONE question: is the derivative
 stack {P, P', P''} legal in tau(b)?  (CDT Remark 2.7.9's 'm such
 functions -> m/(m+1)' explicitly suggests YES; their silence on a
 resulting 5-line zeta(3) proof suggests NO. Decisive: read the
 denominator-arithmetic proofs, Sections 6.6 / 7.)
""")

# zeta(5) target curve under the same template (IF the template is valid):
# rank-k system, conditional overconvergent P5 of type [1..n]^5, decay x2.
print(" zeta(5) target spec (k = usable derivative stack = ODE rank):")
print("  k   m    tau      required x2")
for k in (2, 3, 4, 5, 6, 8, 12):
    m = k + 1
    t = 5.0 * (m * m - 1) / (m * m)
    # need m > (L + 4G/pi)/(L - t) with L = log(8 x2)  => L > (m t + c)/(m-1)
    c = 4 * G / pi
    L = (m * t + c) / (m - 1)
    L = max(L, 5.0)                      # sigma_m condition
    print("  %2d  %2d  %.4f   %10.1f" % (k, m, t, 2.718281828459 ** L / 8))
