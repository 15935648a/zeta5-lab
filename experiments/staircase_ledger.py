"""
M3 opening: the staircase feasibility ledger for zeta(3) in CDT's framework.

CDT bound (univalent form, Thm 2.7.1): m <= lambda / (lambda - tau(b)),
tau(b) = (1/m^2) sum_i (2i-1) sigma_i, types sorted ascending.
Contradiction (irrationality) achievable iff  m > lambda/(lambda - tau),
i.e. lambda > tau * m/(m-1).

Inventory: ONE conditional function P of type t=3 ([1..n]^3 under zeta(3)
rational), plus N companion functions of type 0 (integral coefficients)
that must ALSO be overconvergent at the covered point x1 (the
qualification that disqualifies Apery's A).

Question: for each N, what lambda is required, vs the available ceiling
log 16 = 2.7726 (Caratheodory) and CDT's nested-map ladder
16^(1 - 2^-(n+1)) (n square roots, currently built for Sigma^0={0};
adaptation to Sigma^0={delta} with delta=(sqrt2-1)^8 is the open
conformal problem).
"""
from math import log, sqrt

def tau_staircase(types):
    m = len(types)
    s = sorted(types)
    return sum((2 * i + 1) * s[i] for i in range(m)) / m**2

t = 3.0
print("delta = (sqrt2-1)^8 = %.6f   (covered point, scaled)" % ((sqrt(2) - 1) ** 8))
print("ceiling log16 = %.4f ; e^2.5 = %.3f" % (log(16), 2.718281828 ** 2.5))
print("\n N(type-0 companions)  m   tau(b)   need lambda >   phi'(0) >   verdict vs 16")
for N in range(0, 7):
    types = [0.0] * N + [t]
    m = len(types)
    tau = tau_staircase(types)
    if m == 1:
        print("  %d                    %d   %.4f     (m=1: no bound game)" % (N, m, tau))
        continue
    lam_req = tau * m / (m - 1)
    phi_req = 2.718281828459045 ** lam_req
    verdict = "FEASIBLE" if phi_req < 16 else "dead"
    print("  %d                    %d   %.4f      %.4f        %8.3f    %s"
          % (N, m, tau, lam_req, phi_req, verdict))

print("\n with derivatives of P also in the tuple (each type 3):")
for N in range(0, 7):
    for dP in (1, 2, 3):
        types = [0.0] * N + [t] * dP
        m = len(types)
        if m == 1:
            continue
        tau = tau_staircase(types)
        lam_req = tau * m / (m - 1)
        if 2.718281828 ** lam_req < 16 and dP > 1:
            print("  N=%d, %d copies of type-3: tau=%.3f, need phi'(0) > %.2f  FEASIBLE"
                  % (N, dP, tau, 2.718281828 ** lam_req))

print("\n CDT nested-map ladder (Sigma^0={0} version): phi'_n(0) = 16^(1-2^-(n+1))")
for n in range(0, 6):
    v = 16 ** (1 - 2 ** (-(n + 1)))
    print("  n=%d: %.4f   (log = %.4f)" % (n, v, log(v)))

print("\n zeta(5) comparison: type t=5, N type-0 companions:")
for N in (2, 4, 8, 16, 64):
    types = [0.0] * N + [5.0]
    m = len(types)
    tau = tau_staircase(types)
    lam_req = tau * m / (m - 1)
    print("  N=%2d: tau=%.4f, need phi'(0) > %.2f  %s"
          % (N, tau, 2.718281828 ** lam_req,
             "feasible" if 2.718281828 ** lam_req < 16 else "dead"))
