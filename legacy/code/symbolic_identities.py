"""
Symbolic identities used in Proposition 3.3 of paper/ADT_foundations.md.

  (1) Phi_2(x, c) = Phi_1(-x, c + 1),  discriminants 1 - 4c and -3 - 4c
  (2) disc_x Phi_3 = -(4c+7)^3 (16c^2+4c+7)^2
      Res_x(Phi_3, s - (x + f(x) + f^2(x))) = (s^2 + s + c + 2)^3
  (3) Res_x(Phi_4, s - (x + f(x) + f^2(x) + f^3(x))) = (s^3 + (4c+3)s + 4)^4
  (4) factorisation pattern of Phi_3, Phi_4 for c = -3..3

Requires sympy.  The degree-4 resultant takes a few minutes.
"""
import sympy as sp

x, c, s = sp.symbols("x c s")


def orbit(k):
    pts = [x]
    for _ in range(k - 1):
        pts.append(sp.expand(pts[-1] ** 2 + c))
    return pts


def f_iter(k):
    t = x
    for _ in range(k):
        t = sp.expand(t ** 2 + c)
    return t


Phi1 = x ** 2 - x + c
Phi2 = sp.cancel((f_iter(2) - x) / Phi1)
Phi3 = sp.expand(sp.cancel((f_iter(3) - x) / Phi1))
Phi4 = sp.expand(sp.cancel((f_iter(4) - x) / (Phi1 * Phi2)))

print("(1) Phi_2 =", sp.factor(Phi2), "; Phi_2(x,c) - Phi_1(-x,c+1) =",
      sp.expand(Phi2 - Phi1.subs({x: -x, c: c + 1}, simultaneous=True)))
print("    disc Phi_1 =", sp.discriminant(Phi1, x), "; disc Phi_2 =", sp.factor(sp.discriminant(Phi2, x)))
print("(2) disc_x Phi_3 =", sp.factor(sp.discriminant(Phi3, x)))
print("    Res(Phi_3, s - cycle sum) =", sp.factor(sp.resultant(Phi3, s - sum(orbit(3)), x)))
for cv in range(-3, 4):
    d3 = sorted(sp.degree(q, x) for q, _ in sp.factor_list(Phi3.subs(c, cv))[1])
    d4 = sorted(sp.degree(q, x) for q, _ in sp.factor_list(Phi4.subs(c, cv))[1])
    print(f"(4) c = {cv:2d}: degrees of irreducible factors of Phi_3: {d3}, of Phi_4: {d4}")
print("(3) Res(Phi_4, s - cycle sum) =", sp.factor(sp.resultant(Phi4, s - sum(orbit(4)), x)))
