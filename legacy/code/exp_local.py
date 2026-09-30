"""
Experiment: local (single-prime) statistics versus the Galois / random
mapping predictions.

 (a) mean of P_c(p, k) over primes  vs  tau(k)  (generic c)
 (b) distribution of the number of rational k-cycles, k = 3, 4, 5, vs the
     cycle statistics of the wreath product Z/k wr S_{nu(k)/k}
 (c) cross-probe constraints:  rational 3-cycle of f_c  =>  f_{c+2} has a
     fixed point;  rational 4-cycle  =>  s^3 + (4c+3)s + 4 has a root mod p
 (d) number of periodic points N_c(p) / sqrt(p)  vs  sqrt(pi/2)
"""
import os
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import comb, factorial, pi, sqrt

from adt import SCRATCH, load_primes, mark

DATA = os.path.join(SCRATCH, "data", "primes_c{c}.txt")
GENERIC = [1, 2, 3, 4, 5, 7]
ALLC = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
loc = load_primes(DATA, ALLC)
primes = sorted({p for (p, c) in loc if p > 50})


def tau(k):
    return sum(1 for d in range(1, k + 1) if k % d == 0)


def nu(k):
    # number of points of exact period k for a degree-2 map (Mobius)
    mu = lambda n: 0 if any(n % (q * q) == 0 for q in range(2, n + 1)) else (-1) ** sum(
        1 for q in range(2, n + 1) if n % q == 0 and all(q % r for r in range(2, q)))
    return sum(mu(k // d) * 2 ** d for d in range(1, k + 1) if k % d == 0)


def wreath_law(k):
    """Law of #{blocks fixed pointwise} for a uniform element of Z/k wr S_r."""
    r = nu(k) // k
    # number of permutations of S_r with exactly f fixed points
    der = [1, 0]
    for n in range(2, r + 1):
        der.append((n - 1) * (der[-1] + der[-2]))
    law = Counter()
    for f in range(r + 1):
        pf = Fraction(comb(r, f) * der[r - f], factorial(r))
        for j in range(f + 1):
            law[j] += pf * comb(f, j) * Fraction(1, k) ** j * Fraction(k - 1, k) ** (f - j)
    return law


if __name__ == "__main__":
    print("(a) mean of P_c(p,k) over primes 50 < p < 1e5; prediction tau(k) for generic c")
    print("  k   tau  " + "".join(f"  c={c:<5}" for c in GENERIC))
    for k in range(1, 9):
        row = [sum(mark(loc[(p, c)], k) for p in primes) / len(primes) for c in GENERIC]
        print(f"  {k}   {tau(k)}   " + "".join(f"  {v:7.3f}" for v in row))

    print("\n(b) law of the number of rational k-cycles: data (c=1,2,3,4,5,7 pooled) vs wreath product")
    for k in (3, 4, 5):
        law = wreath_law(k)
        cnt = Counter(loc[(p, c)].get(k, 0) for p in primes for c in GENERIC)
        tot = sum(cnt.values())
        js = sorted(set(cnt) | set(law))
        print(f"  k={k} (r={nu(k)//k}):  " + "  ".join(
            f"j={j}: {cnt[j]/tot:.4f} vs {float(law.get(j, 0)):.4f}" for j in js if j <= 3))

    print("\n(c) cross-probe constraints (violations should be 0)")
    v3 = v4 = n3 = n4 = 0
    for p in primes:
        for c in [-6, -5, -4, -3, -1, 1, 2, 3, 4, 5]:
            if (p, c + 2) in loc and loc[(p, c)].get(3, 0) > 0:
                n3 += 1
                if mark(loc[(p, c + 2)], 1) == 0:
                    v3 += 1
            if loc[(p, c)].get(4, 0) > 0:
                n4 += 1
                roots = sum(1 for s in range(p) if (s * s * s + (4 * c + 3) * s + 4) % p == 0) if p < 3000 else None
                if roots is not None and roots < loc[(p, c)][4]:
                    v4 += 1
    print(f"  3-cycles: {n3} (p,c) with a rational 3-cycle; violations of 'f_(c+2) has a fixed point': {v3}")
    print(f"  4-cycles: {n4} (p,c) with a rational 4-cycle; violations of '#4-cycles <= #roots of cubic' (p<3000): {v4}")

    print("\n(d) mean N_c(p)/sqrt(p) over primes 1e4 < p < 1e5 (random-mapping prediction sqrt(pi/2) = 1.2533)")
    big = [p for p in primes if p > 10 ** 4]
    for c in ALLC:
        m = sum(mark(loc[(p, c)], 1) * 0 + sum(L * n for L, n in loc[(p, c)].items()) / sqrt(p) for p in big) / len(big)
        print(f"  c={c:3d}: {m:.4f}")
