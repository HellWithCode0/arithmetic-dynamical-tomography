"""
Experiment: how much depth is needed?

 (a) For probe sets S of generic constants, the least K such that the
     depth-K fingerprint ADF_{S,K} is injective on the primes in (X/2, X].
     The Poisson-Chebotarev heuristic predicts  K ~ pi(X)^(1/|S|)  (up to
     constants): the "tomographic uncertainty principle" |S| log K ~ log X.
 (b) Theorem target B for the paper's C = {-3..3}: the least K(X) such that
     ADF_{C,K} is injective on odd semiprimes pq <= X.
"""
import os
import sys
from bisect import bisect_right
from collections import defaultdict
from math import log

from adt import SCRATCH, load_primes, truncate, zset_mul, zset_canon

SCR = SCRATCH
INF = 10 ** 9


def least_depth(items, Kmax=INF):
    """items: list of tuples of Z-sets (dicts), one per object.  Return the
    least K making the depth-K fingerprint injective (INF if even the full
    Z-sets collide)."""
    full = defaultdict(int)
    for zs in items:
        full[tuple(zset_canon(z) for z in zs)] += 1
    if max(full.values()) > 1:
        return INF
    lo, hi = 1, max(max(z) for zs in items for z in zs)
    while lo < hi:
        mid = (lo + hi) // 2
        seen = set()
        ok = True
        for zs in items:
            key = tuple(truncate(z, mid) for z in zs)
            if key in seen:
                ok = False
                break
            seen.add(key)
        if ok:
            hi = mid
        else:
            lo = mid + 1
    return lo


if __name__ == "__main__":
    part = sys.argv[1] if len(sys.argv) > 1 else "ab"
    if "a" in part:
        GEN = [1, 2, 3, 4, 5, 7, -4, -5]
        loc = load_primes(os.path.join(SCR, "data", "primes_c{c}.txt"), GEN)
        primes = sorted({p for (p, c) in loc})
        sets = {1: [[2], [5], [-4]], 2: [[2, 5], [1, 4], [3, -5]],
                3: [[2, 5, -4], [1, 4, -5], [3, 7, -4]], 4: [[1, 2, 5, -4], [3, 4, 7, -5]]}
        Xs = [1250, 2500, 5000, 10000, 20000, 40000, 80000]
        print("(a) least injective depth on primes in (X/2, X]")
        print("  |S|  S            " + "".join(f"{X:>8}" for X in Xs) + "   slope d(log K)/d(log #primes)")
        for s, Ss in sets.items():
            for S in Ss:
                Ks, ns = [], []
                for X in Xs:
                    win = primes[bisect_right(primes, X // 2):bisect_right(primes, X)]
                    Ks.append(least_depth([tuple(loc[(p, c)] for c in S) for p in win]))
                    ns.append(len(win))
                fin = [(log(n), log(K)) for n, K in zip(ns, Ks) if K < INF]
                slope = ((fin[-1][1] - fin[0][1]) / (fin[-1][0] - fin[0][0])) if len(fin) > 1 else float("nan")
                print(f"  {s}    {str(S):12} " + "".join(f"{(K if K < INF else 'inf'):>8}" for K in Ks)
                      + f"   {slope:.2f}  (heuristic 1/|S| = {1/s:.2f})")

    if "b" in part:
        C = [-3, -2, -1, 0, 1, 2, 3]
        loc = load_primes(os.path.join(SCR, "data", "full_c{c}.txt"), C)
        primes = sorted({p for (p, c) in loc})
        print("\n(b) Theorem target B for C = {-3..3}: least K making ADF_{C,K} injective on odd semiprimes pq <= X")
        for X in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]:
            items = []
            for i, p in enumerate(primes):
                if p * p > X:
                    break
                for q in primes[i + 1:]:
                    if p * q > X:
                        break
                    items.append(tuple(zset_mul(loc[(p, c)], loc[(q, c)]) for c in C))
            print(f"  X = {X:>8}: {len(items):>7} semiprimes, least depth K(X) = {least_depth(items)}")
