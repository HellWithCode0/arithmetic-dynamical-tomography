"""
Experiment: collisions among PRIMES at infinite depth.

For a single generic probe the random-mapping heuristic predicts that the
number of pairs p < q <= X with Per_c(p) = Per_c(q) grows like X^{1/2}/log^2 X
(so infinitely many), while for two independent probes the expected number
of prime collisions is finite.
"""
import os
import sys
from collections import defaultdict
from itertools import combinations

from adt import SCRATCH, load_primes, zset_canon

DATA = os.path.join(SCRATCH, "data", "primes_c{c}.txt")
PROBES = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
loc = load_primes(DATA, PROBES)
primes = sorted({p for (p, c) in loc})


def collisions(C, X):
    buckets = defaultdict(list)
    for p in primes:
        if p > X:
            break
        buckets[tuple(zset_canon(loc[(p, c)]) for c in C)].append(p)
    pairs = sum(len(v) * (len(v) - 1) // 2 for v in buckets.values())
    ex = [v for v in buckets.values() if len(v) > 1]
    return pairs, ex


if __name__ == "__main__":
    Xs = [1000, 3000, 10000, 30000, 100000]
    print("single probe: number of colliding prime pairs p<q<=X")
    print("   c  " + "".join(f"{X:>9}" for X in Xs))
    for c in PROBES:
        row = [collisions([c], X)[0] for X in Xs]
        print(f"{c:4d}  " + "".join(f"{v:>9}" for v in row))

    print("\nlargest colliding primes, generic single probes (X = 1e5):")
    for c in [1, 2, 3, 5]:
        _, ex = collisions([c], 10 ** 5)
        big = sorted(ex, key=lambda v: -min(v))[:3]
        for v in big:
            print(f"  c={c}: {v}  Per = {zset_canon(loc[(v[0], c)])}")

    print("\npairs of probes: colliding prime pairs up to 1e5")
    gen = [-6, -5, -4, 1, 2, 3, 4, 5, 6, 7]
    tot = 0
    for c1, c2 in combinations(gen, 2):
        pairs, ex = collisions([c1, c2], 10 ** 5)
        tot += pairs
        if pairs:
            print(f"  C={{{c1},{c2}}}: {pairs} pair(s): {ex[:4]}")
    print(f"  total over {len(list(combinations(gen, 2)))} generic probe pairs: {tot}")

    C7 = [-3, -2, -1, 0, 1, 2, 3]
    pairs, ex = collisions(C7, 10 ** 5)
    print(f"\npaper's C={{-3..3}}: colliding prime pairs up to 1e5: {pairs} {ex[:5]}")
