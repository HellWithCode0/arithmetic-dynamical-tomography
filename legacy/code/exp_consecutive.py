"""
Experiment: consecutive constants (the proposal's C = {-3..3}) versus a
7-element set of generic, non-consecutive constants.  Compared on
 (i) least injective depth on primes in (X/2, X],
 (ii) collision probability of depth-K fingerprints of primes,
 (iii) least injective depth on odd semiprimes pq <= 1e5 (Target B).
"""
import os
from bisect import bisect_right
from collections import Counter

from adt import SCRATCH, load_primes, truncate, zset_mul
from exp_depth import least_depth

SETS = {"consecutive {-3..3}": [-3, -2, -1, 0, 1, 2, 3], "spread {-5,-4,1,2,4,5,7}": [-5, -4, 1, 2, 4, 5, 7]}
loc = load_primes(os.path.join(SCRATCH, "data", "primes_c{c}.txt"), sorted({c for S in SETS.values() for c in S}))
primes = sorted({p for (p, c) in loc})
for name, C in SETS.items():
    row = []
    for X in (5000, 20000, 100000):
        w = primes[bisect_right(primes, X // 2):bisect_right(primes, X)]
        row.append(least_depth([tuple(loc[(p, c)] for c in C) for p in w]))
    w = primes[bisect_right(primes, 40000):]
    npairs = len(w) * (len(w) - 1) // 2
    pr = []
    for K in (2, 3, 4, 6):
        cnt = Counter(tuple(truncate(loc[(p, c)], K) for c in C) for p in w)
        pr.append(sum(v * (v - 1) // 2 for v in cnt.values()) / npairs)
    items = [tuple(zset_mul(loc[(p, c)], loc[(q, c)]) for c in C)
             for i, p in enumerate(primes) if p * p <= 10 ** 5 for q in primes[i + 1:] if p * q <= 10 ** 5]
    print(f"{name:26} least depth primes (X=5e3,2e4,1e5): {row};  Pr[collide] at K=2,3,4,6: "
          + ", ".join(f"{x:.1e}" for x in pr) + f";  semiprimes<=1e5: K = {least_depth(items)}")
