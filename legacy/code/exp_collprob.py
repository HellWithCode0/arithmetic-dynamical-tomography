"""
Experiment: per-pair collision probability of depth-K fingerprints of
primes.  Poisson-Chebotarev prediction: for |S| independent generic probes,
Pr[ADF_{S,K}(p) = ADF_{S,K}(q)] ~ kappa_S * K^(-2|S|), i.e. a log-log slope
of -2|S| in K (independent of p, q).  A polylog depth could only suffice if
this decayed faster than every power of K.
"""
import os
from bisect import bisect_right
from collections import Counter
from math import log

from adt import SCRATCH, load_primes, truncate

loc = load_primes(os.path.join(SCRATCH, "data", "primes_c{c}.txt"), [1, 2, 3, 4, 5, -4, -5])
primes = sorted({p for (p, c) in loc})
win = primes[bisect_right(primes, 40000):]           # 40000 < p < 1e5
npairs = len(win) * (len(win) - 1) // 2
Ks = [2, 4, 8, 16, 32, 64]
print(f"{len(win)} primes in (40000, 1e5); {npairs} pairs")
for S in ([2], [5], [2, 5], [1, 4], [2, 5, -4], [1, 4, -5]):
    probs = []
    for K in Ks:
        cnt = Counter(tuple(truncate(loc[(p, c)], K) for c in S) for p in win)
        probs.append(sum(v * (v - 1) // 2 for v in cnt.values()) / npairs)
    sl = [(log(probs[i + 1]) - log(probs[i])) / log(2) if probs[i + 1] > 0 else float('nan') for i in range(len(Ks) - 1)]
    print(f"S={str(S):12} " + " ".join(f"K={K}:{q:.2e}" for K, q in zip(Ks, probs)))
    print(f"{'':15}local slopes {', '.join(f'{s:.2f}' for s in sl)}   (prediction -> {-2*len(S)})")
