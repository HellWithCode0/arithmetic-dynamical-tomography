"""
Experiment: exactly solvable probes on restricted semiprime classes.
Collisions of ADF_{C,oo}(pq) among semiprimes pq <= N for C = {0}, {-2}, {0,-2},
separately for Blum semiprimes (p = q = 3 mod 4) and for all odd semiprimes.
"""
import sys
from collections import defaultdict
import exp_solvable as S
from adt import zset_mul, zset_canon, primes_upto

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
P = [p for p in primes_upto(N // 3) if p > 2]
res = {}
for cls in ("blum", "all"):
    b0, b2, b02 = defaultdict(list), defaultdict(list), defaultdict(list)
    for i, p in enumerate(P):
        if p * p > N:
            break
        for q in P[i + 1:]:
            if p * q > N:
                break
            if cls == "blum" and (p % 4 != 3 or q % 4 != 3):
                continue
            a0, a2 = S.local(p)
            c0, c2 = S.local(q)
            f0, f2 = zset_canon(zset_mul(a0, c0)), zset_canon(zset_mul(a2, c2))
            b0[f0].append((p, q)); b2[f2].append((p, q)); b02[(f0, f2)].append((p, q))
    for name, b in (("{0}", b0), ("{-2}", b2), ("{0,-2}", b02)):
        coll = [v for v in b.values() if len(v) > 1]
        print(f"{cls:5} semiprimes <= {N}, C={name:7}: {sum(map(len, b.values())):7} numbers, "
              f"{len(coll)} colliding fibres; e.g. {coll[:3]}")
