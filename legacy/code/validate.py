"""
Validation checks for the computational core (Observation 2.5 and §4 of
paper/ADT_foundations.md):

  1. CRT synthesis: the periodic Z-set of f_c mod n computed by brute force
     equals the Burnside product of the local Z-sets of its prime factors,
     for all odd squarefree composites n <= 4000 and 300 random ones < 60000,
     and all 7 probes c = -3..3.
  2. Exact formulas for c = 0 and c = -2 (Theorems 4.1, 4.2) against brute
     force for every odd prime < 1e5, in two independent implementations.
  3. The collision 65 ~ 119 for C = {0, -2} (Example 5.5).

Needs the prime data from gen_data.sh (--small suffices).
"""
import os
import random

import exp_solvable as S
from adt import (SCRATCH, load_primes, spf_table, factor_sqfree, fingerprint, run_cycles,
                 zset_canon, primes_upto, per0_prime, per_m2_prime)

C = [-3, -2, -1, 0, 1, 2, 3]
loc = load_primes(os.path.join(SCRATCH, "data", "primes_c{c}.txt"), C)
spf = spf_table(200000)

comp = lambda n: (lambda ps: ps is not None and len(ps) > 1)(factor_sqfree(n, spf))
ns = [n for n in range(3, 4001, 2) if comp(n)]
random.seed(1)
ns += random.sample([n for n in range(4001, 60000, 2) if comp(n)], 300)
bad = 0
for i in range(0, len(ns), 200):
    for n, c, Z in run_cycles("list", ",".join(map(str, ns[i:i + 200])), None, C):
        bad += zset_canon(Z) != fingerprint(factor_sqfree(n, spf), [c], loc)[0]
print(f"1. CRT synthesis vs direct: {len(ns)} moduli x {len(C)} probes, mismatches: {bad}")

b0 = b2 = bf = 0
for p in primes_upto(100000)[1:]:
    b0 += zset_canon(per0_prime(p)) != zset_canon(loc[(p, 0)])
    b2 += zset_canon(per_m2_prime(p)) != zset_canon(loc[(p, -2)])
    if p < 20000:
        a, b = S.local(p)
        bf += (zset_canon(a), zset_canon(b)) != (zset_canon(loc[(p, 0)]), zset_canon(loc[(p, -2)]))
print(f"2. exact formulas, primes < 1e5: c=0 mismatches {b0}, c=-2 mismatches {b2}; "
      f"fast order-based implementation (p < 2e4) mismatches {bf}")

direct = {(n, c): zset_canon(Z) for n, c, Z in run_cycles("list", "65,119", None, [0, -2])}
print("3. Per_c(65) == Per_c(119) for c = 0, -2:",
      all(direct[(65, c)] == direct[(119, c)] for c in (0, -2)), {k: v for k, v in direct.items()})
