"""
Experiment: search for collisions of the full (infinite-depth) fingerprint
ADF_{C,oo}(n) for the paper's probe set C = {-3,...,3} over odd squarefree n.

A collision n ~ m is called primitive if gcd(n, m) = 1 (every collision
n*r ~ m*r obtained by multiplying a collision by a common coprime factor
is non-primitive).
"""
import os
import sys
from collections import defaultdict
from math import gcd

from adt import SCRATCH, load_primes, spf_table, factor_sqfree, fingerprint

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
C = [int(c) for c in sys.argv[2].split(",")] if len(sys.argv) > 2 else [-3, -2, -1, 0, 1, 2, 3]
DATA = os.path.join(SCRATCH, "data", os.environ.get("ADT_PRIMEFILE", "full_c{c}.txt"))
loc = load_primes(DATA, C)
spf = spf_table(N)

# key on the 64-bit hash of the exact fingerprint; candidate collisions are
# re-verified on the exact (unhashed) fingerprints below
hashes = {}
cand = defaultdict(set)
count = skipped = 0
for n in range(3, N + 1, 2):
    ps = factor_sqfree(n, spf)
    if ps is None:
        continue
    if any((p, C[0]) not in loc for p in ps):
        skipped += 1
        continue
    count += 1
    h = hash(fingerprint(ps, C, loc))
    if h in hashes:
        cand[h].update((hashes[h], n))
    else:
        hashes[h] = n
buckets = defaultdict(list)
for h, ns in cand.items():
    for n in sorted(ns):
        buckets[fingerprint(factor_sqfree(n, spf), C, loc)].append(n)
print(f"odd squarefree n <= {N} skipped for lack of local data (a prime factor > data range): {skipped}")
coll = [v for v in buckets.values() if len(v) > 1]
pairs = [(a, b) for v in coll for i, a in enumerate(v) for b in v[i + 1:]]
prim = [(a, b) for a, b in pairs if gcd(a, b) == 1]
print(f"C = {C}: odd squarefree n <= {N} examined: {count}, distinct fingerprints: {len(hashes)}")
print(f"colliding pairs: {len(pairs)}, primitive (coprime) pairs: {len(prim)}")
for a, b in sorted(prim)[:25]:
    print(f"  {a} = {factor_sqfree(a, spf)}  ~  {b} = {factor_sqfree(b, spf)}")
