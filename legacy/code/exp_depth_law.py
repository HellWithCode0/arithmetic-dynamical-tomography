"""
Experiment: the depth-width law (Heuristics 8.1 and 8.3, Figures 2 and 3)
with all odd primes below 1e6 and nine generic probes.

 (1) For many generic probe sets S (|S| = 1..5), the least depth K making
     ADF_{S,K} injective on the primes in (X/2, X], X = 1e6 / 2^j.
     Heuristic 8.3 predicts log K ~ (1/|S|) log #primes; we report the
     least-squares slope of log K against log(#primes in the window).
 (2) Pair-collision fraction at depth K for primes in (5e5, 1e6]
     (Heuristic 8.1 predicts local slope -> -2|S|).
Needs build/data/full_c{c}.txt for c in GEN (see gen_data.sh).
Writes a table to stdout and results/depth_law.json (read by make_figures.py).
"""
import itertools
import json
import math
import os
from bisect import bisect_right
from collections import Counter
from multiprocessing import Pool

from adt import SCRATCH

DATA = os.path.join(SCRATCH, "data", "full_c{c}.txt")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "depth_law.json")
GEN = [1, 2, 3, 4, 5, 6, 7, -4, -5]
XS = [10 ** 6 // 2 ** j for j in range(9, -1, -1)]
KS = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
INF = 10 ** 9
LOC = None


def load(c):
    out = {}
    with open(DATA.format(c=c)) as fh:
        for line in fh:
            head, tail = line.split("|")
            p = int(head.split()[0])
            pairs = tuple(sorted((int(a), int(b)) for a, b in (t.split(":") for t in tail.split())))
            out[p] = (pairs, [L for L, _ in pairs])
    return out


def init():
    global LOC
    LOC = {c: load(c) for c in GEN}


def trunc(p, c, K):
    pairs, lens = LOC[c][p]
    return pairs[:bisect_right(lens, K)]


def least_depth(win, S):
    if max(Counter(tuple(LOC[c][p][0] for c in S) for p in win).values()) > 1:
        return INF
    lo, hi = 1, max(LOC[c][p][1][-1] for p in win for c in S)
    while lo < hi:
        mid = (lo + hi) // 2
        seen, ok = set(), True
        for p in win:
            key = tuple(trunc(p, c, mid) for c in S)
            if key in seen:
                ok = False
                break
            seen.add(key)
        if ok:
            hi = mid
        else:
            lo = mid + 1
    return lo


def task_depth(S):
    primes = sorted(LOC[S[0]])
    Ks, ns = [], []
    for X in XS:
        win = primes[bisect_right(primes, X // 2):bisect_right(primes, X)]
        Ks.append(least_depth(win, S))
        ns.append(len(win))
    pts = [(math.log(n), math.log(K)) for n, K in zip(ns, Ks) if K < INF]
    slope = None
    if len(pts) >= 3:
        mx = sum(x for x, _ in pts) / len(pts)
        my = sum(y for _, y in pts) / len(pts)
        slope = sum((x - mx) * (y - my) for x, y in pts) / sum((x - mx) ** 2 for x, _ in pts)
    return list(S), Ks, ns, slope


def task_coll(S):
    primes = sorted(LOC[S[0]])
    win = primes[bisect_right(primes, 500000):]
    npairs = len(win) * (len(win) - 1) / 2
    res = []
    for K in KS:
        cnt = Counter(tuple(trunc(p, c, K) for c in S) for p in win)
        res.append(sum(m * (m - 1) / 2 for m in cnt.values()) / npairs)
    return list(S), len(win), res


def pick(k, m):
    combos = list(itertools.combinations(GEN, k))
    if len(combos) <= m:
        return combos
    step = len(combos) / m
    return [combos[int(i * step)] for i in range(m)]


if __name__ == "__main__":
    depth_sets = pick(1, 9) + pick(2, 36) + pick(3, 30) + pick(4, 30) + pick(5, 20)
    coll_sets = [(2,), (2, 5), (2, 5, -4)] + [s for s in pick(1, 9) + pick(2, 6) + pick(3, 4)
                                              if s not in [(2,), (2, 5), (2, 5, -4)]]
    with Pool(min(8, os.cpu_count() or 1), initializer=init) as pool:
        depth = pool.map(task_depth, depth_sets, chunksize=1)
        coll = pool.map(task_coll, coll_sets, chunksize=1)
    print("(1) least depth making ADF_{S,K} injective on the primes in (X/2, X], X =", XS)
    print("    #primes per window:", depth[0][2])
    summary = {}
    for S, Ks, ns, slope in depth:
        summary.setdefault(len(S), []).append(slope)
        print(f"  {str(S):22} K = {' '.join(str(K) if K < INF else 'inf' for K in Ks):60} "
              f"slope {'n/a' if slope is None else f'{slope:.3f}'}")
    for s in sorted(summary):
        v = sorted(x for x in summary[s] if x is not None)
        if v:
            q = lambda f: v[min(len(v) - 1, int(f * len(v)))]
            print(f"  |S|={s}: {len(v)}/{len(summary[s])} sets with >= 3 finite points; median LS slope "
                  f"{q(0.5):.3f} (IQR {q(0.25):.3f}-{q(0.75):.3f}); heuristic 1/|S| = {1/s:.3f}")
        else:
            print(f"  |S|={s}: full-depth collisions in almost every window (Heuristic 8.2)")
    print("(2) pair-collision fraction at depth K, primes in (5e5, 1e6]; local slope in brackets")
    for S, n, res in coll:
        parts = []
        for i, (K, f) in enumerate(zip(KS, res)):
            if i and f > 0 and res[i - 1] > 0:
                parts.append(f"K={K}:{f:.2e}({(math.log(f) - math.log(res[i-1])) / (math.log(K) - math.log(KS[i-1])):+.2f})")
            else:
                parts.append(f"K={K}:{f:.2e}")
        print(f"  S={S} (n={n}, heuristic slope -> {-2 * len(S)}): " + " ".join(parts))
    with open(OUT, "w") as fh:
        json.dump({"XS": XS, "KS": KS, "INF": INF,
                   "depth": [{"S": S, "K": Ks, "n": ns, "slope": sl} for S, Ks, ns, sl in depth],
                   "coll": [{"S": S, "n": n, "frac": res} for S, n, res in coll]}, fh)
