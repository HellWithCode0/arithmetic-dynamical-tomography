"""
Experiment: the exactly solvable probe pair C = {0, -2}.

Local data (Theorems 4.1 and 4.2 of the paper):
    Per_0(p)    = [1] + D(m_-),            m_- = odd part of p - 1
    Per_{-2}(p) = E(m_-) v E(m_+),         m_+ = odd part of p + 1
with
    D(m) = sum_{d | m} phi(d)/ord_d(2) [Z/ord_d(2)]
    E(m) = [1] + sum_{d | m, d > 1} phi(d)/(2 ord'_d) [Z/ord'_d],
    ord'_d = least k with 2^k = +-1 (mod d).

We search for collisions ADF(n) = ADF(m) among odd squarefree n <= N.
"""
import sys
from collections import defaultdict
from functools import lru_cache

from adt import zset_mul, zset_canon, spf_table, factor_sqfree, odd_part, ONE

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
SPF = spf_table(2 * N + 2)


def factorize(m):
    f = {}
    while m > 1:
        p = SPF[m]
        f[p] = f.get(p, 0) + 1
        m //= p
    return f


def divisors_with_phi(m):
    out = [(1, 1)]
    for p, e in factorize(m).items():
        new = []
        for d, ph in out:
            pk, phk = 1, 1
            for k in range(e + 1):
                new.append((d * pk, ph * phk))
                phk = (p - 1) * pk if k == 0 else phk * p
                pk *= p
        out = new
    return out


@lru_cache(maxsize=None)
def order2(d, pm=False):
    """ord_d(2), or with pm=True the least k with 2^k = +-1 mod d."""
    if d == 1:
        return 1
    phi = 1
    for p, e in factorize(d).items():
        phi *= (p - 1) * p ** (e - 1)
    o = phi
    for q in factorize(phi):
        while o % q == 0 and pow(2, o // q, d) == 1:
            o //= q
    if not pm:
        return o
    return o // 2 if o % 2 == 0 and pow(2, o // 2, d) == d - 1 else o


@lru_cache(maxsize=None)
def D(m):
    Z = defaultdict(int)
    for d, ph in divisors_with_phi(m):
        o = order2(d)
        Z[o] += ph // o
    return zset_canon(Z)


@lru_cache(maxsize=None)
def E(m):
    Z = defaultdict(int)
    Z[1] += 1
    for d, ph in divisors_with_phi(m):
        if d > 1:
            o = order2(d, True)
            Z[o] += ph // (2 * o)
    return zset_canon(Z)


@lru_cache(maxsize=None)
def local(p):
    mm, mp = odd_part(p - 1), odd_part(p + 1)
    Z0 = defaultdict(int, dict(D(mm)))
    Z0[1] += 1
    Z2 = defaultdict(int, dict(E(mm)))
    for L, c in E(mp):
        Z2[L] += c
    Z2[1] -= 1
    return dict(Z0), dict(Z2)


def fp(ps):
    A, B = ONE, ONE
    for p in ps:
        a, b = local(p)
        A, B = zset_mul(A, a), zset_mul(B, b)
    return zset_canon(A), zset_canon(B)


if __name__ == "__main__":
    seen = {}
    coll = []
    only0 = defaultdict(list)
    count = 0
    for n in range(3, N + 1, 2):
        ps = factor_sqfree(n, SPF)
        if ps is None:
            continue
        count += 1
        f = fp(ps)
        only0[hash(f[0])].append(n)
        if f in seen:
            coll.append((seen[f], n))
        else:
            seen[f] = n
    print(f"odd squarefree n <= {N}: {count}")
    print(f"collisions for C={{0,-2}}: {len(coll)}")
    for a, b in coll[:40]:
        print("  ", a, factor_sqfree(a, SPF), "~", b, factor_sqfree(b, SPF))
    sizes = sorted((len(v) for v in only0.values()), reverse=True)
    print(f"C={{0}} alone: {len(only0)} distinct fingerprints; largest fibres {sizes[:8]}")
