"""
adt.py -- computational core for Arithmetic Dynamical Tomography.

Conventions
-----------
f_c(x) = x^2 + c acting on Z/nZ.
Per_c(n) = the periodic part of the functional graph, viewed as a finite
Z-set (a permutation up to isomorphism).  It is stored as a dict
{cycle length: number of cycles}.  The Burnside-ring product of Z-sets is

    [Z/a] * [Z/b] = gcd(a, b) [Z/lcm(a, b)],

and CRT gives Per_c(ab) = Per_c(a) * Per_c(b) for coprime a, b.

The mark of a Z-set X at k is |Fix(sigma^k)| = sum_{d | k} d * X[d];
for X = Per_c(n) this is exactly the periodic-point count P_c(n, k).
"""
from __future__ import annotations

import os
import subprocess
from collections import defaultdict
from functools import lru_cache
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.environ.get("ADT_SCRATCH", os.path.join(HERE, "..", "build"))
CYCLES_BIN = os.path.join(SCRATCH, "cycles")


# ---------------------------------------------------------------- Z-sets
def zset_mul(A: dict, B: dict) -> dict:
    out = defaultdict(int)
    for a, ca in A.items():
        for b, cb in B.items():
            g = gcd(a, b)
            out[a // g * b] += ca * cb * g
    return dict(out)


def zset_canon(A: dict) -> tuple:
    return tuple(sorted((k, v) for k, v in A.items() if v))


def mark(A: dict, k: int) -> int:
    """P(k) = number of points fixed by sigma^k."""
    return sum(d * c for d, c in A.items() if k % d == 0)


def truncate(A: dict, K: int) -> tuple:
    """Depth-K information: the cycles of length <= K (equivalent to the
    marks P(1..K) by Mobius inversion)."""
    return tuple(sorted((d, c) for d, c in A.items() if d <= K and c))


def zset_from_marks(P: list[int]) -> dict:
    """Inverse of k -> mark(A, k) for k = 1..len(P) (Mobius inversion)."""
    K = len(P)
    exact = [0] * (K + 1)
    for d in range(1, K + 1):
        exact[d] = P[d - 1] - sum(exact[e] for e in range(1, d) if d % e == 0)
    return {d: exact[d] // d for d in range(1, K + 1) if exact[d]}


ONE = {1: 1}


# ---------------------------------------------------------- local data
def parse_line(line: str):
    head, tail = line.split("|")
    n, c, nper = head.split()
    Z = {}
    for tok in tail.split():
        L, m = tok.split(":")
        Z[int(L)] = int(m)
    return int(n), int(c), Z


def run_cycles(mode: str, a, b, cs) -> list:
    args = [CYCLES_BIN, mode] + ([a] if mode == "list" else [str(a), str(b)]) + [str(c) for c in cs]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    return [parse_line(l) for l in out.splitlines() if l.strip()]


def load_primes(path_pattern: str, cs) -> dict:
    """Load {(p, c): Z-set} from files produced by `cycles primes`."""
    data = {}
    for c in cs:
        with open(path_pattern.format(c=c)) as fh:
            for line in fh:
                n, cc, Z = parse_line(line)
                data[(n, cc)] = Z
    return data


# --------------------------------------------------------- arithmetic
def primes_upto(N: int) -> list[int]:
    s = bytearray([1]) * (N + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(s[i * i :: i]))
    return [i for i in range(N + 1) if s[i]]


def spf_table(N: int) -> list[int]:
    spf = list(range(N + 1))
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factor_sqfree(n: int, spf) -> list[int] | None:
    ps = []
    while n > 1:
        p = spf[n]
        n //= p
        if n % p == 0:
            return None
        ps.append(p)
    return ps


def odd_part(m: int) -> int:
    while m % 2 == 0:
        m //= 2
    return m


def mult_order(a: int, m: int) -> int:
    if m == 1:
        return 1
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


# ------------------------------------------------------- fingerprints
def fingerprint(ps, C, local, K=None):
    """ADF_{C,K}(prod ps) synthesised multiplicatively from local data.
    K=None means infinite depth (the full periodic Z-set)."""
    fp = []
    for c in C:
        Z = ONE
        for p in ps:
            Z = zset_mul(Z, local[(p, c)])
        fp.append(zset_canon(Z) if K is None else truncate(Z, K))
    return tuple(fp)


# ------------------------------------- exactly solvable probes c = 0, -2
@lru_cache(maxsize=None)
def D_zset(m: int) -> tuple:
    """Z-set of a -> 2a on Z/m (m odd): for each d | m, phi(d)/ord_d(2)
    cycles of length ord_d(2)."""
    out = defaultdict(int)
    for d in range(1, m + 1):
        if m % d == 0:
            o = mult_order(2, d)
            out[o] += euler_phi(d) // o
    return zset_canon(out)


@lru_cache(maxsize=None)
def E_zset(m: int) -> tuple:
    """Z-set of a -> 2a on (Z/m)/{+-1} (m odd)."""
    seen, out = set(), defaultdict(int)
    reps = {}
    for a in range(m):
        reps[a] = min(a, (-a) % m)
    for a in range(m):
        r = reps[a]
        if r in seen:
            continue
        # walk until a repeat; every class is periodic (doubling is a bijection)
        cyc, x = [], r
        while True:
            cyc.append(x)
            seen.add(x)
            x = reps[(2 * x) % m]
            if x == r:
                break
        out[len(cyc)] += 1
    return zset_canon(out)


@lru_cache(maxsize=None)
def euler_phi(n: int) -> int:
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def per0_prime(p: int) -> dict:
    """Per_0(p) = [1] + D(odd part of p-1)   (Theorem 4.1)."""
    Z = defaultdict(int, dict(D_zset(odd_part(p - 1))))
    Z[1] += 1
    return dict(Z)


def per_m2_prime(p: int) -> dict:
    """Per_{-2}(p) = E(odd(p-1)) v E(odd(p+1)) (wedge at the fixed point)."""
    Z = defaultdict(int, dict(E_zset(odd_part(p - 1))))
    for L, c in E_zset(odd_part(p + 1)):
        Z[L] += c
    Z[1] -= 1
    return dict(Z)
