from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from math import gcd


def primes_below(limit: int) -> list[int]:
    if limit <= 2:
        return []
    sieve = bytearray(b"\x01") * limit
    sieve[:2] = b"\x00\x00"
    for p in range(2, int((limit - 1) ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = b"\x00" * len(sieve[p * p :: p])
    return [p for p in range(limit) if sieve[p]]


def odd_part(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    return n


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors_phi(n: int):
    values = [(1, 1)]
    for p, e in factor(n).items():
        nxt = []
        for d, ph in values:
            pk = 1
            for k in range(e + 1):
                phk = 1 if k == 0 else (p - 1) * p ** (k - 1)
                nxt.append((d * pk, ph * phk))
                pk *= p
        values = nxt
    return values


@lru_cache(maxsize=None)
def order2(d: int, modulo_sign: bool = False) -> int:
    if d == 1:
        return 1
    phi = 1
    for p, e in factor(d).items():
        phi *= (p - 1) * p ** (e - 1)
    order = phi
    for q in factor(phi):
        while order % q == 0 and pow(2, order // q, d) == 1:
            order //= q
    if modulo_sign and order % 2 == 0 and pow(2, order // 2, d) == d - 1:
        return order // 2
    return order


def canon(z: dict[int, int]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((length, count) for length, count in z.items() if count))


def zset_mul(a, b):
    out = defaultdict(int)
    for la, ca in a.items():
        for lb, cb in b.items():
            g = gcd(la, lb)
            out[la // g * lb] += ca * cb * g
    return dict(out)


@lru_cache(maxsize=None)
def local_solvable(p: int):
    minus, plus = odd_part(p - 1), odd_part(p + 1)
    z0 = defaultdict(int)
    for d, phi in divisors_phi(minus):
        order = order2(d)
        z0[order] += phi // order
    z0[1] += 1
    zm2 = defaultdict(int, {1: 1})
    for m in (minus, plus):
        for d, phi in divisors_phi(m):
            if d > 1:
                order = order2(d, True)
                zm2[order] += phi // (2 * order)
    return dict(z0), dict(zm2)


def solvable_fingerprint(prime_factors):
    a = b = {1: 1}
    for p in prime_factors:
        pa, pb = local_solvable(p)
        a, b = zset_mul(a, pa), zset_mul(b, pb)
    return canon(a), canon(b)


def periodic_cycles(n: int, c: int) -> dict[int, int]:
    """Cycle inventory of x -> x^2+c modulo n, by direct enumeration."""
    state = bytearray(n)
    cycles = defaultdict(int)
    for start in range(n):
        if state[start]:
            continue
        path, pos, x = [], {}, start
        while not state[x] and x not in pos:
            pos[x] = len(path)
            path.append(x)
            x = (x * x + c) % n
        if x in pos:
            cycles[len(path) - pos[x]] += 1
        for y in path:
            state[y] = 1
    return dict(cycles)


def marks(z: dict[int, int], depth: int) -> tuple[int, ...]:
    return tuple(sum(d * count for d, count in z.items() if k % d == 0)
                 for k in range(1, depth + 1))
